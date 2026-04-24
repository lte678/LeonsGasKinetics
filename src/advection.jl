# Leon's Gas Kinetics
# Copyright (C) 2026  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

using JET
using Atomix
import AcceleratedKernels as AK

"""
Performs the advection step.
"""
function advect!(sim::SimulationState, mesh::Mesh, config, dt)
    particles = sim.particles

    if config.vrbgk.enabled
        fill!(particles.features.last_collided_side, 0)
        for side in mesh.bc_sides
            @atomic side.vrbgk_incident_sum = 0.0
            @atomic side.vrbgk_incident_count = 0
        end
    end

    # Just move each particle and handle boundary collisions as they occur. No more, no less.
    # The Val objects allow for specialization on the boolean. The assert is surprisingly expensive.
    advect_move!(particles, mesh, config, dt, Val(config.degrees_of_freedom), Val(config.asserts))

    # Remove particles absorbed by open boundaries (marked cell = 0 by take_advection_step).
    compact_deleted_particles!(sim.particles)

    # Insert new particles entering through open boundaries.
    insert_open_boundary_particles!(sim, mesh, config, dt)

    # Recount after parallel advection to avoid write races on cell_part_count.
    fill!(sim.cells.part_count, 0)
    for p in sim.particles
        sim.cells.part_count[p.cell] += 1
    end

    # Apply VRBGK fixes
    if config.vrbgk.enabled
        vrbgk_advection_finalize(particles, mesh, config)
    end

    if config.asserts
        assert_particles_in_mesh(sim.particles, mesh)
        assert_cell_part_count(sim.particles, sim.cells.part_count)
    end
end


function advect_move!(particles, mesh, config, dt, dofs, asserts)
    AK.foreachindex(particles, max_tasks=Threads.nthreads()) do i
        time_remaining = dt
        pos = particles.pos[i]
        vel = particles.vel[i]
        while time_remaining > 0.0
            # By encoding this flag as a type, we can check it when advect_move! is compiled
            if asserts isa Val{true}
                old_cell = particles[i].cell
                old_pos = pos
            end

            cell = mesh.cells[particles[i].cell]
            time_remaining, pos, vel = take_advection_step(particles, i, pos, vel, time_remaining, cell, mesh, config, dofs)

            if asserts isa Val{true} && !cell_contains(mesh.cells[particles[i].cell], pos)
                println("WARNING: Lost particle while moving from $old_pos @ cell $old_cell -> $pos @ cell $(particles[i].cell)")
                break
            end
        end
        particles.pos[i] = pos
        particles.vel[i] = vel
    end
end


# This function dynamically dispatches on the cell type
function take_advection_step(particles::ParticleData, i, pos, vel, time_remaining::Float64, cell, mesh, config, dofs)
    side_i, time_to_intersection = find_exit_face(pos, vel, cell.normals, cell.face_origins, dofs)

    if time_remaining <= time_to_intersection
        if dofs isa Val{2}
            pos += time_remaining * SVector(vel[1], vel[2], 0.0)
        else
            pos += time_remaining * vel
        end
        return 0.0, pos, vel
    else
        time_remaining -= time_to_intersection
        if dofs isa Val{2}
            pos += time_to_intersection * SVector(vel[1], vel[2], 0.0)
        else
            pos += time_to_intersection * vel
        end
    end
        
    # Handle boundary condition
    bc_side_indx = cell.bc_side_idx[side_i]
    if bc_side_indx == 0
        # Connected to another cell
        neighbour = cell.neighbours[side_i]
        particles.cell[i] = neighbour
        if neighbour == 0
            error("Particle attempted to leave cell (x=$pos, v=$vel)")
        end
    else
        bc_side = mesh.bc_sides[bc_side_indx]
        bc = config.boundaries[bc_side.bc_index]
        if variantof(bc) == OpenBoundary
            # Absorb the particle; cell = 0 signals deletion to advect!
            particles.cell[i] = 0
            return 0.0, pos, vel
        else
            p = particles[i]
            @reset p.pos = pos
            @reset p.vel = vel
            particles[i] = handle_boundary(p, bc_side, bc_side_indx, config.species[1], variant(bc), config)
            pos = particles.pos[i]
            vel = particles.vel[i]
        end
    end

    # We handle a single step per call of take_advection_step. Break.
    return time_remaining, pos, vel
end


"""
Remove particles absorbed by open boundaries (marked with cell = 0) by compacting in place.
"""
function compact_deleted_particles!(pdata::ParticleData)
    write_idx = 1
    for read_idx in 1:length(pdata)
        if pdata.cell[read_idx] != 0
            if write_idx != read_idx
                pdata[write_idx] = pdata[read_idx]
            end
            write_idx += 1
        end
    end
    resize!(pdata, write_idx - 1)
end


"""
Insert new particles entering the domain through open boundaries.

For each open boundary side, the expected number of simulation particles crossing inward
during `dt` is:

    N = density * face_area * v_flux * dt / mpf

where v_flux is the half-range Maxwell flux speed (`mean_inflow_flux_speed`). The integer
count is obtained by stochastic rounding of N. Each new particle is placed at a random
position on the face, assigned an inflow-weighted velocity, and then advected for a
uniform-random fraction of `dt` to simulate continuous rather than instantaneous insertion.
"""
function insert_open_boundary_particles!(sim::SimulationState, mesh::Mesh, config::SimulationConfig, dt::Float64)
    species = config.species[1]

    for bc_side in mesh.bc_sides
        bc = config.boundaries[bc_side.bc_index]
        variantof(bc) == OpenBoundary || continue
        ob = variant(bc)  # ::OpenBoundary

        v_flux = mean_inflow_flux_speed(bc_side.normal, ob.temperature, ob.velocity, species.mass)
        expected_sim = ob.density * bc_side.area * v_flux * dt / config.mpf

        # Make sure that the correct number of particles are emitted on average, even with low counts
        n_new = floor(Int, expected_sim + rand())

        for _ in 1:n_new
            pos = random_point_on_face(bc_side.face_vertices)
            vel = sample_inflow_maxwellian(bc_side.normal, ob.temperature, ob.velocity, species.mass)

            if config.vrbgk.enabled
                features = (vr_weight = 1.0, last_collided_side = UInt64(0))
            else
                features = (;)
            end

            p = SingleParticle(pos, vel, UInt64(bc_side.adjacent_cell), features)

            # Advect by a random sub-interval of dt to simulate a continuous flux.
            time_remaining = rand() * dt
            while time_remaining > 0.0
                cell = mesh.cells[p.cell]
                p, time_remaining = take_advection_step(p, time_remaining, cell, mesh, config)
                p.cell == 0 && break  # Left through another open boundary; discard
            end

            p.cell != 0 && insert_particle!(sim.particles, p)
        end
    end
end


function vrbgk_advection_finalize(particles::ParticleData, mesh::Mesh, config::SimulationConfig)
    # Calculate the correction factor for each boundary side.
    for i in 1:length(mesh.bc_sides)
        if variantof(config.boundaries[mesh.bc_sides[i].bc_index]) == DiffuseBoundary
            if mesh.bc_sides[i].vrbgk_incident_count == 0
                @atomic mesh.bc_sides[i].vrbgk_incident_sum = NaN
            else
                @atomic mesh.bc_sides[i].vrbgk_incident_sum /= mesh.bc_sides[i].vrbgk_incident_count
            end
        end
    end

    # Apply the correction to each particle
    for i in 1:length(particles)
        if particles[i].features.last_collided_side != 0
            bc_side = mesh.bc_sides[particles[i].features.last_collided_side]
            # This should be impossible, but we will let it pass with a warning anyway.
            if isnan(bc_side.vrbgk_incident_sum)
                # The current value of vr_weight is best left untouched.
                if !config.silent
                    @printf "Warning: No valid particles collided with side %d\n" particles[i].features.last_collided_side
                end
            else
                p = particles[i]
                p = @set p.features.vr_weight = p.features.last_collided_weight * bc_side.vrbgk_incident_sum
                particles[i] = p
            end
        end
    end
end