# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
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

    AK.foreachindex(particles, max_tasks=Threads.nthreads()) do i
        time_remaining = dt
        p = particles[i]

        while time_remaining > 0.0
            cell = mesh.cells[p.cell]
            p_old = p
            p, time_remaining = take_advection_step(p, time_remaining, cell, mesh, config)

            if config.asserts && !cell_contains(mesh.cells[p.cell], p.pos)
                println("WARNING: Lost particle while moving from $(p_old.pos) @ cell $(p_old.cell) -> $(p.pos) @ cell $(p.cell)")
                p = p_old
                break
            end
        end

        particles[i] = p
    end

    # Recount after parallel advection to avoid write races on cell_part_count.
    fill!(sim.cell_part_count, 0)
    for p in sim.particles
        sim.cell_part_count[p.cell] += 1
    end

    # Apply VRBGK fixes
    if config.vrbgk.enabled
        vrbgk_advection_finalize(particles, mesh, config)
    end

    if config.asserts
        assert_particles_in_mesh(sim.particles, mesh)
        assert_cell_part_count(sim.particles, sim.cell_part_count)
    end
end


# This function dynamically dispatches on the cell type
function take_advection_step(p::SingleParticle, time_remaining::Float64, cell, mesh, config) :: Tuple{SingleParticle, Float64}
    # This is the general purpose 3D code. Right now, some things are still assumed to be 2D
    #side_i, time_to_intersection = find_exit_face(p.pos, p.vel, cell.normals, cell.face_origins)
    side_i, time_to_intersection = find_exit_face(p.pos, p.vel, cell.normals, cell.face_origins, Val(config.spatial_dof))

    if time_remaining < time_to_intersection
        p = @set p.pos += time_remaining * SVector(p.vel[1], p.vel[2], 0.0)
        return p, 0.0
    else
        time_remaining -= time_to_intersection
        p = @set p.pos += time_to_intersection * SVector(p.vel[1], p.vel[2], 0.0)
    end
        
    # Handle boundary condition
    bc_side_indx = cell.bc_side_idx[side_i]
    if bc_side_indx == 0
        # Connected to another cell
        neighbour = cell.neighbours[side_i]
        p = @set p.cell = neighbour
        if p.cell == 0
            error("Particle attempted to leave cell (x=$(p.pos), v=$(p.vel))")
        end
    else
        bc_side = mesh.bc_sides[bc_side_indx]
        bc = config.boundaries[bc_side.bc_index]
        p = handle_boundary(p, bc_side, bc_side_indx, config.species[1], variant(bc), config)
    end

    # Handling a single collision is sufficient. Break.
    return p, time_remaining
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
            # Performing this check only when actually required avoid producing an error message in the very common case of
            # low particle counts. In such a case, the incident_sum is both NaN, but also unused since no particle collided.
            if isnan(bc_side.vrbgk_incident_sum)
                @atomic bc_side.vrbgk_incident_sum = 1.0 # Safe-ish default
                if !config.silent
                    @printf "Warning: No valid particles collided with side %d\n" particles[i].features.last_collided_side
                end
            end
            p = particles[i]
            p = @set p.features.vr_weight *= bc_side.vrbgk_incident_sum
            particles[i] = p
        end
    end
end