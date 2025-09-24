# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

using JET

"""
Performs the advection step.
"""
function advect!(sim::SimulationState, mesh::Mesh, config, dt)
    for i = 1:length(sim.particles)
        # The travel time left for this particle.
        time_remaining = dt
        p = sim.particles[i]
        sim.cell_part_count[p.cell] -= 1
        
        # Keep colliding with further faces until the time_remaining has elapsed.
        while time_remaining > 0.0
            cell = mesh.cells[p.cell]
            p_old = p
            p, time_remaining = take_advection_step(p, time_remaining, cell, config.boundaries, config)

            if config.asserts && !cell_contains(mesh.cells[p.cell], p.pos)
                println("WARNING: Lost particle while moving from $(p_old.pos) @ cell $(p_old.cell) -> $(p.pos) @ cell $(p.cell)")
                p = p_old
                break
            end
        end

        sim.particles[i] = p
        sim.cell_part_count[p.cell] += 1
    end
    
    if config.asserts
        assert_particles_in_mesh(sim.particles, mesh)
        assert_cell_part_count(sim.particles, sim.cell_part_count)
    end
end


# This function dynamically dispatches on the cell type
function take_advection_step(p::SingleParticle, time_remaining::Float64, cell, boundaries, config) :: Tuple{SingleParticle, Float64}
    if all(abs.(p.vel) .< eps(0.0))
        return p, 0.0
    end
    
    # Check for collisions with each face.
    sides = get_sides(cell)
    for side_i = 1:length(sides)
        side = sides[side_i]
        bc_indx = cell.bcs[side_i]
        neighbour = cell.neighbours[side_i]
        hit, intersection = ray_intersects_rect(p.pos, p.vel, side)
        if !hit
            continue
        end
        # Advance the particle position
        time_to_intersection = norm(intersection - p.pos) / sqrt(p.vel[1]^2 + p.vel[2]^2 + p.vel[3]^2)
        if time_remaining < time_to_intersection + eps(0.0)
            p = @set p.pos += time_remaining * p.vel
            return p, 0.0
        else
            time_remaining -= time_to_intersection
            p = @set p.pos = intersection
        end
        
        # Handle boundary condition
        if bc_indx == 0
            # Connected to another cell
            p = @set p.cell = neighbour
            if p.cell == 0
                error("Particle attempted to leave cell (x=$(p.pos), v=$(p.vel))")
            end
        else
            bc = boundaries[bc_indx]
            n = side_normal(side)
            p = handle_boundary(p, n, config.species[1], bc)
        end
 
        # Handling a single collision is sufficient. Break.
        return p, time_remaining
    end

    error("Particle (x=$(p.pos), v=$(p.vel)) does not intersect any faces.")
end