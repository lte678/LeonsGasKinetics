# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
Takes the initial_state and propagates it in time until the stopping condition is reached.
"""
function run_simulation!(initial_state::SimulationState, mesh::Mesh, config)
    state = initial_state
    dt = config.dt
    sampling_start_time = (1.0 - config.sample_fraction) * config.t_end
    n_cells = length(mesh.cells)

    # Avoid reallocation inside loop
    moments = [MomentAccumulator() for _ in 1:n_cells]
    moments_avg = [MomentAccumulator() for _ in 1:n_cells]
    # flow_vars = [FlowProperties() for _ in 1:length(mesh.cells)]

    particle_reordering = zeros(Int, length(state.particles))
    _, t = Base.Sort.make_scratch(nothing, eltype(particle_reordering), length(particle_reordering))
    sorted_particles = ParticleData(length(state.particles))

    iteration = 1
    while state.time < config.t_end
        # Gracefully handle the final timestep. Make sure we dont make a zero length step afterwards.
        if state.time + dt > config.t_end
            dt = config.t_end - state.time
            state.time = config.t_end
        else
            state.time += dt
        end
        
        # Advect the particles. This also updates the particle cell information.
        advect!(state, mesh, config, dt)


        # Resort the particles. This is to make the BGK collision routine much simpler and to improve cache locality
        # for the moment calculation for example.
        sortperm!(particle_reordering, state.particles.cell; scratch=t)
        sorted_particles.pos .= state.particles.pos[particle_reordering]
        sorted_particles.vel .= state.particles.vel[particle_reordering]
        sorted_particles.cell .= state.particles.cell[particle_reordering]

        # Calculate moments
        accumulate_moments!(moments, sorted_particles)
        
        # Add the moments to the time average
        if state.time > sampling_start_time
            if state.time - dt < sampling_start_time && !config.silent
                @printf "Starting sampling.\n"
            end
            for i = 1:length(moments)
                add_moment!(moments_avg[i], moments[i])
            end
        end
    
        # Perform the collision step
        particle_start_idx = 1
        for i in 1:n_cells
            flow_vars = calc_flow_properties(moments[i], config, mesh.cells[i].volume)
            particle_x = @view sorted_particles.pos[particle_start_idx:particle_start_idx + state.cell_part_count[i] - 1]
            particle_v = @view sorted_particles.vel[particle_start_idx:particle_start_idx + state.cell_part_count[i] - 1]
            config.collision_operator(particle_x, particle_v, config, flow_vars, dt)
            particle_start_idx += state.cell_part_count[i]
        end

        # Clear moments
        clear_moments!(moments)
        
        # Swap sorting array
        state.particles.pos .= sorted_particles.pos
        state.particles.vel .= sorted_particles.vel
        state.particles.cell .= sorted_particles.cell

        if iteration % 1000 == 0 && !config.silent
            @printf "[Iteration = %6d]\n" iteration
        end

        iteration += 1
    end

    return map((m, c) -> calc_flow_properties(m, config, c.volume), moments_avg, mesh.cells)
end