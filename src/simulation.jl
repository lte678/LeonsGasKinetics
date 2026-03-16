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
    if config.vrbgk.enabled
        cell_moments = [VRBGKMomentAccumulator() for _ in 1:n_cells]
    else
        cell_moments = [MomentAccumulator() for _ in 1:n_cells]
    end
    cell_accumulators = [DefaultDict{Symbol, Averager}(() -> Averager()) for _ in 1:n_cells]

    particle_reordering = zeros(Int, length(state.particles))
    _, t = Base.Sort.make_scratch(nothing, eltype(particle_reordering), length(particle_reordering))
    particle_data_scratch = ParticleData(length(state.particles); vrbgk_enabled=config.vrbgk.enabled)

    iteration = 1
    reporter = ProgressReporter(config.report_interval)
    while state.time < config.t_end
        # Gracefully handle the final timestep. Make sure we dont make a zero length step afterwards.
        if state.time + dt > config.t_end
            dt = config.t_end - state.time
            state.time = config.t_end
        else
            state.time += dt
        end
        
        # Advect the particles. This also updates the particle cell information.
        timed_region(state.perf_counters, :advection) do
            advect!(state, mesh, config, dt)
        end

        # Re-sort the particles. This is to make the BGK collision routine much simpler and to improve cache locality
        # for the moment calculation for example.
        timed_region(state.perf_counters, :sorting) do
            sortperm!(particle_reordering, state.particles.cell; scratch=t)
            reorder!(state.particles, particle_reordering; scratch=particle_data_scratch)
        end

        # Calculate moments
        timed_region(state.perf_counters, :accumulate_moments) do
            accumulate_moments!(cell_moments, state.particles)
        end
    
        # Perform the collision step
        timed_region(state.perf_counters, :collision) do
            particle_start_idx = 1
            for i in 1:n_cells
                flow_vars = calc_flow_properties(cell_moments[i], config, mesh.cells[i].volume)
                cell_particles = particle_view(state.particles, particle_start_idx, particle_start_idx + state.cell_part_count[i] - 1)

                # Perform collision on the current cell's particles.
                config.collision_operator(cell_particles, cell_accumulators[i], config, flow_vars, dt)
                
                particle_start_idx += state.cell_part_count[i]
            end
        end

        # Add the moments to the time average
        if state.time > sampling_start_time
            if state.time - dt < sampling_start_time && !config.silent
                @printf "Starting sampling.\n"
            end
            for i in 1:n_cells
                flow_vars = calc_flow_properties(cell_moments[i], config, mesh.cells[i].volume)
                add_sample!(cell_accumulators[i][:u_x] , flow_vars.velocity[1])
                add_sample!(cell_accumulators[i][:u_y] , flow_vars.velocity[2])
                add_sample!(cell_accumulators[i][:u_z] , flow_vars.velocity[3])
                add_sample!(cell_accumulators[i][:T_x], flow_vars.temperature[1])
                add_sample!(cell_accumulators[i][:T_y], flow_vars.temperature[2])
                add_sample!(cell_accumulators[i][:T_z], flow_vars.temperature[3])
                add_sample!(cell_accumulators[i][:density], flow_vars.density)
                add_sample!(cell_accumulators[i][:sim_part_count], flow_vars.sim_particle_count)
            end
        end

        # Clear moments
        clear_moments!(cell_moments)

        # Print report every config.report_interval seconds.
        config.silent || maybe_report!(reporter, iteration, state.time)

        iteration += 1
    end

    return cell_accumulators
end

"""
Load configuration, run simulation, and return postprocessed volume data without 
writing HDF5 output.
"""
function run_simulation_from_config(config_path::String; enable_asserts=false)
    @assert isfile(config_path) "Config file not found: $config_path"
    
    config = TOML.parsefile(config_path)
    config_dir = dirname(config_path)
    
    # Create mesh
    meshfile_path = joinpath(config_dir, config["meshfile"])
    mesh = mesh_from_h5(meshfile_path)
    
    # Build simulation configuration
    sim_config = sim_config_from_config(config, config_dir, enable_asserts, mesh.bc_names)
    
    # Initialize simulation
    sim = SimulationState(
        ParticleData(; vrbgk_enabled=sim_config.vrbgk.enabled),
        Vector{UInt32}(),
        0.0,
        PerformanceCounters()
    )
    initialize_simulation!(sim, mesh, sim_config, config["initialization"])
    
    volume_samples = run_simulation!(sim, mesh, sim_config)
    
    # Postprocess to get macroscopic quantities (velocity, temperature, etc.)
    volume_output = postprocess(volume_samples, sim_config, mesh)
    
    return volume_output, sim_config, mesh
end


mutable struct ProgressReporter
    last_report_time::Float64
    report_interval::Float64
end

ProgressReporter(interval::Float64) = ProgressReporter(time(), interval)

function maybe_report!(reporter::ProgressReporter, iteration::Int, sim_time::Float64)
    now = time()
    if now - reporter.last_report_time >= reporter.report_interval
        @printf "[Iteration = %6d] t = %8.3fμs\n" iteration sim_time * 1e6
        reporter.last_report_time = now
    end
end