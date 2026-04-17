# Leon's Gas Kinetics
# Copyright (C) 2026  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

module LeonsGasKinetics

using Printf
using TOML
using HDF5
using LIKWID
using InteractiveUtils
using ArgParse
using Accessors
using JET
using DataStructures

include("constants.jl")
include("statistics.jl")
include("distributions/maxwellian.jl")
include("distributions/inflow_maxwellian.jl")
include("distributions/wall_maxwellian.jl")
include("particle_data.jl")
include("vrbgk.jl")
include("performance_counters.jl")
include("simulation_state.jl")
include("mesh.jl")
include("config.jl")
include("boundaries/boundary.jl")
include("initialization.jl")
include("tracking.jl")
include("moments.jl")
include("assertions.jl")
include("models/atomic.jl")
include("collisions/bgk.jl")
include("advection.jl")
include("simulation.jl")
include("output/postprocessing.jl")
include("output/hdf5_output.jl")
include("output/output.jl")


# To ease getting performance metrics.
function precompile(mesh, sim_config)
    # Prepare particle data.
    part_counts = zeros(UInt32, length(mesh.cells))
    part_counts[1] = 2
    part_data = ParticleData(; vrbgk_enabled=sim_config.vrbgk.enabled)
    if sim_config.vrbgk.enabled
        feature_data = (; :vr_weight => 1.0)
    else
        feature_data = (;)
    end
    
    insert_particle!(
        part_data,
        SingleParticle(mesh.cells[1].barycenter, [100.0, 100.0, 100.0], 1, feature_data)
    )
    insert_particle!(
        part_data,
        SingleParticle(mesh.cells[1].barycenter, [-100.0, -100.0, -100.0], 1, feature_data)
    ) 
    
    sim = SimulationState(
        part_data,
        part_counts,
        0.0,
        PerformanceCounters()
    )
    sim_config = @set sim_config.t_end = sim_config.dt*100
    sim_config = @set sim_config.silent = true
    sim_config = @set sim_config.output_interval = 0  # No outputs
    run_simulation!(sim, mesh, sim_config)
end


function run(args)
    arg_set = ArgParseSettings(description="Leon's Gas Kinetics  Copyright (C) 2026 Leon Teichroeb. This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE.")
    @add_arg_table! arg_set begin
        "--asserts"
            help = "Enables extra assertions (slow!)"
            action = :store_true
        "--profile"
            help = "Print how compute time was spent in the simulation"
            action = :store_true
        "config_file"
            help = "Path to simulation configuration file"
            required = true
        "output_dir"
            help = "Output path"
    end
    args = parse_args(args, arg_set)
    if args["output_dir"] === nothing
        args["output_dir"] = joinpath(dirname(args["config_file"]), "output") 
    end
    if !ispath(args["output_dir"])
        mkdir(args["output_dir"])
    end

    println("Leon's Gas Kinetics  Copyright (C) 2026 Leon Teichroeb")
    println("This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE.")
    println()

    println("Executing on CPU with $(Threads.nthreads()) threads.")
    config_path = args["config_file"]
    if !isfile(config_path)
        error("Failed to find config file $config_path")
    end
    config = TOML.parsefile(config_path)
    
    # Create the mesh
    meshfile_path = joinpath(dirname(config_path), config["meshfile"])
    mesh = mesh_from_h5(meshfile_path)
    
    # Build the simulation configuration structure
    sim_config = sim_config_from_config(config, dirname(config_path), args["output_dir"], args["asserts"], mesh.bc_names)

    # This printf is used to precompile.
    if !sim_config.silent
        @printf "Pre-compiling...\n"
    end
    precompile(mesh, sim_config)

    # Prepare the simulation structure
    sim = SimulationState(
        ParticleData(; vrbgk_enabled=sim_config.vrbgk.enabled),
        Vector{UInt32}(),
        0.0,  # Start time
        PerformanceCounters()
    )
    
    # Create and insert the initial particles
    timed_region(sim.perf_counters, :initialization) do
        initialize_simulation!(sim, mesh, sim_config, config["initialization"])
    end

    # Output simulation stats
    print_particle_data_info(sim.particles)

    elapsed = @elapsed begin
        # For profiling the average flops
        #volume_samples = @perfmon "FLOPS_SP" run_simulation!(sim, mesh, sim_config)
        volume_samples = run_simulation!(sim, mesh, sim_config)
    end

    @printf "Simulation finished in %.1f seconds.\n" elapsed
    
    write_flow_state(volume_samples, sim_config.t_end, mesh, sim_config)

    if args["profile"]
        println()
        print_performance_stats(sim.perf_counters)
    end
    # Benchmark
    #@printf "Profiling"
    #sim = SimulationState(ParticleData(), [], 0.0)
    #initialize_simulation!(sim, mesh, sim_config, config["initialization"])
    #print(@time run_simulation!(sim, mesh, sim_config))
end

export run_simulation!, run_simulation_from_config

end