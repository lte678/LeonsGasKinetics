# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

using Printf
using TOML
using HDF5
using LIKWID
using InteractiveUtils
using ArgParse
using Accessors
using JET

include("constants.jl")
include("sampling/maxwellian.jl")
include("particle_data.jl")
include("simulation_state.jl")
include("mesh.jl")
include("boundary.jl")
include("config.jl")
include("initialization.jl")
include("tracking.jl")
include("moments.jl")
include("assertions.jl")
include("models/atomic.jl")
include("collisions/bgk.jl")
include("advection.jl")
include("simulation.jl")
include("output.jl")


# To ease getting performance metrics.
function precompile(mesh, sim_config)
    part_counts = zeros(length(mesh.cells))
    part_counts[1] = 1
    sim = SimulationState(
        ParticleData([mesh.cells[1].barycenter], [[1000.0, 1000.0, 1000.0]], [1]),
        part_counts,
        0.0
    )
    sim_config = @set sim_config.t_end = sim_config.dt*100
    sim_config = @set sim_config.silent = true
    run_simulation!(sim, mesh, sim_config)
end


function (@main)(args)
    arg_set = ArgParseSettings(description="Leon's Gas Kinetics  Copyright (C) 2025 Leon Teichroeb. This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE.")
    @add_arg_table! arg_set begin
        "--asserts"
            help = "Enables extra assertions (slow!)"
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

    println("Leon's Gas Kinetics  Copyright (C) 2025 Leon Teichroeb")
    println("This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE.")
    println()
    config_path = args["config_file"]
    config = TOML.parsefile(config_path)
    
    # Create the mesh
    meshfile_path = joinpath(dirname(config_path), config["meshfile"])
    mesh = mesh_from_h5(meshfile_path)
    
    # Build the simulation configuration structure
    sim_config = sim_config_from_config(config, dirname(config_path), args["asserts"], mesh.bc_names)

    # This printf is used to precompile.
    @printf "Pre-compiling...\n"
    precompile(mesh, sim_config)

    # Prepare the simulation structure
    sim = SimulationState(
        ParticleData(),
        [],
        0.0  # Start time
    )
    
    # Create and insert the initial particles
    initialize_simulation!(sim, mesh, sim_config, config["initialization"])

    # Output simulation stats
    print_particle_data_info(sim.particles)


    #moments = calc_moments(sim.particles, length(sim.mesh.cells))
    #flow_vars = map(m -> calc_flow_properties(m, specie["mass"]), moments)
    # code_native(accumulate_moments!, (Vector{MomentAccumulator}, ParticleData))
    # @perfmon "FLOPS_SP" calc_moments(sim.particles, length(mesh.cells))

    #println(@report_opt run_simulation!(sim, mesh, sim_config))
    elapsed = @elapsed begin
        #averages = @perfmon "FLOPS_SP" run_simulation!(sim, mesh, sim_config)
        averages = run_simulation!(sim, mesh, sim_config)
    end

    @printf "Simulation finished in %.1f seconds.\n" elapsed
    output_file = @sprintf "%s_DSMCState_%08.4f.h5" sim_config.project_name sim_config.t_end
    @printf "Writing flow state to %s...\n" output_file
    
    write_macro_vals_hdf5(
        joinpath(args["output_dir"], output_file),
        sim_config,
        averages,
        sim_config.t_end,
        true
    )

    # Benchmark
    #@printf "Profiling"
    #sim = SimulationState(ParticleData(), [], 0.0)
    #initialize_simulation!(sim, mesh, sim_config, config["initialization"])
    #print(@time run_simulation!(sim, mesh, sim_config))
end
