# Leon's Gas Kinetics
# Copyright (C) 2026  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
Retrieves a value from `dict[key]`. Validates that the value is within `options`.
Returns the value as a `Symbol`.

# Arguments
- `dict`: The dictionary to search.
- `key`: The key to lookup.
- `options`: A collection of strings representing valid choices.
- `default`: (Optional) The default string if the key is missing.
"""
function get_option(dict, key, options::AbstractVector{String}; default=nothing)
    val = get(dict, key, default)

    if isnothing(val)
        error("Missing required configuration key: \"$key\". Must be one of $options")
    end
    if !(val ∈ options)
        error("Unknown option \"$val\" for \"$key\". Must be ∈ $options.")
    end
    
    return Symbol(val)
end

function get_option(dict, key, options::AbstractVector{<:Number}; default=nothing)
    val = get(dict, key, default)

    if isnothing(val)
        error("Missing required configuration key: \"$key\". Must be one of $options")
    end
    if !(val ∈ options)
        error("Unknown option \"$val\" for \"$key\". Must be ∈ $options.")
    end
    
    return val
end


struct SymmetricBoundary end
struct ReflectiveBoundary end
struct DiffuseBoundary
    accommodation :: Float64
    temperature :: Float64
    # Global coordinates
    velocity :: SVector{3, Float64}
end
struct OpenBoundary
    temperature :: Float64
    density :: Float64
    # Global coordinates
    velocity :: SVector{3, Float64}
end

@sumtype Boundary(SymmetricBoundary, ReflectiveBoundary, DiffuseBoundary, OpenBoundary)


struct VRBGKConfig
    enabled :: Bool
    ref_temperature :: Float64
    ref_density :: Float64
    adaptive_equilibrium :: Bool
    adaptive_smoothing_factor :: Float64
end


struct SimulationConfig{T1<:Function}
    # Species definition
    species :: Vector{SpeciesConfig}
    # Boundaries. In the same order as the numerical index in the mesh sides.
    boundaries :: Vector{Boundary}
    # Function to use as the collision operator
    collision_operator :: T1
    # Macro particle factor
    mpf :: Float64
    # Final time of simulation
    t_end :: Float64
    # Timestep size
    dt :: Float64
    # Sampling fraction
    sample_fraction :: Float64
    # Number of iterations between file outputs
    output_interval :: Int64
    # Reporting inveral for console output
    report_interval :: Float64
    # Project name
    project_name :: String
    # Meshfile
    mesh_file :: String
    # Output path
    output_path :: String
    # Whether to print output
    silent :: Bool
    # Whether to enable asserts or not
    asserts :: Bool
    # Settings for VRBGK
    vrbgk :: VRBGKConfig
    # Sample accumulation mode: :continuous (entire sim) or :per_interval (reset on each output)
    accumulation_mode :: Symbol
    # 3D, 2D with z=0
    degrees_of_freedom :: Int64
end


function sim_config_from_config(config, config_dir, output_path, asserts) :: SimulationConfig
    species = convert(Dict{String, Dict{String, Float64}}, config["species"])
    if length(species) > 1
        error("Multi-species flow is not supported yet.")
    end
    species = map(species_from_config, keys(species), values(species))

    accumulation_mode = get_option(
        config["output"],
        "accumulation",
        ["continuous", "per-interval"],
        default="continuous"
    )

    degrees_of_freedom = get_option(
        config, "degrees_of_freedom", [2, 3], default=3
    )

    coll_op = get_coll_op(
        get_option(config["dsmc"], "collision_operator", ["bgk", "none"])
    )

    if haskey(config, "denoise") && config["denoise"]["enabled"]
        vrbgk_config = VRBGKConfig(
            true,
            Float64(config["denoise"]["T_ref"]),
            Float64(config["denoise"]["n_ref"]),
            get(config["denoise"], "adaptive_equilibrium", false),
            get(config["denoise"], "adaptive_smoothing_factor", 0.9),
        )
    else
        vrbgk_config = VRBGKConfig(false, 0.0, 0.0, false, 0.9)
    end

    # TODO: The boundaries are loaded separately. This should be improved.
    return SimulationConfig(
        species,
        Vector{Boundary}(),
        coll_op,
        Float64(config["dsmc"]["mpf"]),
        Float64(config["timestep"]["tend"]),
        Float64(config["timestep"]["dt"]),
        get(config["output"], "sample_fraction", 1.0),
        get(config["output"], "output_interval", 0),
        get(config["output"], "report_interval", 5.0),
        config["name"],
        joinpath(config_dir, config["meshfile"]),
        output_path,
        false,
        asserts,
        vrbgk_config,
        accumulation_mode,
        degrees_of_freedom
    )
end


function boundaries_from_config(config, boundary_order)
    boundaries = Vector{Boundary}()
    for bc in boundary_order
        bc_idx = findfirst(b -> lowercase(b["identifier"]) == bc, config["boundary"])
        if bc_idx === nothing
            error("BC '$bc' missing from simulation config!")
        end
        push!(boundaries, Boundary(boundary_from_config(config["boundary"][bc_idx])))
    end
    return boundaries
end


function boundary_from_config(config)
    type = get_option(config, "type", ["reflective", "symmetric", "diffuse", "open"], default="reflective")
    if type == :reflective
        return ReflectiveBoundary()
    elseif type == :symmetric
        return SymmetricBoundary()
    elseif type == :diffuse
        return DiffuseBoundary(
            config["accommodation"],
            config["temperature"],
            get(config, "velocity", SVector(0.0, 0.0, 0.0)),
        )
    elseif type == :open
        return OpenBoundary(
            Float64(config["temperature"]),
            Float64(config["density"]),
            SVector{3,Float64}(get(config, "velocity", [0.0, 0.0, 0.0])),
        )
    end
end


function get_coll_op(operator_name)
    if operator_name == :bgk
        return bgk_collision!
    elseif operator_name == :none
        return (particle_data, cell_data, samples, config, flow_variables, dt) -> ()
    else
        error("Unknown DSMC collision operator \"$operator_name\"")
    end
end