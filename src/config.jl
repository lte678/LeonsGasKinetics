# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
    # Project name
    project_name :: String
    # Meshfile
    mesh_file :: String
    # Whether to print output
    silent :: Bool
    # Whether to enable asserts or not
    asserts :: Bool
end


function sim_config_from_config(config, config_dir, asserts, bc_order)
    species = convert(Dict{String, Dict{String, Float64}}, config["species"])
    if length(species) > 1
        error("Multi-species flow is not supported yet.")
    end
    species = map(species_from_config, keys(species), values(species))
    boundaries = Vector{Boundary}()
    for bc in bc_order
        bc_idx = findfirst(b -> b["identifier"] == bc, config["boundary"])
        if bc_idx === nothing
            error("BC $bc missing from simulation config!")
        end
        push!(boundaries, Boundary(boundary_from_config(config["boundary"][bc_idx])))
    end

    return SimulationConfig(
        species,
        boundaries,
        coll_op_from_config(config["dsmc"]["collision_operator"]),
        config["dsmc"]["mpf"],
        config["timestep"]["tend"],
        config["timestep"]["dt"],
        get(config["output"], "sample_fraction", 1.0),
        config["name"],
        joinpath(config_dir, config["meshfile"]),
        false,
        asserts,
    )
end

function boundary_from_config(config)
    type = get(config, "type", "reflective")
    if type == "reflective"
        return ReflectiveBoundary()
    elseif type == "symmetric"
        return SymmetricBoundary()
    elseif type == "diffuse"
        return DiffuseBoundary(
            config["accommodation"],
            config["temperature"],
            config["velocity"],
        )
    else
        error("Unknown boundary type \"$type\"")
    end
end


function coll_op_from_config(operator_name)
    if operator_name == "bgk"
        return bgk_collision!
    elseif operator_name == "none"
        return (part_x, part_v, config, flow_variables, dt) -> ()
    else
        error("Unknown DSMC collision operator \"$operator_name\"")
    end
end