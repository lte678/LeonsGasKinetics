# Leon's Gas Kinetics
# Copyright (C) 2026  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
Per-cell data structure with optional feature fields.

This mirrors the design of `ParticleData` but for cell-level data.
Features are conditionally included based on enabled simulation features
(e.g., VRBGK variance reduction).
"""
struct CellData{Features<:NamedTuple}
    # Number of particles in each cell
    part_count :: Vector{UInt32}
    # Data structures related to gated features, e.g. variance reduction, etc.
    features :: Features
end


"""
Construct a CellData object with optional feature fields.

# Arguments
- `n_cells`: Number of cells in the mesh
- `vrbgk_enabled`: Whether VRBGK variance reduction is enabled

# Returns
A `CellData` object with appropriate feature fields based on the enabled features.
"""
function CellData(n_cells; vrbgk_enabled=false)
    feature_fields = []
    feature_data = []
    # VRBGK-specific per-cell fields for noise reduction
    if vrbgk_enabled
        push!(feature_fields, :vrbgk_ref_temperature)
        push!(feature_data, Vector{Float64}(undef, n_cells))

        push!(feature_fields, :vrbgk_ref_density)
        push!(feature_data, Vector{Float64}(undef, n_cells))
    end

    return CellData(
        zeros(UInt32, n_cells),
        NamedTuple{Tuple(feature_fields)}(feature_data),
    )
end


"""
Resize the CellData to a new number of cells.

This is primarily used during initialization when the mesh size is determined.
"""
function Base.resize!(cdata::CellData, n_cells::Integer)
    resize!(cdata.part_count, n_cells)
    fill!(cdata.part_count, 0)
    
    # Also resize all feature-specific data structures
    for field in cdata.features
        resize!(field, n_cells)
        fill!(field, 0.0)  # Initialize to zero for numeric types
    end
    return cdata
end


"""
Get the number of cells in the CellData.
"""
Base.length(cdata::CellData) = length(cdata.part_count)


"""
Create a scratch CellData with the same structure as the given CellData.
"""
function make_scratch(cdata::CellData)
    return CellData(
        similar(cdata.part_count),
        map(v -> similar(v), cdata.features),
    )
end


"""
Reset all cell data to zero.
"""
function reset!(cdata::CellData)
    fill!(cdata.part_count, 0)
    for field in cdata.features
        fill!(field, 0.0)
    end
    return cdata
end