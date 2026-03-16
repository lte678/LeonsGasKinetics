# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

using StaticArrays
using KernelAbstractions


struct SpeciesConfig
    mass :: Float64
    # Diameter at T_ref
    d_ref :: Float64
    # Reference temperature for diameter
    T_ref :: Float64
    # TODO
    omega :: Float64
    name :: String
end


function species_from_config(name, config)
    return SpeciesConfig(
        config["mass"],
        config["d_ref"],
        config["T_ref"],
        config["omega"],
        name
    )
end


struct ParticleData{Features<:NamedTuple}
    pos :: Vector{SVector{3, Float64}}
    vel :: Vector{SVector{3, Float64}}
    # Cell that the particle lives in. Must be kept in sync with pos_x/y/z
    cell  :: Vector{UInt64}

    # Data structures related to gated features, e.g. variance reduction, chemistry, etc.
    features :: Features
end


function ParticleData(; vrbgk_enabled=false)
    feature_fields = []
    feature_data = []
    if vrbgk_enabled
        push!(feature_fields, :vr_weight, :last_collided_side)
        push!(feature_data, Vector{Float64}(), Vector{UInt64}())
    end

    return ParticleData(
        Vector{SVector{3,Float64}}(),
        Vector{SVector{3,Float64}}(),
        Vector{UInt64}(),
        NamedTuple{Tuple(feature_fields)}(feature_data),
    )
end


function ParticleData(n_particles; args...)
    pdata = ParticleData(; args...)
    resize!(pdata, n_particles)
    return pdata
end


struct SingleParticle{Features<:NamedTuple}
    pos :: SVector{3, Float64}
    vel :: SVector{3, Float64}
    cell :: UInt64
    features :: Features
end

function SingleParticle(pos::AbstractVector{<:Real}, vel::AbstractVector{<:Real}, id::Integer, features::NamedTuple)
    SingleParticle(
        SVector{3, Float64}(pos),
        SVector{3, Float64}(vel),
        UInt64(id),
        features
    )
end

# SingleParticle(pos, vel, cell; features...) = return SingleParticle(pos, cell, cell, features)

Base.copy(p::SingleParticle) = SingleParticle(p.pos, p.vel, p.cell, p.features)


function Base.resize!(pdata::ParticleData, n_particles::Integer)
    resize!(pdata.pos, n_particles)
    resize!(pdata.vel, n_particles)
    resize!(pdata.cell, n_particles)
    # Also resize all feature-specific data structures
    for field in pdata.features
        resize!(field, n_particles)
    end
    return pdata
end


function Base.getindex(data::ParticleData, i)
    feature_data = map(v -> v[i], data.features)
    return SingleParticle(
        data.pos[i],
        data.vel[i],
        data.cell[i],
        feature_data,
    )
end

function Base.setindex!(data::ParticleData, particle::SingleParticle, i)
    data.pos[i] = particle.pos
    data.vel[i] = particle.vel
    data.cell[i] = particle.cell

    for (k, v) in pairs(particle.features)
        data.features[k][i] = v
    end
    
    return data
end

"""
Copy entire ParticleData
"""
function Base.copyto!(dst::ParticleData, src::ParticleData)
    dst.pos   .= src.pos
    dst.vel   .= src.vel
    dst.cell  .= src.cell
    for fname in keys(src.features)
        dst.features[fname] .= src.features[fname]
    end
    return dst
end

"""
Reorder particles into `dst` using permutation `perm`.
"""
function reorder!(dst::ParticleData, src::ParticleData, perm::AbstractVector{<:Integer})
    dst.pos  .= view(src.pos,  perm)
    dst.vel  .= view(src.vel,  perm)
    dst.cell .= view(src.cell, perm)
    for fname in keys(src.features)
        dst.features[fname] .= view(src.features[fname], perm)
    end
    return dst
end

"""
Reorder particles in pdata according to the permutation. Uses pre-allocated scratch to avoid allocation.
"""
function reorder!(pdata::ParticleData, perm::AbstractVector{<:Integer}; scratch::ParticleData)
    @assert scratch !== pdata
    reorder!(scratch, pdata, perm)
    copyto!(pdata, scratch)
    return pdata
end


Base.length(particles::ParticleData) = length(particles.pos)

function Base.iterate(pdata::ParticleData, i=1)
    i > length(pdata) && return nothing
    return (pdata[i], i + 1)
end

Base.keys(pdata::ParticleData) = keys(pdata.pos)

KernelAbstractions.get_backend(pdata::ParticleData) = get_backend(pdata.pos)

"""
Insert a particle into the provided `ParticleData` object.
"""
function insert_particle!(pdata::ParticleData, p::SingleParticle)
    resize!(pdata, length(pdata) + 1)
    pdata[length(pdata)] = p
end


function make_scratch(pdata::ParticleData)
    return ParticleData(
        similar(pdata.pos),
        similar(pdata.vel),
        similar(pdata.cell),
        map(v -> similar(v), pdata.features),
    )
end


"""
Sanity-check to make sure that the particle data is coherent.
"""
function check_particle_data(part_data :: ParticleData)
    @assert length(part_data.pos) == length(part_data.vel)
    @assert length(part_data.pos) == length(part_data.cell)
    for feature_arr in part_data.features
        @assert length(part_data.pos) == length(feature_arr)
    end
end


function print_particle_data_info(part_data :: ParticleData)
    @printf "Simulation contains %.3e particles.\n" length(part_data.pos) 
    @printf "Datatype is %s. Struct uses %.2e bytes of memory.\n" eltype(eltype(part_data.pos)) Base.summarysize(part_data)
end


struct ParticleDataView
    data :: ParticleData
    offset :: UInt
    stop_index :: UInt
end

function particle_view(data::ParticleData, start_index, end_index)
    @assert end_index <= length(data)
    # end_index may be one smaller than start_index if the view is zero-length.
    @assert end_index + 1 >= start_index
    # The offset needs to be one smaller than the start index. Example: Starting at 1 requires adding 0 to index.
    return ParticleDataView(data, start_index - 1, end_index)
end

Base.length(view::ParticleDataView) = Int(view.stop_index - view.offset)

function Base.getindex(view::ParticleDataView, i::Integer)
    @boundscheck 1 <= i <= length(view) || throw(BoundsError(view, i))
    return view.data[view.offset + i]
end

function Base.setindex!(view::ParticleDataView, particle::SingleParticle, i::Integer)
    @boundscheck 1 <= i <= length(view) || throw(BoundsError(view, i))
    view.data[view.offset + i] = particle
    return view
end

function Base.iterate(view::ParticleDataView, i=1)
    i > length(view) && return nothing
    return (view[i], i + 1)
end

Base.eltype(::ParticleDataView) = SingleParticle
Base.firstindex(::ParticleDataView) = 1
Base.lastindex(view::ParticleDataView) = length(view)