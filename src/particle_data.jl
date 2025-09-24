# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

using StaticArrays


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


mutable struct ParticleData
    pos :: Vector{SVector{3, Float64}}
    vel :: Vector{SVector{3, Float64}}

    # Cell that the particle lives in. Must be kept up to date with pos_x/y/z
    cell  :: Vector{UInt}
end


struct SingleParticle
    pos :: SVector{3, Float64}
    vel :: SVector{3, Float64}
    cell :: UInt
end


function ParticleData()
    return ParticleData(
        [],
        [],
        [],        
    )
end


function ParticleData(n_particles)
    return ParticleData(
        fill(SVector(0.0, 0.0, 0.0), n_particles),
        fill(SVector(0.0, 0.0, 0.0), n_particles),
        zeros(n_particles),        
    )
end


function pos(data::ParticleData, i)
    return data.pos[i]
end


function vel(data::ParticleData, i)
    return data.vel[i]
end


function Base.getindex(data::ParticleData, i)
    return SingleParticle(
        pos(data, i),
        vel(data, i),
        data.cell[i],
    )
end

function Base.setindex!(data::ParticleData, particle::SingleParticle, i)
    data.pos[i] = particle.pos
    data.vel[i] = particle.vel
    data.cell[i] = particle.cell
end


Base.length(particles::ParticleData) = length(particles.pos)


"""
    insert_particle!(pdata::ParticleData, position::AbstractVector{<:Number}, velocity::AbstractVector{<:Number}, cell)

Insert a particle into the provided `ParticleData` object.

# Arguments
- `pdata::ParticleData`: The particle data structure to insert into
- `position::AbstractVector{<:Number}`: A 3-element vector containing the (x, y, z) position
- `velocity::AbstractVector{<:Number}`: A 3-element vector containing the (x, y, z) velocity
- `cell`: The integer index of the cell this particle is located in.

# Examples
```julia
pdata = ParticleData()
insert_particle!(pdata, [1.0, 2.0, 3.0], [0.1, 0.2, 0.3], 1)
```
"""
function insert_particle!(pdata::ParticleData, position::AbstractVector{<:Number}, velocity::AbstractVector{<:Number}, cell)
    # Validate input vectors
    if length(position) != 3
        throw(ArgumentError("Position vector must have exactly 3 elements (x, y, z)"))
    end
    if length(velocity) != 3
        throw(ArgumentError("Velocity vector must have exactly 3 elements (x, y, z)"))
    end
    
    # Insert position components
    push!(pdata.pos, position)
    
    # Insert velocity components
    push!(pdata.vel, velocity)

    # Insert cell
    push!(pdata.cell, cell)
end


"""
    Sanity-check to make sure that the particle data is coherent.
"""
function check_particle_data(part_data :: ParticleData)
    @assert length(part_data.pos) == length(part_data.vel)
    @assert length(part_data.pos) == length(part_data.cell)
end


function print_particle_data_info(part_data :: ParticleData)
    @printf "Simulation contains %.3e particles.\n" length(part_data.pos) 
    @printf "Datatype is %s. Struct uses %.2e bytes of memory.\n" eltype(part_data.pos) Base.summarysize(part_data)
end
