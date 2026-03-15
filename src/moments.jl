# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# This file contains functions for the calculation of moments, such as velocity and temperature.
# Sampling averaging is considered a separate task and is handled in averaging.jl

using StaticArrays

struct FlowProperties
    velocity :: SVector{3, Float64}
    temperature :: SVector{3, Float64}
    mean_temperature :: Float64
    density :: Float64
    sim_particle_count :: Float64
end

function FlowProperties()
    return FlowProperties(SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), 0.0, 0.0, 0.0)
end


mutable struct MomentAccumulator
    c_i  :: SVector{3, Float64}
    c_ii :: SVector{3, Float64}

    # Count may be any sort of weight used to normalize the moments by M/count.
    count :: Float64
    # Number of averaging steps.
    samples :: Int
end

function MomentAccumulator()
    return MomentAccumulator([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.0, 0)
end


# Variance reduced moments
mutable struct VRBGKMomentAccumulator
    # Sum of (1-W) * v_i
    c_i  :: SVector{3, Float64}
    # Sum of (1-W) * v_i²        
    c_ii :: SVector{3, Float64}
    # Sum of W (variance reduction weights)
    vr_sum :: Float64

    # Simulation particle count
    count :: Int64
    # Number of averaging steps.               
    samples::Int64
end

function VRBGKMomentAccumulator()
    return VRBGKMomentAccumulator([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.0, 0, 0)
end


function add_moment!(moment_a::MomentAccumulator, moment_b::MomentAccumulator)
    moment_a.c_i += moment_b.c_i
    moment_a.c_ii += moment_b.c_ii
    moment_a.count += moment_b.count
    moment_a.samples += moment_b.samples
end


function add_moment!(moment_a::VRBGKMomentAccumulator, moment_b::VRBGKMomentAccumulator)
    moment_a.c_i += moment_b.c_i
    moment_a.c_ii += moment_b.c_ii
    moment_a.vr_sum += moment_b.vr_sum
    moment_a.count += moment_b.count
    moment_a.samples += moment_b.samples
end


# 3.3 GFLOPS single core (last test)
function accumulate_moments!(moments::Vector{MomentAccumulator}, particles::ParticleData)
    check_particle_data(particles)

    @inbounds for i = 1:length(particles)        
        cell_i = particles.cell[i]
        v = particles.vel[i]
        
        moments[cell_i].count += 1
        moments[cell_i].c_ii  += v.*v
        moments[cell_i].c_i   += v
    end

    for i = 1:length(moments)
        moments[i].samples += 1
    end
end

# Variance-reduced moment calculation
function accumulate_moments!(moments::Vector{VRBGKMomentAccumulator}, particles::ParticleData)
    check_particle_data(particles)

    @inbounds for i in 1:length(particles)
        cell_i = particles.cell[i]
        v = particles.vel[i]
        w = particles.features.vr_weight[i]
        
        vr_factor = 1.0 - w
        
        moments[cell_i].count += 1
        moments[cell_i].vr_sum += w
        moments[cell_i].c_ii += vr_factor .* v.*v
        moments[cell_i].c_i += vr_factor .* v
    end
    
    @inbounds for i in 1:length(moments)
        moments[i].samples += 1
    end
end


function clear_moments!(moments::Vector{MomentAccumulator})
    for i = 1:length(moments)
        moments[i].c_i = SVector(0.0, 0.0, 0.0)
        moments[i].c_ii = SVector(0.0, 0.0, 0.0)
        moments[i].count = 0.0
        moments[i].samples = 0
    end
end

function clear_moments!(moments::Vector{VRBGKMomentAccumulator})
    for i = 1:length(moments)
        moments[i].c_i = SVector(0.0, 0.0, 0.0)
        moments[i].c_ii = SVector(0.0, 0.0, 0.0)
        moments[i].vr_sum = 0.0
        moments[i].count = 0.0
        moments[i].samples = 0
    end
end



"""
Calculates macroscopic flow properties based on the mean velocity, squared velocity and particle number.
`c_i` is the sum of velocities, while `c_ii` is the sum of square velocities.
"""
function calc_flow_properties(count, c_i, c_ii, config, cell_volume)
    temperature = (config.species[1].mass / BOLTZMANN) * (c_ii - c_i.^2)
    sim_part_count = count
    density = sim_part_count * config.mpf / cell_volume
    
    return FlowProperties(
        c_i,
        temperature,
        sum(temperature) / 3.0,
        density,
        sim_part_count,
    )
end

function calc_flow_properties(moments :: MomentAccumulator, config, cell_volume)
    if moments.count < eps(Float64)
        return FlowProperties()
    end
    calc_flow_properties(
        moments.count / moments.samples,
        moments.c_i / moments.count,
        moments.c_ii / moments.count,
        config,
        cell_volume
    )
end

# Variance-reduced flow properties
function calc_flow_properties(moments::VRBGKMomentAccumulator, config, cell_volume)
    mass = config.species[1].mass
    
    # Global equilibrium properties from config
    n_eq = config.vrbgk.ref_density      # Reference number density [1/m³]
    T_eq = config.vrbgk.ref_temperature  # Reference temperature [K]
    
    # Physical particles represented by equilibrium distribution
    count_eq = n_eq * cell_volume
    
    # Particle count
    count_vr = moments.count - moments.vr_sum + count_eq

    # Velocity
    c_i = moments.c_i ./ moments.count
    c_ii = moments.c_ii ./ moments.count .+ (count_eq / count_vr) * BOLTZMANN * T_eq / mass
    
    # Thermal velocity squared
    temperature = (config.species[1].mass / BOLTZMANN) * (c_ii - c_i.^2)
    
    # Number density
    density = count_vr / cell_volume
    sim_part_count = moments.count / moments.samples

    return FlowProperties(
        c_i,
        temperature,
        sum(temperature) / 3.0,
        density,
        sim_part_count,
    )
end