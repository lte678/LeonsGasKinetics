# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

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


function add_moment!(moment_a :: MomentAccumulator, moment_b :: MomentAccumulator)
    moment_a.c_i += moment_b.c_i
    moment_a.c_ii += moment_b.c_ii
    moment_a.count += moment_b.count
    moment_a.samples += moment_b.samples
end


# 3.3 GFLOPS single core (last test)
function accumulate_moments!(moments :: Vector{MomentAccumulator}, particles :: ParticleData)
    check_particle_data(particles)

    @inbounds for i = 1:length(particles)        
        cell_i = particles.cell[i]
        v = vel(particles, i)
        
        moments[cell_i].count += 1
        moments[cell_i].c_ii  += v.*v
        moments[cell_i].c_i   += v
    end

    for i = 1:length(moments)
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


function calc_flow_properties(moments :: MomentAccumulator, config, cell_volume)
    if moments.count < eps(Float64)
        return FlowProperties()
    end
    
    velocity = moments.c_i / moments.count
    c2 = moments.c_ii / moments.count
    temperature = (config.species[1].mass / BOLTZMANN) * (c2 - velocity.^2)
    sim_part_count = moments.count / moments.samples
    density = sim_part_count * config.mpf / cell_volume
    return FlowProperties(
        velocity,
        temperature,
        sum(temperature) / 3.0,
        density,
        sim_part_count,
    )
end
