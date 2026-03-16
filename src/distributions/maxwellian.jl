# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

using StaticArrays

"""
Generate a single Maxwell-Boltzmann distributed velocity vector.
"""
sample_maxwellian(temperature, bulk_velocity, mass) = bulk_velocity .+ sqrt(BOLTZMANN * temperature / mass) * randn(3)


"""
Generate Maxwell-Boltzmann distributed velocities.
"""
function sample_maxwellian(temperature, bulk_velocity, mass, n_particles)
    velocities = Vector{SVector{3, Float64}}(undef, n_particles)
    sample_maxwellian!(velocities, temperature, bulk_velocity, mass)
    return velocities
end


"""
Generate Maxwell-Boltzmann distributed velocities.
"""
function sample_maxwellian!(velocities, temperature, bulk_velocity, mass)   
    # Thermal velocity (most probable speed in Maxwellian distribution)
    v_th = sqrt(BOLTZMANN * temperature / mass)
    for i in 1:length(velocities)
        # Total velocity = bulk velocity + thermal velocity
        velocities[i] = bulk_velocity .+ v_th * randn(3)
    end
end


"""
Calculate probability density of maxwellian distribution.
"""
function maxwellian(temperature, bulk_velocity, mass, particle_velocity::AbstractArray{<:Number})
    v_sq = sum((particle_velocity - bulk_velocity).^2)
    return (mass / (2 * PI * BOLTZMANN * temperature))^1.5 * exp(-(mass*v_sq)/(2*BOLTZMANN*temperature))
end


"""
Calculate probability density of maxwellian distribution.
"""
function maxwellian(temperature, bulk_velocity, mass, particle_velocity)
    p = Vector{Float64}(undef, length(particle_velocity))
    maxwellian!(p, temperature, bulk_velocity, mass, particle_velocity)
    return p
end


"""
Calculate probability density of maxwellian distribution.
"""
function maxwellian!(probability, temperature, bulk_velocity, mass, particle_velocity)
    for i in 1:length(probability)
        probability[i] = maxwellian(temperature, bulk_velocity, mass, particle_velocity[i])
    end
end