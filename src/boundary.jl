# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

using StaticArrays
using LightSumTypes
using Atomix


handle_boundary(particle::SingleParticle, normal::SVector{3}, species::SpeciesConfig, boundary::Boundary, config::SimulationConfig) =
    handle_boundary(particle, normal, species, variant(boundary), config)

"""
Handles a "symmetric" boundary.
This is equivalent to a periodic boundary where the cells are self-connected.
Therefore, only one cell's width is permitted in the symmetry direction.
"""
function handle_boundary(
    particle::SingleParticle,
    side::BoundarySide,
    side_idx::UInt32,
    species::SpeciesConfig,
    boundary::SymmetricBoundary,
    config::SimulationConfig
) :: SingleParticle
    # `SingleParticle` is immutable, since this change will not propagate into the `ParticleData` collection 
    particle = @set particle.pos -= 2*(particle.pos ⋅ side.normal)*side.normal
    return particle
end

function handle_boundary(
    particle::SingleParticle,
    side::BoundarySide,
    side_idx::UInt32,
    species::SpeciesConfig,
    boundary::ReflectiveBoundary,
    config::SimulationConfig
) :: SingleParticle
    post_particle = @set particle.vel -= 2*side.normal*dot(side.normal, particle.vel)
    return post_particle
end


function handle_boundary(
    particle::SingleParticle,
    side::BoundarySide,
    side_idx::UInt32,
    species::SpeciesConfig,
    boundary::DiffuseBoundary,
    config::SimulationConfig
) :: SingleParticle
    new_local_vel = sample_wall_distribution(boundary.temperature, species.mass)
    
    tang1, tang2 = perpendicular_pair(side.normal)
    new_vel = new_local_vel[1]*tang1 + new_local_vel[2]*tang2 + new_local_vel[3]*side.normal

    new_vel += boundary.velocity
    particle = @set particle.vel = new_vel

    # Update the variance reduction weight if applicable
    if haskey(particle.features, :vr_weight)
        if particle.features.last_collided_side == 0
            @atomic side.vrbgk_incident_sum += particle.features.vr_weight
            @atomic side.vrbgk_incident_count += UInt32(1)
        end
        particle = @set particle.features.last_collided_side = side_idx
        # Does not contain number density, since this is added at the end of the advection step.
        particle = @set particle.features.vr_weight =
            sqrt(boundary.temperature / config.vrbgk.ref_temperature) *
            maxwellian(config.vrbgk.ref_temperature, zeros(3)         , species.mass, new_vel) /
            maxwellian(        boundary.temperature, boundary.velocity, species.mass, new_vel)
        vrbgk_check_weight(particle.features.vr_weight)
    end

    return particle
end


function perpendicular_pair(n::SVector{3,Float64})
    # Duff et al. 2017 - robust across all normals
    sign_n = copysign(1.0, n[3])
    a = -1.0 / (sign_n + n[3])
    b = n[1] * n[2] * a

    tang1 = SVector(1.0 + sign_n * n[1]^2 * a,  sign_n * b,         -sign_n * n[1])
    tang2 = SVector(b,                             n[2]^2 * a + sign_n, -n[2])
    return tang1, tang2
end


function sample_wall_distribution(wall_temp, mass_ic)
    # Thermal velocity
    cmr = sqrt(BOLTZMANN * wall_temp / mass_ic)
    return SVector(
        cmr * randn(),
        cmr * randn(),
        cmr * sqrt(2*randexp())
    )
end