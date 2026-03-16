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
    vmag2 = particle.vel[1]^2 + particle.vel[2]^2 + particle.vel[3]^2
    new_local_vel = sample_wall_distribution(vmag2, boundary.temperature, boundary.accommodation, species.mass)
    # TODO: Produce perpendicular vector pair according to normal (critical!)
    new_vel = new_local_vel[1]*[0, 1, 0] + new_local_vel[2]*[0, 0, 1] + new_local_vel[3]*side.normal
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
    end

    return particle
end


"""
    wall_distribution(spec_id, velo_square, wall_temp, trans_acc, mass_ic) -> Vector{Float64}

Compute the post-wall-collision velocity vector for a reflected particle.
The reflection model uses the wall temperature `wall_temp` and the
translational accommodation coefficient `trans_acc`.  The kinetic energy
is adjusted stochastically and the outgoing direction is uniformly
randomised in the tangent plane.

# Arguments
- `spec_id`    : species identifier (unused except for dispatch)
- `velo_square`: squared magnitude of the incoming velocity
- `wall_temp`  : wall temperature (K)
- `trans_acc`  : translational accommodation coefficient ∈ [0,1]
- `mass_ic`    : particle mass (kg)

# Returns
3-component velocity vector (m/s) in the wall-local frame
(tang₁, tang₂, normal).
"""
function sample_wall_distribution(velo_square, wall_temp, trans_acc, mass_ic)
    # Helper: Box-Muller transform
    randn() = sqrt(-log(rand()))

    etra_old = 0.5 * mass_ic * velo_square

    velo_crad = randn()
    velo_cz   = randn()
    fak_d     = velo_crad^2 + velo_cz^2

    etra_new = etra_old + trans_acc * (BOLTZMANN * wall_temp * fak_d - etra_old)
    cmr      = sqrt(2.0 * etra_new / (mass_ic * fak_d))
    phi      = 2.0 * PI * rand()

    return SVector(
        cmr * velo_crad * cos(phi),   # tang1
        cmr * velo_crad * sin(phi),   # tang2
        cmr * velo_cz                  # normal
    )
end