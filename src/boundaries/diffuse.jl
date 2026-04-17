using Atomix


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