function handle_boundary(
    particle::SingleParticle,
    side::BoundarySide,
    side_idx,
    species::SpeciesConfig,
    boundary::ReflectiveBoundary,
    config::SimulationConfig
) :: SingleParticle
    post_particle = @set particle.vel -= 2*side.normal*dot(side.normal, particle.vel)
    return post_particle
end