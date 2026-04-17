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