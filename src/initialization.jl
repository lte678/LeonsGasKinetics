# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

function initialize_simulation!(sim::SimulationState, mesh::Mesh, sim_config::SimulationConfig, initializations)
    if length(sim.cell_part_count) < length(mesh.cells)
        sim.cell_part_count = zeros(length(mesh.cells))
    end
    
    for init_config in initializations
        species = init_config["species"]
        println("Inserting particles of species \"$species\"")
        initialize!(sim, mesh, sim_config, init_config)
    end

    if sim_config.asserts
        assert_particles_in_mesh(sim.particles, mesh)
    end
end

"""
    initialize!(sim::Simulation, init_dict::Dict)

Initialize particles in the simulation domain based on the provided initialization dictionary.

This function creates particles according to the specified density, temperature, and velocity
parameters. The number of particles in each cell is determined by the cell volume, density,
and macro particle factor (mpf). Particles are placed at random positions within each cell
and assigned velocities according to the specified velocity distribution.

# Arguments
- `sim::Simulation`: The simulation structure containing particles, mesh, and mpf
- `init_dict::Dict`: Initialization parameters with the following keys:
    - `density::Float64`: Number density of particles [m^-3]
    - `temperature::Float64`: Temperature [K]
    - `velocity::Vector{Float64}`: Bulk velocity vector [m/s]
    - `velocity_dist::String`: Velocity distribution type ("maxwell" for Maxwell-Boltzmann)

# Examples
```julia
init_params = Dict(
    "density" => 1.3e18,
    "temperature" => 280.0,
    "velocity" => [1.0, 0.0, 0.0],
    "velocity_dist" => "maxwell"
)
initialize!(sim, init_params)
```
"""
function initialize!(sim::SimulationState, mesh::Mesh, sim_config::SimulationConfig, init_dict::Dict)
    # Extract initialization parameters
    species = init_dict["species"]
    density = init_dict["density"]
    temperature = init_dict["temperature"]
    bulk_velocity = get(init_dict, "velocity", [0.0, 0.0, 0.0])
    velocity_dist = get(init_dict, "velocity_dist", "maxwell")
    
    # Process each cell in the mesh
    for (cell_i, cell) in enumerate(mesh.cells)
        # Calculate cell volume using numerical integration
        volume = cell_volume(cell)

        # Calculate number of real particles in this cell
        num_real_particles = density * volume
        
        # Calculate number of macro particles based on mpf
        num_macro_particles = round(Int, num_real_particles / sim_config.mpf)
        # Skip cells with no particles
        if num_macro_particles <= 0
            continue
        end
        
        # Generate random positions within the cell using rejection sampling
        positions = _intialize_positions(variant(cell), num_macro_particles)
        
        # Generate velocities based on distribution
        if velocity_dist == "maxwell"
            mass = sim_config.species[findfirst(s -> s.name == species, sim_config.species)].mass
            velocities = sample_maxwellian(temperature, bulk_velocity, mass, num_macro_particles)
        else
            throw(ArgumentError("Unsupported velocity distribution: $velocity_dist"))
        end
        
        # Insert particles into the simulation
        for i in 1:num_macro_particles
            insert_particle!(sim.particles, positions[i], velocities[i], cell_i)
            sim.cell_part_count[cell_i] += 1
        end
    end
    
    if sim_config.asserts
        assert_particles_in_mesh(sim.particles, mesh)
        assert_cell_part_count(sim.particles, sim.cell_part_count)
    end

    return sim
end


"""
Generate random positions within a hexahedral cell using rejection sampling.
"""
function _intialize_positions(cell::Hexahedron, n_particles::Int)
    positions = Vector{Vector{Float64}}(undef, n_particles)
    
    # Find maximum Jacobian determinant for rejection sampling
    max_jac_det = cell_max_jacobian(cell)
    
    particle_count = 0
    while particle_count < n_particles
        # Generate random point in reference space [-1, 1]^3
        xi = 2.0 * rand() - 1.0
        eta = 2.0 * rand() - 1.0
        zeta = 2.0 * rand() - 1.0
        
        # Calculate Jacobian determinant at this point
        J = cell_jacobian(cell, [xi, eta, zeta])
        jac_det = abs(det(J))
        
        # Rejection sampling based on Jacobian determinant
        if rand() * max_jac_det <= jac_det
            # Transform to global coordinates
            pos = cell_to_glob(cell, [xi, eta, zeta])
            particle_count += 1
            positions[particle_count] = [pos[1], pos[2], pos[3]]
        end
    end
    
    return positions
end
