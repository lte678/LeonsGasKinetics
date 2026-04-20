# Leon's Gas Kinetics
# Copyright (C) 2026  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

function initialize_simulation!(sim::SimulationState, mesh::Mesh, sim_config::SimulationConfig, initializations)
    # Initialize per-cell reference temperature and density to global values
    if sim_config.vrbgk.enabled
        fill!(sim.cells.features.vrbgk_ref_temperature, sim_config.vrbgk.ref_temperature)
        fill!(sim.cells.features.vrbgk_ref_density, sim_config.vrbgk.ref_density)
    end

    # Insert particles
    for init_config in initializations
        species = init_config["species"]

        if !sim_config.silent
            println("Inserting particles of species \"$species\"")
        end
        initial_particles!(sim, mesh, sim_config, init_config)
    end

    if sim_config.asserts
        assert_particles_in_mesh(sim.particles, mesh)
    end
end

"""
Initialize particles in the simulation domain based on the provided initialization dictionary.

This function creates particles according to the specified density, temperature, and velocity
parameters. The number of particles in each cell is determined by the cell volume, density,
and macro particle factor (mpf). Particles are placed at random positions within each cell
and assigned velocities according to the specified velocity distribution.
"""
function initial_particles!(sim::SimulationState, mesh::Mesh, sim_config::SimulationConfig, init_dict::Dict)
    # Extract initialization parameters
    species = init_dict["species"]
    density = init_dict["density"]
    temperature = init_dict["temperature"]
    bulk_velocity = get(init_dict, "velocity", [0.0, 0.0, 0.0])
    velocity_dist = get(init_dict, "velocity_dist", "maxwell")
    
    total_inserted = 0

    # Process each cell in the mesh
    for (cell_i, cell) in enumerate(mesh.cells)
        # Calculate cell volume using numerical integration
        volume = cell_volume(cell)

        # Calculate number of real particles in this cell
        num_real_particles = density * volume
        
        # Calculate number of macro particles based on mpf
        num_particles = round(Int, num_real_particles / sim_config.mpf)

        # Skip cells with no particles
        if num_particles <= 0
            continue
        end
        
        # Generate random positions within the cell using rejection sampling
        positions = _intialize_positions(variant(cell), num_particles, sim_config.degrees_of_freedom)
        
        # Generate velocities based on distribution
        if velocity_dist == "maxwell"
            mass = sim_config.species[findfirst(s -> s.name == species, sim_config.species)].mass
            velocities = sample_maxwellian(temperature, bulk_velocity, mass, num_particles)
        else
            throw(ArgumentError("Unsupported velocity distribution: $velocity_dist"))
        end
        
        # Insert particles into the simulation
        for i in 1:num_particles
            feature_fields = []
            if sim_config.vrbgk.enabled
                # Ratio of probabilities. Used for importance sampling. See VRBGK papers.
                vr_weight = (sim_config.vrbgk.ref_density .* maxwellian(sim.cells.features.vrbgk_ref_temperature[cell_i], [0.0, 0.0, 0.0], mass, velocities[i])) ./
                            (                     density .* maxwellian(temperature                        , bulk_velocity  , mass, velocities[i]))
                push!(feature_fields, :vr_weight => vr_weight)
            end

            p = SingleParticle(
                positions[i],
                velocities[i],
                cell_i,
                (; feature_fields...)  # Convert to named tuple
            )
            insert_particle!(sim.particles, p)
            sim.cells.part_count[cell_i] += 1
        end

        total_inserted += num_particles
    end
    
    if sim_config.asserts
        assert_particles_in_mesh(sim.particles, mesh)
        assert_cell_part_count(sim.particles, sim.cells.part_count)
    end

    if !sim_config.silent
        println("Inserted $total_inserted particles.")
    end

    return sim
end


"""
Generate random positions within a hexahedral cell using rejection sampling.
"""
function _intialize_positions(cell::Hexahedron, n_particles::Int, degrees_of_freedom)
    positions = Vector{Vector{Float64}}(undef, n_particles)
    
    # Find maximum Jacobian determinant for rejection sampling
    max_jac_det = cell_max_jacobian(cell)
    
    particle_count = 0
    while particle_count < n_particles
        # Generate random point in reference space [-1, 1]^3
        xi = 2.0 * rand() - 1.0
        eta = 2.0 * rand() - 1.0
        if degrees_of_freedom == 2
            zeta = 0.0
        else
            zeta = 2.0 * rand() - 1.0
        end

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
