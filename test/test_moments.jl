# Leon's Gas Kinetics
# Copyright (C) 2026  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

using LeonsGasKinetics
using Test
using Statistics
using StaticArrays
using Random

# Import necessary functions and constants from LeonsGasKinetics
import LeonsGasKinetics: sample_maxwellian, BOLTZMANN, ParticleData, 
                          MomentAccumulator, accumulate_moments!, 
                          calc_flow_properties, add_moment!, clear_moments!

# AI DISCLAIMER: These tests are AI generated.
# I have checked them and they look quite good. They actually test what they are supposed to.

# Set random seed for reproducibility
Random.seed!(43)

@testset "MomentAccumulator Tests" begin
    # Common test parameters
    mass = 4.65e-26  # Mass of argon atom (kg)
    cell_volume = 1.0e-15  # m^3
    mpf = 1.0e10  # Macro particle factor
    
    # Create a minimal simulation config
    species = [LeonsGasKinetics.SpeciesConfig(mass, 3.66e-10, 273.0, 0.81, "Argon")]
    vrbgk_config = LeonsGasKinetics.VRBGKConfig(false, 0.0, 0.0, false, 0.9)
    sim_config = LeonsGasKinetics.SimulationConfig(
        species,
        Vector{LeonsGasKinetics.Boundary}(),
        (pdata, samples, config, flow_variables, dt) -> (),
        mpf,
        1.0,  # t_end
        1e-6,  # dt
        1.0,  # sample_fraction
        1,  # output_interval
        1.0,  # report_interval
        "test",
        "test.msh",
        "output",
        true,  # silent
        false,  # asserts
        vrbgk_config,
        :continuous,
        3  # degrees_of_freedom
    )

    @testset "Test 1: Maxwell-Boltzmann Distribution" begin
        # Test with known bulk velocity and temperature
        bulk_velocity = SVector(100.0, 200.0, 300.0)  # m/s
        temperature = 300.0  # K
        n_particles = 10000
        
        # Generate Maxwell-Boltzmann distributed velocities
        velocities = sample_maxwellian(temperature, bulk_velocity, mass, n_particles)
        
        # Create particle data
        particles = ParticleData(n_particles)
        for i in 1:n_particles
            particles.vel[i] = velocities[i]
            particles.pos[i] = SVector(0.0, 0.0, 0.0)
            particles.cell[i] = 1
        end
        
        # Create moment accumulator
        moments = [MomentAccumulator()]
        
        # Accumulate moments
        accumulate_moments!(moments, particles)
        
        # Calculate flow properties
        flow_props = calc_flow_properties(moments[1], sim_config, cell_volume)
        
        # Test velocity - should match bulk velocity within statistical error
        # Standard error of mean = sigma / sqrt(n)
        # For MB distribution, sigma = sqrt(kB * T / m)
        sigma = sqrt(BOLTZMANN * temperature / mass)
        velocity_error = 3 * sigma / sqrt(n_particles)  # 3-sigma confidence
        
        @test isapprox(flow_props.velocity, bulk_velocity, atol=velocity_error)
        
        # Test temperature - should match input temperature within statistical error
        # Standard error of variance = sqrt(2) * sigma^2 / sqrt(n)
        temp_error = 3 * sqrt(2) * temperature / sqrt(n_particles)  # 3-sigma confidence
        
        @test isapprox(flow_props.temperature, SVector(temperature, temperature, temperature), rtol=temp_error/temperature)
        
        # Test mean temperature
        @test isapprox(flow_props.mean_temperature, temperature, rtol=temp_error/temperature)
        
        # Test density
        expected_density = n_particles * mpf / cell_volume
        @test isapprox(flow_props.density, expected_density, rtol=1e-10)
        
        # Test particle count
        @test isapprox(flow_props.sim_particle_count, n_particles, rtol=1e-10)
    end

    @testset "Test 2: Linear Interpolation of Two MB Distributions" begin
        # Distribution 1
        bulk_velocity_1 = SVector(50.0, 100.0, 150.0)  # m/s
        temperature_1 = 250.0  # K
        n_particles_1 = 5000
        
        # Distribution 2
        bulk_velocity_2 = SVector(150.0, 200.0, 250.0)  # m/s
        temperature_2 = 350.0  # K
        n_particles_2 = 5000
        
        total_particles = n_particles_1 + n_particles_2
        weight_1 = n_particles_1 / total_particles
        weight_2 = n_particles_2 / total_particles
        
        # Generate velocities from both distributions
        velocities_1 = sample_maxwellian(temperature_1, bulk_velocity_1, mass, n_particles_1)
        velocities_2 = sample_maxwellian(temperature_2, bulk_velocity_2, mass, n_particles_2)
        
        # Create particle data
        particles = ParticleData(total_particles)
        for i in 1:n_particles_1
            particles.vel[i] = velocities_1[i]
            particles.pos[i] = SVector(0.0, 0.0, 0.0)
            particles.cell[i] = 1
        end
        for i in 1:n_particles_2
            idx = n_particles_1 + i
            particles.vel[idx] = velocities_2[i]
            particles.pos[idx] = SVector(0.0, 0.0, 0.0)
            particles.cell[idx] = 1
        end
        
        # Create moment accumulator
        moments = [MomentAccumulator()]
        
        # Accumulate moments
        accumulate_moments!(moments, particles)
        
        # Calculate flow properties
        flow_props = calc_flow_properties(moments[1], sim_config, cell_volume)
        
        # Expected velocity is weighted mean
        expected_velocity = weight_1 * bulk_velocity_1 + weight_2 * bulk_velocity_2
        
        # Expected temperature calculation
        # For each axis: T = (m/kB) * (<v^2> - <v>^2)
        # <v^2> = w1 * <v1^2> + w2 * <v2^2>
        # <v> = w1 * <v1> + w2 * <v2>
        # For MB: <v_i^2> = u_i^2 + kB*T/m
        expected_temp_array = [0.0, 0.0, 0.0]
        for axis in 1:3
            mean_v_sq_1 = bulk_velocity_1[axis]^2 + BOLTZMANN * temperature_1 / mass
            mean_v_sq_2 = bulk_velocity_2[axis]^2 + BOLTZMANN * temperature_2 / mass
            mean_v_sq = weight_1 * mean_v_sq_1 + weight_2 * mean_v_sq_2
            mean_v = expected_velocity[axis]
            expected_temp_array[axis] = (mass / BOLTZMANN) * (mean_v_sq - mean_v^2)
        end
        expected_temperature = SVector(expected_temp_array[1], expected_temp_array[2], expected_temp_array[3])
        
        # We use 5-sigma confidence for this test, since there is probably some additional variance introduced 
        # Statistical error estimation
        sigma_1 = sqrt(BOLTZMANN * temperature_1 / mass)
        sigma_2 = sqrt(BOLTZMANN * temperature_2 / mass)
        velocity_error = 5 * sqrt(weight_1^2 * sigma_1^2 + weight_2^2 * sigma_2^2) / sqrt(total_particles)
        
        # Test velocity
        @test isapprox(flow_props.velocity, expected_velocity, atol=velocity_error)
        
        # Test temperature
        temp_error = 5 * sqrt(2) * sqrt(weight_1^2 * temperature_1^2 + weight_2^2 * temperature_2^2) / sqrt(total_particles)
        @test isapprox(flow_props.temperature, expected_temperature, rtol=temp_error/expected_temperature[1])
        
        # Test density
        expected_density = total_particles * mpf / cell_volume
        @test isapprox(flow_props.density, expected_density, rtol=1e-10)
    end

    @testset "Test 3: Anisotropic Temperature Distribution" begin
        # Create distribution with different temperatures in each axis
        bulk_velocity = SVector(0.0, 0.0, 0.0)  # m/s
        temperatures = SVector(200.0, 300.0, 400.0)  # K - different for each axis
        n_particles = 10000
        
        # Generate velocities with anisotropic temperature
        velocities = Vector{SVector{3, Float64}}(undef, n_particles)
        for i in 1:n_particles
            v_x = bulk_velocity[1] + sqrt(BOLTZMANN * temperatures[1] / mass) * randn()
            v_y = bulk_velocity[2] + sqrt(BOLTZMANN * temperatures[2] / mass) * randn()
            v_z = bulk_velocity[3] + sqrt(BOLTZMANN * temperatures[3] / mass) * randn()
            velocities[i] = SVector(v_x, v_y, v_z)
        end
        
        # Create particle data
        particles = ParticleData(n_particles)
        for i in 1:n_particles
            particles.vel[i] = velocities[i]
            particles.pos[i] = SVector(0.0, 0.0, 0.0)
            particles.cell[i] = 1
        end
        
        # Create moment accumulator
        moments = [MomentAccumulator()]
        
        # Accumulate moments
        accumulate_moments!(moments, particles)
        
        # Calculate flow properties
        flow_props = calc_flow_properties(moments[1], sim_config, cell_volume)
        
        # Test velocity - should be zero (within statistical error)
        # Calculate statistical error for each axis
        for axis in 1:3
            sigma = sqrt(BOLTZMANN * temperatures[axis] / mass)
            velocity_error = 3 * sigma / sqrt(n_particles)
            @test isapprox(flow_props.velocity[axis], 0.0, atol=velocity_error)
        end
        
        # Test anisotropic temperatures
        for axis in 1:3
            sigma = sqrt(BOLTZMANN * temperatures[axis] / mass)
            temp_error = 3 * sqrt(2) * temperatures[axis] / sqrt(n_particles)
            @test isapprox(flow_props.temperature[axis], temperatures[axis], rtol=temp_error/temperatures[axis])
        end
        
        # Test mean temperature
        expected_mean_temp = sum(temperatures) / 3.0
        @test isapprox(flow_props.mean_temperature, expected_mean_temp, rtol=0.05)
    end

    @testset "Test 4: Multiple Cells" begin
        # Test moment accumulation across multiple cells
        n_cells = 3
        n_particles_per_cell = 2000
        
        particles = ParticleData(n_cells * n_particles_per_cell)
        moments = [MomentAccumulator() for _ in 1:n_cells]
        
        # Different bulk velocities for each cell
        bulk_velocities = [
            SVector(100.0, 0.0, 0.0),
            SVector(0.0, 100.0, 0.0),
            SVector(0.0, 0.0, 100.0)
        ]
        temperature = 300.0
        
        # Generate particles for each cell
        for cell in 1:n_cells
            start_idx = (cell - 1) * n_particles_per_cell + 1
            velocities = sample_maxwellian(temperature, bulk_velocities[cell], mass, n_particles_per_cell)
            
            for i in 1:n_particles_per_cell
                idx = start_idx + i - 1
                particles.vel[idx] = velocities[i]
                particles.pos[idx] = SVector(0.0, 0.0, 0.0)
                particles.cell[idx] = cell
            end
        end
        
        # Accumulate moments
        accumulate_moments!(moments, particles)
        
        # Check each cell
        sigma = sqrt(BOLTZMANN * temperature / mass)
        velocity_error = 3 * sigma / sqrt(n_particles_per_cell)
        temp_error = 3 * sqrt(2) * temperature / sqrt(n_particles_per_cell)
        
        for cell in 1:n_cells
            flow_props = calc_flow_properties(moments[cell], sim_config, cell_volume)
            
            @test isapprox(flow_props.velocity, bulk_velocities[cell], atol=velocity_error)
            
            @test isapprox(flow_props.temperature, SVector(temperature, temperature, temperature), rtol=temp_error/temperature)
        end
    end

    @testset "Test 5: Empty MomentAccumulator" begin
        # Test behavior with no particles
        moments = [MomentAccumulator()]
        particles = ParticleData(0)
        
        accumulate_moments!(moments, particles)
        flow_props = calc_flow_properties(moments[1], sim_config, cell_volume)
        
        # Should return default FlowProperties
        @test flow_props.velocity == SVector(0.0, 0.0, 0.0)
        @test flow_props.temperature == SVector(0.0, 0.0, 0.0)
        @test flow_props.mean_temperature == 0.0
        @test flow_props.density == 0.0
        @test flow_props.sim_particle_count == 0.0
    end

    @testset "Test 6: Moment Accumulation Over Multiple Samples" begin
        # Test that accumulating moments over multiple samples works correctly
        n_samples = 5
        n_particles_per_sample = 1000
        bulk_velocity = SVector(50.0, 100.0, 150.0)
        temperature = 300.0
        
        moments = [MomentAccumulator()]
        
        for sample in 1:n_samples
            velocities = sample_maxwellian(temperature, bulk_velocity, mass, n_particles_per_sample)
            
            particles = ParticleData(n_particles_per_sample)
            for i in 1:n_particles_per_sample
                particles.vel[i] = velocities[i]
                particles.pos[i] = SVector(0.0, 0.0, 0.0)
                particles.cell[i] = 1
            end
            
            accumulate_moments!(moments, particles)
        end
        
        flow_props = calc_flow_properties(moments[1], sim_config, cell_volume)
        
        # Should still match expected values
        total_particles = n_samples * n_particles_per_sample
        sigma = sqrt(BOLTZMANN * temperature / mass)
        velocity_error = 3 * sigma / sqrt(total_particles)
        temp_error = 3 * sqrt(2) * temperature / sqrt(total_particles)
        
        @test isapprox(flow_props.velocity, bulk_velocity, atol=velocity_error)
        
        @test isapprox(flow_props.temperature, SVector(temperature, temperature, temperature), rtol=temp_error/temperature)
        
        # Check that samples counter is correct
        @test moments[1].samples == n_samples
    end

    @testset "Test 7: add_moment! Function" begin
        # Test the add_moment! function
        n_particles_1 = 1000
        n_particles_2 = 1500
        
        # Create two sets of particles
        velocities_1 = sample_maxwellian(300.0, SVector(100.0, 0.0, 0.0), mass, n_particles_1)
        velocities_2 = sample_maxwellian(300.0, SVector(0.0, 100.0, 0.0), mass, n_particles_2)
        
        # Accumulate separately
        moments_1 = [MomentAccumulator()]
        particles_1 = ParticleData(n_particles_1)
        for i in 1:n_particles_1
            particles_1.vel[i] = velocities_1[i]
            particles_1.pos[i] = SVector(0.0, 0.0, 0.0)
            particles_1.cell[i] = 1
        end
        accumulate_moments!(moments_1, particles_1)
        
        moments_2 = [MomentAccumulator()]
        particles_2 = ParticleData(n_particles_2)
        for i in 1:n_particles_2
            particles_2.vel[i] = velocities_2[i]
            particles_2.pos[i] = SVector(0.0, 0.0, 0.0)
            particles_2.cell[i] = 1
        end
        accumulate_moments!(moments_2, particles_2)
        
        # Add moments together
        add_moment!(moments_1[1], moments_2[1])
        
        # Create combined particles and accumulate
        total_particles = n_particles_1 + n_particles_2
        particles_combined = ParticleData(total_particles)
        for i in 1:n_particles_1
            particles_combined.vel[i] = velocities_1[i]
            particles_combined.pos[i] = SVector(0.0, 0.0, 0.0)
            particles_combined.cell[i] = 1
        end
        for i in 1:n_particles_2
            idx = n_particles_1 + i
            particles_combined.vel[idx] = velocities_2[i]
            particles_combined.pos[idx] = SVector(0.0, 0.0, 0.0)
            particles_combined.cell[idx] = 1
        end
        
        moments_combined = [MomentAccumulator()]
        accumulate_moments!(moments_combined, particles_combined)
        
        # Results should be the same (except samples which add up)
        @test moments_1[1].c_i ≈ moments_combined[1].c_i
        @test moments_1[1].c_ii ≈ moments_combined[1].c_ii
        @test moments_1[1].count ≈ moments_combined[1].count
        # samples add up when using add_moment! (2 separate accumulations vs 1 combined)
        @test moments_1[1].samples == 2
        @test moments_combined[1].samples == 1
    end

    @testset "Test 8: clear_moments! Function" begin
        # Test the clear_moments! function
        n_particles = 1000
        velocities = sample_maxwellian(300.0, SVector(100.0, 0.0, 0.0), mass, n_particles)
        
        particles = ParticleData(n_particles)
        for i in 1:n_particles
            particles.vel[i] = velocities[i]
            particles.pos[i] = SVector(0.0, 0.0, 0.0)
            particles.cell[i] = 1
        end
        
        moments = [MomentAccumulator()]
        accumulate_moments!(moments, particles)
        
        # Verify moments are not zero
        @test moments[1].count > 0
        @test moments[1].samples > 0
        @test any(moments[1].c_i .!= 0.0)
        @test any(moments[1].c_ii .!= 0.0)
        
        # Clear moments
        clear_moments!(moments)
        
        # Verify moments are reset
        @test moments[1].c_i == SVector(0.0, 0.0, 0.0)
        @test moments[1].c_ii == SVector(0.0, 0.0, 0.0)
        @test moments[1].count == 0.0
        @test moments[1].samples == 0
    end

    @testset "Test 9: High Temperature Distribution" begin
        # Test with high temperature to ensure numerical stability
        bulk_velocity = SVector(0.0, 0.0, 0.0)
        temperature = 10000.0  # K - very high temperature
        n_particles = 5000
        
        velocities = sample_maxwellian(temperature, bulk_velocity, mass, n_particles)
        
        particles = ParticleData(n_particles)
        for i in 1:n_particles
            particles.vel[i] = velocities[i]
            particles.pos[i] = SVector(0.0, 0.0, 0.0)
            particles.cell[i] = 1
        end
        
        moments = [MomentAccumulator()]
        accumulate_moments!(moments, particles)
        flow_props = calc_flow_properties(moments[1], sim_config, cell_volume)
        
        # Temperature should be correct even at high values
        temp_error = 3 * sqrt(2) * temperature / sqrt(n_particles)
        @test isapprox(flow_props.temperature, SVector(temperature, temperature, temperature), rtol=temp_error/temperature)
    end
end