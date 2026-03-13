using Test
using LeonsGasKinetics

@testset "Simulation Test" begin
    @testset "BGK Couette" begin
        # Test shear flow between moving plates
        BASE_FOLDER = dirname(dirname(pathof(LeonsGasKinetics)))
        WALL_VEL = 36.0
        
        config_path = joinpath(BASE_FOLDER, "test", "configs", "couette.toml")

        volume_output, sim_config, mesh = run_simulation_from_config(config_path; enable_asserts=true)

        # Extract velocity profile along y-axis (assuming 2D or 3D with walls at y=0 and y=H)
        x = [2*(cell.barycenter[1] - 0.5) for cell in mesh.cells]
        x_sortidx = sortperm(x)
        x = x[x_sortidx]
        y_vel = volume_output["Total_VeloY"][x_sortidx]
        y_vel_expected = WALL_VEL .* x
        @test all(abs.(y_vel_expected[3:end-2] .- y_vel[3:end-2]) .< 7.5e-2 * WALL_VEL)
        @test all(abs.(volume_output["Total_VeloX"]) .< 1.0)
        @test all(abs.(volume_output["Total_VeloZ"]) .< 2.5)
        # The temperature profile is not exactly known
        @test all(volume_output["Total_TempTransX"] .< 320.0)
        @test all(volume_output["Total_TempTransX"] .> 250.0)
        @test all(abs.(volume_output["Total_NumberDensity"] .- 1.3e18) .< 1e17)
    end
end