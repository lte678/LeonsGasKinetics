using Test
using LeonsGasKinetics
using Statistics

@testset "Simulation Test" begin
    @testset "BGK Couette" begin
        # Test shear flow between moving plates
        BASE_FOLDER = dirname(dirname(pathof(LeonsGasKinetics)))
        config_path = joinpath(BASE_FOLDER, "test", "configs", "couette.toml")
        volume_output, sim_config, mesh = run_simulation_from_config(config_path; enable_asserts=true)

        # Extract velocity profile along y-axis (assuming 2D or 3D with walls at y=0 and y=H)
        x = [2*(cell.barycenter[1] - 0.5) for cell in mesh.cells]
        x_sortidx = sortperm(x)
        x = x[x_sortidx]
        y_vel = volume_output["Total_VeloY"][x_sortidx]
 
        @test 13.5 < cov(x, y_vel) < 14.0
        @test all(abs.(volume_output["Total_VeloX"]) .< 5.0)
        @test all(abs.(volume_output["Total_VeloZ"]) .< 5.0)
        # The temperature profile is not exactly known
        @test all(volume_output["Total_TempTransX"] .< 320.0)
        @test all(volume_output["Total_TempTransX"] .> 250.0)
        @test all(abs.(volume_output["Total_NumberDensity"] .- 1.3e18) .< 1e17)
    end

    @testset "VR + BGK Couette" begin
        # Test shear flow between moving plates
        BASE_FOLDER = dirname(dirname(pathof(LeonsGasKinetics)))
        config_path = joinpath(BASE_FOLDER, "test", "configs", "couette_vrbgk.toml")
        volume_output, sim_config, mesh = run_simulation_from_config(config_path; enable_asserts=true)

        # Extract velocity profile along y-axis (assuming 2D or 3D with walls at y=0 and y=H)
        x = [2*(cell.barycenter[1] - 0.5) for cell in mesh.cells]
        x_sortidx = sortperm(x)
        x = x[x_sortidx]
        y_vel = volume_output["Total_VeloY"][x_sortidx]

        @test abs(cov(x, y_vel) - 0.136) < 0.002
        @test all(abs.(volume_output["Total_VeloX"]) .< 0.05)
        @test all(abs.(volume_output["Total_VeloZ"]) .< 0.05)
        @test all(volume_output["Total_TempTransX"] .< 280.1)
        @test all(volume_output["Total_TempTransX"] .> 279.9)
        @test all(abs.(volume_output["Total_NumberDensity"] .- 1.3e18) .< 0.02e18)
    end
end