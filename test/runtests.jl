using LeonsGasKinetics

using Test
using Statistics
using StaticArrays
using LinearAlgebra


@testset "Tracking Tests" begin
    # 1. Define a unit cube (0,0,0) to (1,1,1)
    # Normals must point INWARD according to function requirements
    normals = SVector{6, SVector{3, Float64}}(
        SVector(1.0, 0.0, 0.0),  # Left face (x=0)
        SVector(-1.0, 0.0, 0.0), # Right face (x=1)
        SVector(0.0, 1.0, 0.0),  # Bottom face (y=0)
        SVector(0.0, -1.0, 0.0), # Top face (y=1)
        SVector(0.0, 0.0, 1.0),  # Back face (z=0)
        SVector(0.0, 0.0, -1.0)  # Front face (z=1)
    )

    # Points on each face
    vertices = SVector{6, SVector{3, Float64}}(
        SVector(0.0, 0.5, 0.5), # Left
        SVector(1.0, 0.5, 0.5), # Right
        SVector(0.5, 0.0, 0.5), # Bottom
        SVector(0.5, 1.0, 0.5), # Top
        SVector(0.5, 0.5, 0.0), # Back
        SVector(0.5, 0.5, 1.0)  # Front
    )

    @testset "Single collision" begin
        origin = SVector(0.5, 0.5, 0.5)
        dir = SVector(1.0, 0.0, 0.0)
        
        idx, dist = find_exit_face(origin, dir, normals, vertices)
        
        # Should hit the Right face (index 2)
        @test idx == 2
        # Distance from 0.5 to 1.0 is 0.5
        @test dist ≈ 0.5

        idx, dist = find_exit_face(origin, -dir, normals, vertices)
        
        # Should hit Left face (index 1)
        @test idx == 1
        @test dist ≈ 0.5
    end

    @testset "Double collision" begin
        origin = SVector(0.5, 0.5, 0.5)
        # Moving towards top-right corner
        dir = normalize(SVector(1.0, 1.0, 0.0))
        
        idx, dist = find_exit_face(origin, dir, normals, vertices)
        
        # It should hit either index 2 (Right) or 4 (Top) 
        # because the distance is the same
        @test idx ∈ (2, 4)
        @test dist ≈ 0.5 / dir[1]
    end

    @testset "Test magnitude" begin
        origin = SVector(0.5, 0.5, 0.5)
        dir = SVector(0.0, 0.0, 100.0)
        
        idx, dist = find_exit_face(origin, dir, normals, vertices)
        
        # The distance needs to express itself in terms of the length of the direction
        @test idx == 6
        @test dist ≈ 0.5 / 100.0
    end

    @testset "Near boundary" begin
        # Start very close to the right wall, moving right
        origin = SVector(1.0 - eps(1.0), 0.5, 0.5)
        dir = SVector(1.0, 0.0, 0.0)
        
        idx, dist = find_exit_face(origin, dir, normals, vertices)
        
        @test idx == 2
        @test abs(dist) < 1e-10
    end

    @testset "In boundary" begin
        # Start very close to the right wall, moving right
        origin = SVector(1.0, 0.5, 0.5)
        dir = SVector(1.0, 0.0, 0.0)
        
        idx, dist = find_exit_face(origin, dir, normals, vertices)
        
        @test idx == 2
        @test abs(dist) < 1e-10
    end

    @testset "Past boundary" begin
        # Start past the right wall, but moving right.
        # This is of course undesirable, but we need a hit anyway
        origin = SVector(1.0 + eps(1.0), 0.5, 0.5)
        dir = SVector(1.0, 0.0, 0.0)
        
        idx, dist = find_exit_face(origin, dir, normals, vertices)
        
        @test idx == 2
        @test abs(dist) < 1e-10
    end
end

@testset "Simulation Test" begin
    @testset "BGK Couette" begin
        # Test shear flow between moving plates
        BASE_FOLDER = dirname(dirname(pathof(LeonsGasKinetics)))
        config_path = joinpath(BASE_FOLDER, "test", "configs", "couette.toml")
        volume_output, sim_config, mesh = run_simulation_from_config(config_path, "output"; enable_asserts=true)

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
        volume_output, sim_config, mesh = run_simulation_from_config(config_path, "output"; enable_asserts=true)

        # Extract velocity profile along y-axis (assuming 2D or 3D with walls at y=0 and y=H)
        x = [2*(cell.barycenter[1] - 0.5) for cell in mesh.cells]
        x_sortidx = sortperm(x)
        x = x[x_sortidx]
        y_vel = volume_output["Total_VeloY"][x_sortidx]

        @test abs(cov(x, y_vel) - 0.136) < 0.005
        @test all(abs.(volume_output["Total_VeloX"]) .< 0.05)
        @test all(abs.(volume_output["Total_VeloZ"]) .< 0.05)
        @test all(volume_output["Total_TempTransX"] .< 280.1)
        @test all(volume_output["Total_TempTransX"] .> 279.9)
        @test all(abs.(volume_output["Total_NumberDensity"] .- 1.3e18) .< 0.02e18)
    end
end