using StaticArrays

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


@testset "Hexahedron Cell Containment" begin
    @testset "Basic internal/external" begin
        cube = TestHexahedron()
        
        @test cell_contains(cube, SVector(0.5, 0.5, 0.5)) == true
        @test cell_contains(cube, SVector(1.5, 0.5, 0.5)) == false
        @test cell_contains(cube, SVector(-0.1, 0.5, 0.5)) == false
    end

    @testset "Face boundaries and epsilon offsets" begin
        cube = TestHexahedron()
        
        # Test directions: +x, -x, +y, -y, +z, -z
        directions = [
            SVector(1.0, 0.0, 0.0), SVector(-1.0, 0.0, 0.0),
            SVector(0.0, 1.0, 0.0), SVector(0.0, -1.0, 0.0),
            SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -1.0)
        ]
        
        # Center of each face
        face_centers = [
            SVector(1.0, 0.5, 0.5), SVector(0.0, 0.5, 0.5),
            SVector(0.5, 1.0, 0.5), SVector(0.5, 0.0, 0.5),
            SVector(0.5, 0.5, 1.0), SVector(0.5, 0.5, 0.0)
        ]

        for i in 1:6
            p_on = face_centers[i]
            dir = directions[i]
            
            # Precisely on the surface
            @test cell_contains(cube, p_on)
            
            # Inside by 2*epsilon
            p_in = p_on - 2 * dir .* eps(dot(dir, p_on))
            @test cell_contains(cube, p_in)
            
            # Outside by 2*epsilon
            p_out = p_on + 2 * dir .* eps(dot(dir, p_on))
            @test cell_contains(cube, p_out)
        end
    end

    @testset "Edge and corner" begin
        cube = TestHexahedron()

        # Edge test: Edge is defined along [0, 1, 1] (X-axis at Y=1, Z=1)
        # Test just outside the edge in the diagonal direction
        edge_point_out = SVector(0.5, 1.0 + 1e-12, 1.0 + 1e-12)
        @test cell_contains(cube, edge_point_out) == false

        # Corner test: Corner at [1, 1, 1]
        # Test slightly outside the corner in all three dimensions
        corner_point_out = SVector(1.0, 1.0, 1.0) + fill(1e-12, 3)
        @test cell_contains(cube, corner_point_out) == false
        
        # Test corner exactly
        @test cell_contains(cube, SVector(1.0, 1.0, 1.0)) == true
    end
end