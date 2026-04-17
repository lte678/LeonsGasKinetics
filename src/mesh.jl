# Leon's Gas Kinetics
# Copyright (C) 2026  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

using LinearAlgebra
using StaticArrays
using HDF5
using Printf
using LightSumTypes
using Atomix


Vertex = SVector{3, Float64}

mutable struct BoundarySide
    bc_index :: Int16
    normal :: SVector{3, Float64}
    face_vertices :: SVector{4, SVector{3, Float64}}
    area :: Float64
    adjacent_cell :: UInt32
    # TODO: Only include this if VRBGK is enabled.
    @atomic vrbgk_incident_sum :: Float64
    @atomic vrbgk_incident_count :: UInt32
end

struct Hexahedron
    vertices :: SVector{8, Vertex}
    normals :: SVector{6, SVector{3, Float64}}
    # Contains one vertex per face. A arbitrary vertex on the face 
    face_origins :: SVector{6, SVector{3, Float64}}
    barycenter :: SVector{3, Float64}
    # Boundary conditions (zero if none)
    bc_side_idx :: SVector{6, UInt32}
    # Neighbour cell for side.
    neighbours :: SVector{6, UInt32}
    volume :: Float64
end

@sumtype Cell(Hexahedron)

struct Mesh
    cells :: Vector{Cell}
    bc_names :: Vector{String}
    bc_sides :: Vector{BoundarySide}
end


function get_sides(vertices :: AbstractVector{SVector{3, Float64}})
    SVector(
        SVector(vertices[1], vertices[4], vertices[3], vertices[2]),
        SVector(vertices[1], vertices[2], vertices[6], vertices[5]),
        SVector(vertices[2], vertices[3], vertices[7], vertices[6]),
        SVector(vertices[3], vertices[4], vertices[8], vertices[7]),
        SVector(vertices[1], vertices[5], vertices[8], vertices[4]),
        SVector(vertices[5], vertices[6], vertices[7], vertices[8])
    )
end

get_sides(cell :: Hexahedron) = get_sides(cell.vertices)
get_sides(cell :: Cell) = get_sides(variant(cell))


function _plane_side(r, p1, p2, p3)
    normal = cross(p1 - p2, p3 - p2)
    return dot(normal, r - p2)
end


function cell_contains(cell :: Hexahedron, r)
    faces = get_sides(cell)
    if _plane_side(r, faces[5][1], faces[5][2], faces[5][3]) < -COLLISION_TOL ||
       _plane_side(r, faces[3][1], faces[3][2], faces[3][3]) < -COLLISION_TOL
       return false
    end
    if _plane_side(r, faces[1][1], faces[1][2], faces[1][3]) < -COLLISION_TOL ||
       _plane_side(r, faces[6][1], faces[6][2], faces[6][3]) < -COLLISION_TOL
       return false
    end
    if _plane_side(r, faces[2][1], faces[2][2], faces[2][3]) < -COLLISION_TOL ||
       _plane_side(r, faces[4][1], faces[4][2], faces[4][3]) < -COLLISION_TOL
       return false
    end

    return true
end

cell_contains(cell :: Cell, r) = cell_contains(variant(cell), r)


function side_normal(side)
    n = cross(side[1] - side[2], side[3] - side[2])
    return n / norm(n)
end


function assert_normal_direction(cell :: Cell)
    for side in get_sides(cell)
        if dot(side_normal(side), cell.barycenter - side[1]) < 0.0
            error("Side normal of cell at $(cell.barycenter) is pointing outwards!")
        end
    end
end


# Transformation from local xi-cell space to global euler xyz space.
function cell_to_glob(cell :: Hexahedron, xi)
    # Vertices are in CGNS ordering.
    v = cell.vertices
    r = (1 - xi[1])*(1 - xi[2])*(1 - xi[3]) .* v[1] +
        (1 + xi[1])*(1 - xi[2])*(1 - xi[3]) .* v[2] +
        (1 + xi[1])*(1 + xi[2])*(1 - xi[3]) .* v[3] +
        (1 - xi[1])*(1 + xi[2])*(1 - xi[3]) .* v[4] +
        (1 - xi[1])*(1 - xi[2])*(1 + xi[3]) .* v[5] +
        (1 + xi[1])*(1 - xi[2])*(1 + xi[3]) .* v[6] +
        (1 + xi[1])*(1 + xi[2])*(1 + xi[3]) .* v[7] +
        (1 - xi[1])*(1 + xi[2])*(1 + xi[3]) .* v[8]
    return r / 8
end

cell_to_glob(cell :: Cell, xi) = cell_to_glob(variant(cell), xi)


function cell_jacobian(cell :: Hexahedron, xi)
    # Vertices are in CGNS ordering.
    v = cell.vertices
    # Column 1: dr / dxi[1]
    J_1 = ((1 - xi[2])*(1 - xi[3]) .* (v[2] - v[1]) +
           (1 + xi[2])*(1 - xi[3]) .* (v[3] - v[4]) +
           (1 - xi[2])*(1 + xi[3]) .* (v[6] - v[5]) +
           (1 + xi[2])*(1 + xi[3]) .* (v[7] - v[8])) / 8
    # Column 2: dr / dxi[2]
    J_2 = ((1 - xi[1])*(1 - xi[3]) .* (v[4] - v[1]) +
           (1 + xi[1])*(1 - xi[3]) .* (v[3] - v[2]) +
           (1 - xi[1])*(1 + xi[3]) .* (v[8] - v[5]) +
           (1 + xi[1])*(1 + xi[3]) .* (v[7] - v[6])) / 8
    # Column 3: dr / dxi[3]
    J_3 = ((1 - xi[1])*(1 - xi[2]) .* (v[5] - v[1]) +
           (1 + xi[1])*(1 - xi[2]) .* (v[6] - v[2]) +
           (1 + xi[1])*(1 + xi[2]) .* (v[7] - v[3]) +
           (1 - xi[1])*(1 + xi[2]) .* (v[8] - v[4])) / 8
    
    return [J_1 J_2 J_3]
end

cell_jacobian(cell :: Cell, xi) = cell_jacobian(variant(cell), xi)


"""
Find maximum Jacobian determinant in the cell for rejection sampling.
"""
function cell_max_jacobian(cell :: Hexahedron)
    max_jac_det = 0.0    
    for xi_1 in [-1, 1]
        for xi_2 in [-1, 1]
            for xi_3 in [-1, 1]
                J = cell_jacobian(cell, [xi_1, xi_2, xi_3])
                max(max_jac_det, abs(det(J)))
            end
        end
    end
    return max_jac_det
end

cell_max_jacobian(cell :: Cell) = cell_max_jacobian(variant(cell))


"""
Calculate the volume of a hexahedral cell using numerical integration.
"""
function cell_volume(cell :: Hexahedron)
    cell_volume = 0.0
    for xi_1 in [-1, 1]
        for xi_2 in [-1, 1]
            for xi_3 in [-1, 1]
                J = cell_jacobian(cell, [xi_1, xi_2, xi_3])
                cell_volume += abs(det(J))
            end
        end
    end
    return cell_volume
end

cell_volume(cell :: Cell) = cell_volume(variant(cell))

function face_quad_area(face::SVector{4, SVector{3, Float64}})
    v1, v2, v3, v4 = face[1], face[2], face[3], face[4]
    a1 = 0.5 * norm(cross(v2 - v1, v3 - v1))
    a2 = 0.5 * norm(cross(v3 - v1, v4 - v1))
    return a1 + a2
end


"""
Imports meshes of the HOPR format.
https://hopr.readthedocs.io/en/latest/userguide/meshformat.html
"""
function mesh_from_h5(path)
    # Open HDF5 file
    h5file = h5open(path, "r")
    
    ngeo    = read_attribute(h5file, "Ngeo")[1]
    nelems  = read_attribute(h5file, "nElems")[1]
    if ngeo != 1
        error("Cannot load curved geometry (Ngeo=$ngeo).")
    end
    @printf "Loading mesh with %d cells.\n" nelems

    # Read data arrays
    elem_info = read(h5file, "ElemInfo")
    # This is no longer available in PyHOPE
    # barycenters = read(h5file, "ElemBarycenters")
    node_coords = read(h5file, "NodeCoords")
    side_info = read(h5file, "SideInfo")
    
    # Initialize mesh and cells
    cells = Cell[]
    bc_sides = BoundarySide[]
    
    # Process each element
    for elem_id in 1:nelems
        # Get element info
        elem_type = elem_info[1, elem_id]
        zone = elem_info[2, elem_id]
        offset_side = elem_info[3, elem_id]
        last_side = elem_info[4, elem_id]
        offset_node = elem_info[5, elem_id]
        last_node = elem_info[6, elem_id]
        
        # Determine element type based on the element type code

        if elem_type == 108  # Hexahedron
            # Load vertices so they are in CGNS ordering.
            vertices = Vector{Vertex}(undef, last_node - offset_node)
            vertices[1] = node_coords[:, offset_node + 1]
            vertices[2] = node_coords[:, offset_node + 2]
            vertices[3] = node_coords[:, offset_node + 4]
            vertices[4] = node_coords[:, offset_node + 3]
            vertices[5] = node_coords[:, offset_node + 5]
            vertices[6] = node_coords[:, offset_node + 6]
            vertices[7] = node_coords[:, offset_node + 8]
            vertices[8] = node_coords[:, offset_node + 7]

            # Calculate barycenter
            barycenter = sum(vertices) ./ 8
            if abs(barycenter[3]) > 1e-10
                error("Element $elem_id is not centered on z=0 (z=$(barycenter[3]))")
            end

            # Read boundary conditions
            bcs = Vector{UInt32}(undef, last_side - offset_side)
            neighbours = Vector{UInt32}(undef, last_side - offset_side)
            sides = get_sides(vertices)
            for i in 1:(last_side - offset_side)
                side_idx = offset_side + i
                if side_info[5, side_idx] == 0
                    bcs[i] = 0
                else
                    push!(
                        bc_sides,
                        BoundarySide(side_info[5, side_idx], side_normal(sides[i]), sides[i], face_quad_area(sides[i]), UInt32(elem_id), 0.0, 0)
                    )
                    bcs[i] = length(bc_sides)
                end
                neighbours[i] = side_info[3, side_idx]
            end
            normals = map(side_normal, sides)
            # Make sure that S1 and S6 are the Z- and Z+ faces.
            # TODO: Disable this check if in 3D case.
            if abs(normals[2][3]) > 1e-12 || abs(normals[3][3]) > 1e-12 || abs(normals[4][3]) > 1e-12 || abs(normals[5][3]) > 1e-12
                error("Element $elem_id side ordering is incorrect for 2D mesh.")
            end

            face_origins = map(v -> v[1], sides)
            
            cell = Hexahedron(vertices, normals, face_origins, barycenter, bcs, neighbours, 0.0)
            cell = @set cell.volume = cell_volume(cell)
            push!(cells, Cell(cell))
        else
            error("Invalid element type '$elem_type'")
        end
    end
    
    bc_names = map(lowercase, map(strip, read(h5file, "BCNames")[:, 1]))

    close(h5file)

    for cell in cells
        assert_normal_direction(cell)
    end

    # Return mesh
    return Mesh(cells, bc_names, bc_sides)
end


"""
Return a uniformly-distributed random point on the quadrilateral face `face`.

Arguments:
- `face`: SVector{4, SVector{3,Float64}} four corner vertices of the face

Returns an SVector{3,Float64} position.
"""
function random_point_on_face(face::SVector{4, SVector{3,Float64}}) :: SVector{3,Float64}
    # Split quadrilateral into two triangles: (1,2,3) and (1,3,4)
    # We use vertex 1 as the common vertex to handle non-planar quads consistently
    a, b, c, d = face[1], face[2], face[3], face[4]
    
    # Compute twice the area of each triangle (magnitude of cross product)
    # No need to divide by 2 or multiply by 0.5 since we only need the ratio
    area123 = norm(cross(b - a, c - a))
    area134 = norm(cross(c - a, d - a))
    
    # Select triangle with probability proportional to its area
    if rand() * (area123 + area134) < area123
        # Sample uniformly from triangle (a, b, c) using barycentric coordinates
        # Generate random barycentric coordinates: u ≥ 0, v ≥ 0, u+v ≤ 1
        r1, r2 = rand(), rand()
        s = sqrt(r1)
        u = 1 - s      # barycentric coordinate for vertex b
        v = s * r2     # barycentric coordinate for vertex c
        # w = 1 - u - v (coordinate for vertex a) is implicit
        return a + u*(b - a) + v*(c - a)
    else
        # Sample uniformly from triangle (a, c, d)
        r1, r2 = rand(), rand()
        s = sqrt(r1)
        u = 1 - s      # barycentric coordinate for vertex c
        v = s * r2     # barycentric coordinate for vertex d
        return a + u*(c - a) + v*(d - a)
    end
end


export Cell, Hexahedron
export cell_contains
export random_point_on_face