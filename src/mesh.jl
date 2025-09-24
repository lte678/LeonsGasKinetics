# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

using LinearAlgebra
using StaticArrays
using HDF5
using Printf
using LightSumTypes

Vertex = SVector{3, Float64}

struct Hexahedron
    vertices :: SVector{8, Vertex}
    barycenter :: SVector{3, Float64}
    # Boundary conditions
    bcs :: SVector{6, Int16}
    # Neighbour cell for side.
    neighbours :: SVector{6, UInt}
    volume :: Float64
end

@sumtype Cell(Hexahedron)

struct Mesh
    cells :: Vector{Cell}
    bc_names :: Vector{String}
end


function get_sides(cell :: Hexahedron)
    SVector(
        SVector(cell.vertices[1], cell.vertices[4], cell.vertices[3], cell.vertices[2]),
        SVector(cell.vertices[1], cell.vertices[2], cell.vertices[6], cell.vertices[5]),
        SVector(cell.vertices[2], cell.vertices[3], cell.vertices[7], cell.vertices[6]),
        SVector(cell.vertices[3], cell.vertices[4], cell.vertices[8], cell.vertices[7]),
        SVector(cell.vertices[1], cell.vertices[5], cell.vertices[8], cell.vertices[4]),
        SVector(cell.vertices[5], cell.vertices[6], cell.vertices[7], cell.vertices[8])
    )
end

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
    barycenters = read(h5file, "ElemBarycenters")
    node_coords = read(h5file, "NodeCoords")
    side_info = read(h5file, "SideInfo")
    
    # Initialize mesh and cells
    cells = Cell[]
    
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
            barycenter = barycenters[:, elem_id]

            # Read boundary conditions
            bcs = Vector{Int16}(undef, last_side - offset_side)
            neighbours = Vector{UInt}(undef, last_side - offset_side)
            for i in 1:(last_side - offset_side)
                side_idx = offset_side + i
                bcs[i] = side_info[5, side_idx]
                neighbours[i] = side_info[3, side_idx]
            end
            
            cell = Hexahedron(vertices, barycenter, bcs, neighbours, 0.0)
            cell = @set cell.volume = cell_volume(cell)
            push!(cells, Cell(cell))
        else
            error("Invalid element type '$elem_type'")
        end
    end
    
    bc_names = read(h5file, "BCNames")[:, 1]

    close(h5file)

    for cell in cells
        assert_normal_direction(cell)
    end

    # Return mesh
    return Mesh(cells, bc_names)
end
