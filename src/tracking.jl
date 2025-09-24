# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

using StaticArrays

# Tests whether a ray intersects with a planar face.
function ray_intersects_rect(origin::SVector{3, Float64}, dir::SVector{3, Float64}, rect::SVector{4, SVector{3, Float64}})
    # Split quad into two triangles
    tri1 = @SVector [rect[1], rect[2], rect[3]]
    tri2 = @SVector [rect[1], rect[3], rect[4]]

    # Test intersection with first triangle
    hit, intersection = ray_intersects_tri(origin, dir, tri1)
    if hit
        return true, intersection
    end

    # Test intersection with second triangle
    hit, intersection = ray_intersects_tri(origin, dir, tri2)
    return hit, intersection
end


# Tests whether a ray intersects with a triangle.
function ray_intersects_tri(origin::SVector{3}, dir::SVector{3}, tri::SVector{3, Vertex})
    e1 = tri[2] - tri[1]
    e2 = tri[3] - tri[1]

    ray_cross_e2 = cross(dir, e2)
    det = dot(e1, ray_cross_e2)

    # Check if ray is parallel to triangle
    #if abs(det) < eps(Float64)
    # Check that we are colliding with the front-face
    if det > -eps(Float64)
        return false, @SVector [0.0, 0.0, 0.0]
    end
    
    inv_det = 1.0 / det
    s = origin - tri[1]
    u = inv_det * dot(s, ray_cross_e2)
    
    if u < -COLLISION_TOL || u > 1.0 + COLLISION_TOL
        return false, @SVector [0.0, 0.0, 0.0]
    end

    s_cross_e1 = cross(s, e1)
    v = inv_det * dot(dir, s_cross_e1)

    if v < -COLLISION_TOL || u + v > 1.0 + COLLISION_TOL
        return false, @SVector [0.0, 0.0, 0.0]
    end

    t = inv_det * dot(e2, s_cross_e1)

    # This ensures that particles only hit what is in front of them.
    # We weaken this restriction slightly so that they only need to be moving in the direction of the boundary.
    #if t > eps(Float64)
    #    intersection_point = origin + dir * t
    #    return intersection_point
    #else
    #    return nothing
    #end
    return true, origin + dir * t
end