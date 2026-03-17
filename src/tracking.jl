# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

using StaticArrays

function find_exit_face(origin::SVector{3,T}, dir::SVector{3,T}, normals::SVector{6,SVector{3,T}}, vertices::SVector{6,SVector{3,T}}) where T
    hit_distance = @MVector zeros(6)

    for i in 1:6
        # n should point INWARD from the cell
        n = normals[i]
        p0 = vertices[i]
        
        denom = dot(n, dir)
        
        # Only hit faces the particle is moving towards (outward normal)
        # The result of this may be negative if we are slightly past the face. We should still collide with
        # such faces, since it means we missed the collision at some other point, maybe due to rounding errors.
        # On the way back, the direction of the normal will prevent it from hitting the face as it passes back
        # back through.
        hit_distance[i] = ifelse(denom < 0.0, dot(p0 - origin, n) / denom, Inf)
    end
    
    hit_index = argmin(hit_distance)
    return hit_index, hit_distance[hit_index]
end