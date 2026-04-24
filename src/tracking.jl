# Leon's Gas Kinetics
# Copyright (C) 2026  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

using StaticArrays


function find_exit_face_3d(origin::AbstractArray, dir::AbstractArray, normals::SVector{6,SVector{3,T}}, vertices::SVector{6,SVector{3,T}}) where T
    best_dist = Inf
    best_idx  = 0

    @inbounds for i in 1:6
        n = normals[i]
        p0 = vertices[i]
        # dot(n, dir)
        denom = n[1]*dir[1] + n[2]*dir[2] + n[3]*dir[3]
        # dot(p0 - origin, n)
        normal_dist = n[1]*(p0[1] - origin[1]) + n[2]*(p0[2] - origin[2]) +  n[3]*(p0[3] - origin[3])
        d = ifelse(denom < 0.0, normal_dist / denom, Inf)
        if d < best_dist
            best_dist = d
            best_idx  = i
        end
    end

    return best_idx, best_dist
end


"""
Specialized routine that calculates the exit face in the 2D case.
The z-axis is completely neglected in this case. None of the first 4 face normals may have a non-zero z-component.
"""
function find_exit_face_2d(origin::AbstractArray, dir::AbstractArray, normals::SVector{6,SVector{3,T}}, vertices::SVector{6,SVector{3,T}}) where T
    best_dist = Inf
    best_idx  = 0

    @inbounds for i in 2:5
        n = normals[i]
        p0 = vertices[i]
        # dot(n, dir)
        denom = n[1]*dir[1] + n[2]*dir[2]
        # dot(p0 - origin, n)
        normal_dist = n[1]*(p0[1] - origin[1]) + n[2]*(p0[2] - origin[2])
        d = ifelse(denom < 0.0, normal_dist / denom, Inf)
        if d < best_dist
            best_dist = d
            best_idx  = i
        end
    end

    return best_idx, best_dist
end


find_exit_face(origin, dir, normals, vertices) = find_exit_face_3d(origin, dir, normals, vertices)

function find_exit_face(origin, dir, normals, vertices, degrees_of_freedom::Val{3})
    return find_exit_face_3d(origin, dir, normals, vertices)
end

function find_exit_face(origin, dir, normals, vertices, degrees_of_freedom::Val{2})
    return find_exit_face_2d(origin, dir, normals, vertices)
end

export find_exit_face, find_exit_face_2d
