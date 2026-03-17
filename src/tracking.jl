# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

using StaticArrays

function find_exit_face(origin::SVector{3,T}, dir::SVector{3,T}, normals::SVector{6,SVector{3,T}}, vertices::SVector{6,SVector{3,T}}) where T
    best_dist = Inf
    best_idx  = 0

    @inbounds for i in 1:6
        n = normals[i]
        p0 = vertices[i]
        denom = dot(n, dir)
        d = ifelse(denom < 0.0, dot(p0 - origin, n) / denom, Inf)
        if d < best_dist
            best_dist = d
            best_idx  = i
        end
    end

    return best_idx, best_dist

end