# Leon's Gas Kinetics
# Copyright (C) 2026  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

using StaticArrays
using LightSumTypes

handle_boundary(particle::SingleParticle, normal::SVector{3}, species::SpeciesConfig, boundary::Boundary, config::SimulationConfig) =
    handle_boundary(particle, normal, species, variant(boundary), config)


function perpendicular_pair(n::SVector{3,Float64})
    # Duff et al. 2017 - robust across all normals
    sign_n = copysign(1.0, n[3])
    a = -1.0 / (sign_n + n[3])
    b = n[1] * n[2] * a

    tang1 = SVector(1.0 + sign_n * n[1]^2 * a,  sign_n * b,         -sign_n * n[1])
    tang2 = SVector(b,                             n[2]^2 * a + sign_n, -n[2])
    return tang1, tang2
end


include("symmetric.jl")
include("reflective.jl")
include("diffuse.jl")
