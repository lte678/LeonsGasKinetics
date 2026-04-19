# Leon's Gas Kinetics
# Copyright (C) 2026  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

using StaticArrays
using LightSumTypes

handle_boundary(particle::SingleParticle, normal::SVector{3}, species::SpeciesConfig, boundary::Boundary, config::SimulationConfig) =
    handle_boundary(particle, normal, species, variant(boundary), config)

include("symmetric.jl")
include("reflective.jl")
include("diffuse.jl")
