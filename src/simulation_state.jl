# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

mutable struct SimulationState
    # Particle information
    particles :: ParticleData
    # Save the number of particles per cell
    cell_part_count :: Vector{UInt32}
    # Current simulation time
    time :: Float64
end