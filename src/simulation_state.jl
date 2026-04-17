# Leon's Gas Kinetics
# Copyright (C) 2026  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

mutable struct SimulationState{Features}
    # Particle information
    particles :: ParticleData{Features}
    # Save the number of particles per cell
    cell_part_count :: Vector{UInt32}
    # Current simulation time
    time :: Float64
    perf_counters :: PerformanceCounters
end