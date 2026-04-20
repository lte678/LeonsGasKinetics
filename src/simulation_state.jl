# Leon's Gas Kinetics
# Copyright (C) 2026  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

  mutable struct SimulationState{ParticleFeatures, CellFeatures}
      # Particle information
      particles :: ParticleData{ParticleFeatures}
      # Per-cell data (particle counts and optional feature data)
      cells :: CellData{CellFeatures}
      # Current simulation time
      time :: Float64
      perf_counters :: PerformanceCounters
  end
