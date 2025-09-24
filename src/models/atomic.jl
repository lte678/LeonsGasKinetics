# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

function dynamic_viscosity(s::SpeciesConfig)
    # From PICLas
    return 30.0*sqrt(s.mass*BOLTZMANN*s.T_ref/PI) /
        (4.0*(4.0 - 2.0*s.omega) * (6.0 - 2.0*s.omega)*s.d_ref^2)
end