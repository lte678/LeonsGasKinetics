# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

using StaticArrays

"""
Generate Maxwell-Boltzmann distributed velocities.
"""
function sample_maxwellian(temperature, bulk_velocity, mass, n_particles)
    velocities = Vector{SVector{3, Float64}}(undef, n_particles)
    sample_maxwellian!(velocities, temperature, bulk_velocity, mass)
    return velocities
end


"""
Generate Maxwell-Boltzmann distributed velocities.
"""
function sample_maxwellian!(velocities, temperature, bulk_velocity, mass)   
    # Thermal velocity (most probable speed in Maxwellian distribution)
    v_th = sqrt(BOLTZMANN * temperature / mass)
    for i in 1:length(velocities)
        # Total velocity = bulk velocity + thermal velocity
        velocities[i] = bulk_velocity .+ v_th * randn(3)
    end
end
