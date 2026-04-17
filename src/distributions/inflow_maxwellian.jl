# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Abramowitz & Stegun 7.1.26 polynomial approximation, max error ≈ 1.5e-7.
function _erf(x::Float64) :: Float64
    t = 1.0 / (1.0 + 0.3275911 * abs(x))
    poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))))
    return copysign(1.0 - poly * exp(-x^2), x)
end


"""
Mean speed of inward-crossing particles through a boundary face from a Maxwellian reservoir.

This is the half-range Maxwell flux speed (Chapman–Enskog):

    v_thermal = sqrt(2 k_B T / m)
    a         = (u · n_in) / v_thermal          # projected speed ratio
    v_flux    = v_thermal / (2√π) * (exp(-a²) + a√π (1 + erf(a)))

Arguments are the same as `sample_inflow_maxwellian`. Returns v_flux [m/s].
"""
function mean_inflow_flux_speed(
    inward_normal::SVector{3,Float64},
    temperature::Float64,
    bulk_velocity::SVector{3,Float64},
    mass::Float64,
) :: Float64
    v_thermal = sqrt(2.0 * BOLTZMANN * temperature / mass)
    a = dot(bulk_velocity, inward_normal) / v_thermal
    return v_thermal / (2.0 * sqrt(PI)) * (exp(-a^2) + a * sqrt(PI) * (1.0 + _erf(a)))
end


"""
Sample the velocity of a particle entering the domain through a boundary face from a
Maxwellian reservoir, conditioned on the particle actually crossing the face inward.

The inward-crossing velocity distribution is not a plain Maxwellian: it is the full
Maxwellian f(v) weighted by the inward normal flux |v · n_in|, where n_in is the
inward-pointing unit normal of the face.  Sampling this distribution correctly ensures
that the flux of particles through the face matches the Chapman–Enskog prediction for
the given reservoir state (temperature, density, bulk velocity).

Arguments:
- `inward_normal`: unit normal pointing into the domain at the face (SVector{3,Float64})
- `temperature`  : reservoir temperature [K]
- `bulk_velocity`: reservoir bulk velocity [m/s]  (SVector{3,Float64})
- `mass`         : particle mass [kg]

Returns an SVector{3,Float64} velocity drawn from the flux-weighted distribution.
"""
function sample_inflow_maxwellian(
    inward_normal::SVector{3,Float64},
    temperature::Float64,
    bulk_velocity::SVector{3,Float64},
    mass::Float64,
) :: SVector{3,Float64}
    error("sample_inflow_maxwellian: not yet implemented")
end
