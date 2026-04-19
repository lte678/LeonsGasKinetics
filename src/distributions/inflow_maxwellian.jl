# Leon's Gas Kinetics
# Copyright (C) 2026  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

using Random
using LinearAlgebra
using StaticArrays


"""
Calculate the error function at x.

Source: Abramowitz & Stegun 7.1.26 polynomial approximation, max error ≈ 1.5e-7.  
Source: https://en.wikipedia.org/wiki/Error_function
"""
function erf(x::Float64) :: Float64
    t = 1.0 / (1.0 + 0.3275911 * abs(x))
    poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))))
    return copysign(1.0 - poly * exp(-x^2), x)
end


"""
Mean speed of inward-crossing particles through a boundary face from a Maxwellian reservoir.
This quantity can be used to calculate the rate of particles crossing a boundary.

This is the half-range Maxwell flux speed (Chapman–Enskog):

v_flux    = v_thermal / (2sqrt(π)) * (exp(-a^2) + a*sqrt(π) (1 + erf(a)))
with normalized velocity a = dot(u, n_in) / v_thermal
"""
function mean_inflow_flux_speed(
    normal::SVector{3,Float64},
    temperature::Float64,
    bulk_velocity::SVector{3,Float64},
    mass::Float64,
) :: Float64
    v_thermal = sqrt(2.0 * BOLTZMANN * temperature / mass)
    a = dot(bulk_velocity, normal) / v_thermal
    return v_thermal / (2.0 * sqrt(PI)) * (exp(-a^2) + a * sqrt(PI) * (1.0 + erf(a)))
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
- `normal`       : unit normal pointing into the domain at the face (Vector{3,Float64})
- `temperature`  : reservoir temperature [K]
- `bulk_velocity`: reservoir bulk velocity [m/s]  (Vector{3,Float64})
- `mass`         : particle mass [kg]

Returns an SVector{3,Float64} velocity drawn from the flux-weighted distribution.

Source: A.L. Garcia, W. Wagner, "_Generation of the Maxwellian inflow distribution_". 2006
"""
function sample_inflow_maxwellian(
    normal::SVector{3,Float64},
    temperature::Float64,
    bulk_velocity::SVector{3,Float64},
    mass::Float64,
)::SVector{3,Float64}
    # Thermal speed v_T = sqrt(2 * k_B * T / m)
    v_t = sqrt(2.0 * BOLTZMANN * temperature / mass)

    # Speed ratio a = (V · n) / v_T
    a = dot(normal, bulk_velocity) / v_t
    
    # Sample the normal velocity as z_star using rejection sampling
    z_star = zero(Float64)
    if a < 0.0 # Envelope 2 (for a < 0)
        z_a    = 0.5 * (a - sqrt(abs2(a) + 2.0))
        beta_a = a - (1.0 - a) * (a - z_a)
        
        # Probability of selecting the exponential branch vs the uniform branch
        denom = exp(-abs2(beta_a)) + 2.0 * (a - z_a) * (a - beta_a) * exp(-abs2(z_a))
        p1    = exp(-abs2(beta_a)) / denom
        
        while true
            if rand() < p1
                # Exponential branch
                z_star = -sqrt(abs2(beta_a) - log(rand()))
                if (a - z_star) / (-z_star) > rand()
                    break
                end
            else
                # Uniform branch
                z_star = beta_a + (a - beta_a) * rand()
                if ((a - z_star) / (a - z_a)) * exp(abs2(z_a) - abs2(z_star)) > rand()
                    break
                end
            end
        end
    else # Envelope 4 (for a >= 0)
        p1 = 1.0 / (2.0 * a * sqrt(PI) + 1.0)
        while true
            if rand() < p1
                z_star = -sqrt(-log(rand()))
            else
                z_star = randn() / sqrt(2.0)
            end
            
            # Acceptance probability
            acc = iszero(a) ? 1.0 : (a - z_star) / a
            if acc > rand()
                break
            end
        end
    end
    
    # Final velocity assembly
    t1, t2 = perpendicular_pair(normal)
    v = bulk_velocity - v_t * z_star * normal + (v_t / sqrt(2.0)) * (t1 * randn() + t2 * randn())
    return v 
end