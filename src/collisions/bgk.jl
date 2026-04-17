# Leon's Gas Kinetics
# Copyright (C) 2026  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
using Random
using InteractiveUtils

function swap!(array, idx1, idx2)
    tmp = array[idx1]
    array[idx1] = array[idx2]
    array[idx2] = tmp
end


function bgk_collision!(pdata, samples, config::SimulationConfig, flow_vars::FlowProperties, dt)
    n_part = length(pdata)
    if n_part < 2
        add_sample!(samples[:relaxation_rate], 0.0)
        return
    end

    if length(config.species) > 1
        error("\"BGK\" collision operator only supports single species.")
    end
    spec = config.species[1]

    if config.vrbgk.enabled
        vr_mass = vrbgk_calculate_mass(pdata)
        vr_mean = vr_mass / n_part
    end

    # Calculate the relax probability
    dyn_visc = dynamic_viscosity(spec)
    relax_freq = flow_vars.density*BOLTZMANN*spec.T_ref^(spec.omega + 0.5)*flow_vars.mean_temperature^(-spec.omega + 0.5) / dyn_visc
    relax_probability = 1.0 - exp(-dt*relax_freq)
    
    # Select particles to be relaxed
    n_relaxed = 0
    for i = 1:n_part
        if rand() < relax_probability
            # Relax the particle
            # swap!(pdata, n_relaxed + 1, i)
            particle = pdata[i]
            particle = @set particle.vel = sample_maxwellian(flow_vars.mean_temperature, flow_vars.velocity, spec.mass)

            if config.vrbgk.enabled
                f_eq = maxwellian(config.vrbgk.ref_temperature, [0.0, 0.0, 0.0]   , spec.mass, particle.vel)
                f    = maxwellian(flow_vars.mean_temperature  , flow_vars.velocity, spec.mass, particle.vel)
                
                particle = @set particle.features.vr_weight = vr_mean * f_eq / f
            end
            # Update particle
            pdata[i] = particle
            n_relaxed += 1
        end
    end

    if config.vrbgk.enabled
        vrbgk_conserve_mass(pdata, vr_mass; silent=config.silent)
    end

    # Update cell averages
    add_sample!(samples[:dynamic_viscosity], dyn_visc)
    add_sample!(samples[:relaxation_rate], n_relaxed / n_part)
end
