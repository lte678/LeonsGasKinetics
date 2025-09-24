# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
using Random

function swap!(array, idx1, idx2)
    tmp = array[idx1]
    array[idx1] = array[idx2]
    array[idx2] = tmp
end


function bgk_collision!(particle_x, particle_v, config::SimulationConfig, flow_vars::FlowProperties, dt)
    if length(particle_x) < 2
        return
    end
    if length(config.species) > 1
        error("\"bgk\" collision operator only supports single species.")
    end
    spec = config.species[1]
    n_part = length(particle_x)
    @assert length(particle_v) == n_part

    # Calculate the relax probability
    dyn_visc = dynamic_viscosity(spec)
    relax_freq = flow_vars.density*BOLTZMANN*spec.T_ref^(spec.omega + 0.5)*flow_vars.mean_temperature^(-spec.omega + 0.5) / dyn_visc
    relax_probability = 1.0 - exp(-dt*relax_freq)
    
    # Select particles to be relaxed and move them to the beginning of the particle array.
    # After this, 1:n_relaxed should be relaxed and n_relaxed+1: should be left alone.    
    n_relaxed = 0
    for i = 1:n_part
        if rand() < relax_probability
            swap!(particle_x, n_relaxed + 1, i)
            swap!(particle_v, n_relaxed + 1, i)
            n_relaxed += 1
        end
    end

    # Do relaxation
    if n_relaxed > 0
        sample_v = @view particle_v[1:n_relaxed - 1]
        sample_maxwellian!(sample_v, flow_vars.mean_temperature, flow_vars.velocity, spec.mass)
    end
end
