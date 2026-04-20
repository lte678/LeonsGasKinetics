# TODO: This sucks. Lets see if we can get rid of it in the future.
const _weight_drift_warned = Ref(false)

function vrbgk_check_weight(weight::Float64)
    if isnan(weight) || weight < 1e-3 || weight > 1e3
        error("Weight of particle is extreme, w = $weight")
    end
end

function vrbgk_calculate_mass(particles)
    vr_sum = 0.0
    for i in 1:length(particles)
        vr_sum += particles[i].features.vr_weight
    end
    return vr_sum
end

"""
This function pushes the variance reduction weights from local to global reference.
In other words, w=f_eq,loc/f --> w=f_eq,glob/f.
"""
function vrbgk_local_to_global!(particles, sim_config, cell_data)
    mass = sim_config.species[1].mass
    T_glob = sim_config.vrbgk.ref_temperature
    n_glob = sim_config.vrbgk.ref_density
    
    @inbounds for i = 1:length(particles)   
        cell_i = particles.cell[i]
        v = particles.vel[i]

        T_loc = cell_data.features.vrbgk_ref_temperature[cell_i]
        n_loc = cell_data.features.vrbgk_ref_density[cell_i]
        particles.features.vr_weight[i] *=
          (n_glob * maxwellian(T_glob, SVector(0.0, 0.0, 0.0), mass, v)) /
          (n_loc  * maxwellian(T_loc , SVector(0.0, 0.0, 0.0), mass, v))
    end
end


"""
This function pushes the variance reduction weights from global to local reference.
In other words, w=f_eq,glob/f --> w=f_eq,loc/f.
"""
function vrbgk_global_to_local!(particles, sim_config, cell_data)
    mass = sim_config.species[1].mass
    T_glob = sim_config.vrbgk.ref_temperature
    n_glob = sim_config.vrbgk.ref_density
    
    @inbounds for i = 1:length(particles)   
        cell_i = particles.cell[i]
        v = particles.vel[i]

        T_loc = cell_data.features.vrbgk_ref_temperature[cell_i]
        n_loc = cell_data.features.vrbgk_ref_density[cell_i]
        particles.features.vr_weight[i] *=
          (n_loc  * maxwellian(T_loc , SVector(0.0, 0.0, 0.0), mass, v)) /
          (n_glob * maxwellian(T_glob, SVector(0.0, 0.0, 0.0), mass, v))
    end
end


function vrbgk_conserve_mass(particles, target_vr_sum::Float64; silent::Bool=false)
    new_vr_sum = vrbgk_calculate_mass(particles)

    alpha = target_vr_sum / new_vr_sum

    for i in 1:length(particles)
        p = particles[i]
        particles[i] = @set p.features.vr_weight *= alpha
    end

    if !silent && !_weight_drift_warned[] && abs(new_vr_sum - target_vr_sum) / target_vr_sum > 0.01
        println("Large drift in weights!")
        _weight_drift_warned[] = true
    end
end