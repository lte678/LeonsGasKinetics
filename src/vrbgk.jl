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