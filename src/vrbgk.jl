function vrbgk_check_weight(weight::Float64)
    if isnan(weight) || weight < 1e-3 || weight > 1e3
        error("Weight of particle is extreme, w = $weight")
    end
end

function vrbgk_conserve_mass(particles::ParticleDataView, target_vr_sum::Float64; silent::Bool=false)
    new_vr_sum = 0.0
    for i in 1:length(particles)
        new_vr_sum += particles[i].features.vr_weight
    end

    alpha = target_vr_sum / new_vr_sum

    for i in 1:length(particles)
        p = particles[i]
        particles[i] = @set p.features.vr_weight *= alpha
    end

    if !silent && abs(new_vr_sum - target_vr_sum) / target_vr_sum > 0.01
        println("Large drift in weights!")
    end
end