function vrbgk_check_weight(weight::Float64)
    if isnan(weight) || weight < 1e-3 || weight > 1e3
        error("Weight of particle is extreme, w = $weight")
    end
end