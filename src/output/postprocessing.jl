"""
Converts raw samples into macroscopic cell values.
`accumulators` is a an array of length N_cell containing an ordered dict of `Averager`
"""
function postprocess(accumulators, config, mesh)
    means = [OrderedDict(k => mean(v) for (k, v) in acc) for acc in accumulators]

    output = OrderedDict{String, Vector{Float64}}(
        "Total_VeloX" => [acc[:u_x] for acc in means],
        "Total_VeloY" => [acc[:u_y] for acc in means],
        "Total_VeloZ" => [acc[:u_z] for acc in means],
        "Total_TempTransX" => [acc[:T_x] for acc in means],
        "Total_TempTransY" => [acc[:T_y] for acc in means],
        "Total_TempTransZ" => [acc[:T_z] for acc in means],
        "Total_NumberDensity" => [acc[:density] for acc in means],
        "Total_SimPartNum" => [acc[:sim_part_count] for acc in means],
        "Total_TempTransMean" => [(acc[:T_x] + acc[:T_y] + acc[:T_z])/3.0 for acc in means],
    )

    if haskey(means[1], :dynamic_viscosity) 
        output["BGK_Viscosity"] = [acc[:dynamic_viscosity] for acc in means]
    end
    
    if haskey(means[1], :relaxation_rate)
        output["BGK_MeanRelaxationRate"] = [acc[:relaxation_rate] for acc in means]
    end

    if haskey(means[1], :vr_weight)
        output["VRBGK_MeanImportance"] = [acc[:vr_weight] for acc in means]
    end

    return output
end
