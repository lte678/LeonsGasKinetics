"""
Converts raw samples into macroscopic cell values.
`accumulators` is a an array of length N_cell containing an ordered dict of `Averager`
"""
function postprocess(accumulators, config, mesh)
    means = [OrderedDict(k => mean(v) for (k, v) in acc) for acc in accumulators]

    output = OrderedDict{String, Vector{Float64}}()

    flow_vars = map(
        (mean, c) -> calc_flow_properties(
            mean[:count],
            [mean[:c_x], mean[:c_y], mean[:c_z]],
            [mean[:c_xx], mean[:c_yy], mean[:c_zz]],
            config,
            c.volume
        ),
        means, mesh.cells
    )

    output = merge(output, OrderedDict{String, Vector{Float64}}(
        "Total_VeloX" => [fp.velocity[1] for fp in flow_vars],
        "Total_VeloY" => [fp.velocity[2] for fp in flow_vars],
        "Total_VeloZ" => [fp.velocity[3] for fp in flow_vars],
        "Total_TempTransX" => [fp.temperature[1] for fp in flow_vars],
        "Total_TempTransY" => [fp.temperature[2] for fp in flow_vars],
        "Total_TempTransZ" => [fp.temperature[3] for fp in flow_vars],
        "Total_NumberDensity" => [fp.density for fp in flow_vars],
        "Total_SimPartNum" => [fp.sim_particle_count for fp in flow_vars],
        "Total_TempTransMean" => [fp.mean_temperature for fp in flow_vars]
    ))

    if haskey(means[1], :dynamic_viscosity) 
        output["BGK_Viscosity"] = [acc[:dynamic_viscosity] for acc in means]
    end
    
    if haskey(means[1], :relaxation_rate)
        output["BGK_MeanRelaxationRate"] = [acc[:relaxation_rate] for acc in means]
    end

    return output
end