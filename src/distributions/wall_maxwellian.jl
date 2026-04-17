function sample_wall_distribution(wall_temp, mass_ic)
    # Thermal velocity
    cmr = sqrt(BOLTZMANN * wall_temp / mass_ic)
    return SVector(
        cmr * randn(),
        cmr * randn(),
        cmr * sqrt(2*randexp())
    )
end