using LeonsGasKinetics

"""
Create a minimal simulation config.

Fields that are not explicitly set should not be assumed to contain valid data.
"""
function TestSimulationConfig()
    argon = LeonsGasKinetics.SpeciesConfig(4.65e-26, 3.66e-10, 273.0, 0.81, "Argon")
    species = [argon]
    vrbgk_config = LeonsGasKinetics.VRBGKConfig(false, 0.0, 0.0, false, 0.9)
    return LeonsGasKinetics.SimulationConfig(
        species,
        Vector{LeonsGasKinetics.Boundary}(),
        (pdata, samples, config, flow_variables, dt) -> (),
        100.0, # mpf
        1.0,  # t_end
        1e-6,  # dt
        1.0,  # sample_fraction
        1,  # output_interval
        1.0,  # report_interval
        "test",
        "test.msh",
        "output",
        true,  # silent
        false,  # asserts
        vrbgk_config,
        :continuous,
        3  # degrees_of_freedom
    )
end


function TestHexahedron(volume = 1)
    len = cbrt(volume)
    normals = SVector{6, SVector{3, Float64}}(
        [1, 0, 0],  # +x
        [-1, 0, 0], # -x
        [0, 1, 0],  # +y
        [0, -1, 0], # -y
        [0, 0, 1],  # +z
        [0, 0, -1]  # -z
    )

    origins = SVector{6, SVector{3, Float64}}(
        [len, 0.0, 0.0], # +x face
        [0.0, 0.0, 0.0], # -x face
        [0.0, len, 0.0], # +y face
        [0.0, 0.0, 0.0], # -y face
        [0.0, 0.0, len], # +z face
        [0.0, 0.0, 0.0]  # -z face
    )

    verts = SVector{8, SVector{3, Float64}}(
        [0.0, 0.0, 0.0],
        [len, 0.0, 0.0],
        [len, len, 0.0],
        [0.0, len, 0.0],
        [0.0, 0.0, len],
        [len, 0.0, len],
        [len, len, len],
        [0.0, len, len],
    )
    
    return Hexahedron(
        verts, normals, origins, 
        SVector(0.5, 0.5, 0.5) * len, # barycenter
        zeros(SVector{6, UInt32}), zeros(SVector{6, UInt32}), len^3
    )
end