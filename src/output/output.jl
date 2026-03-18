function write_flow_state(volume_samples, t, mesh, sim_config)
    volume_output = postprocess(volume_samples, sim_config, mesh)
    
    output_file = @sprintf "%s_DSMCState_%012.8f.h5" sim_config.project_name t
    @printf "Writing flow state to %s...\n" output_file
        
    write_volume_data_hdf5(
        joinpath(sim_config.output_path, output_file),
        sim_config,
        volume_output,
        t,
        true
    )
end