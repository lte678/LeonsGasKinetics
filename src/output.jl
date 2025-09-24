# Leon's Gas Kinetics
# Copyright (C) 2025  Leon Teichroeb
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

using HDF5

"""
Routines to output the simulation result to an HDF5 file.
"""


function write_piclas_attribute(file, name, value)
    if typeof(value) == String
        # For some reason piclas uses this ridiculous string format.
        str_dtype = HDF5.API.h5t_copy(HDF5.API.H5T_C_S1)
        HDF5.API.h5t_set_size(str_dtype, 255)
        HDF5.API.h5t_set_strpad(str_dtype, HDF5.API.H5T_STR_SPACEPAD)
        str_dtype = HDF5.Datatype(str_dtype)
        str_dspace = HDF5.Dataspace(HDF5.API.h5s_create_simple(1, [1], [1]))
        
        attr = create_attribute(file, name, str_dtype, str_dspace)
        write_attribute(attr, str_dtype, rpad(value, 255))
    elseif typeof(value) <: AbstractArray{String}
        str_dtype = HDF5.API.h5t_copy(HDF5.API.H5T_C_S1)
        HDF5.API.h5t_set_size(str_dtype, 255)
        HDF5.API.h5t_set_strpad(str_dtype, HDF5.API.H5T_STR_SPACEPAD)
        str_dtype = HDF5.Datatype(str_dtype)
        str_dspace = HDF5.Dataspace(HDF5.API.h5s_create_simple(1, [length(value)], [length(value)]))
        attr = create_attribute(file, name, str_dtype, str_dspace)
        data = Vector{UInt8}(undef, 255 * length(value))
        for (i, s) in enumerate(value)
            bytes = codeunits(s)
            n = min(length(bytes), 255)
            start_idx = (i-1)*255 + 1
            data[start_idx:start_idx+n-1] .= bytes[1:n]
            # Fill remainder with spaces if needed
            if n < 255
                data[start_idx+n:start_idx+254] .= UInt8(' ')
            end
        end
        HDF5.API.h5a_write(attr, str_dtype, pointer(data))
    else
        write_attribute(file, name, [value])
    end
end


function write_macro_vals_hdf5(path::String, config::SimulationConfig, macro_vals::Vector{FlowProperties}, time::Float64, compat_mode)
    h5file = h5open(path, "w")

    if compat_mode
        write_attr_func = write_piclas_attribute
    else
        write_attr_func = write_attribute
    end

    if compat_mode
        write_attr_func(h5file, "Program", "PICLas")
        write_attr_func(h5file, "Piclas_Version", "3.6.0")
        write_attr_func(h5file, "Piclas_VersionInt", 30600)
        write_attr_func(h5file, "File_Version", 3.6)
    else
        write_attr_func(h5file, "Program", "LeonsGasKinetics")
    end
    write_attr_func(h5file, "File_Type", "DSMCState")
    write_attr_func(h5file, "Project_Name", config.project_name)
    write_attr_func(h5file, "MeshFile", relpath(config.mesh_file, dirname(path)))
    write_attr_func(h5file, "NSpecies", length(config.species))
    write_attr_func(h5file, "Time", time)

    var_names = [
        "Total_VeloX",
        "Total_VeloY",
        "Total_VeloZ",
        "Total_TempTransX",
        "Total_TempTransY",
        "Total_TempTransZ",
        "Total_NumberDensity",
        "Total_SimPartNum",
        "Total_TempTransMean",
    ]

    # Write dataset
    write_attr_func(h5file, "VarNamesAdd", var_names)
    h5_array = zeros(9, length(macro_vals))
    for i = 1:length(macro_vals)
        h5_array[1:3, i] = macro_vals[i].velocity
        h5_array[4:6, i] = macro_vals[i].temperature
        h5_array[7, i] = macro_vals[i].density
        h5_array[8, i] = macro_vals[i].sim_particle_count
        h5_array[9, i] = macro_vals[i].mean_temperature
    end
    write_dataset(h5file, "ElemData", h5_array)

    close(h5file)
end