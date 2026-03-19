using LeonsGasKinetics
using Profile
using PProf


BASE_FOLDER = dirname(dirname(pathof(LeonsGasKinetics)))
config_path = joinpath(BASE_FOLDER, "examples", "couette", "couette.toml")

Profile.clear()
@profile run_simulation_from_config(config_path, "output")
pprof()