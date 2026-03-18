## How to run
First, install [gmsh](https://gmsh.info) and [PyHOPE](https://hopr-framework.github.io/PyHOPE).

```
gmsh knudsen_pump.geo
pyhope pyhope.ini
cd 'project root'
julia --project=. src/main.jl examples/knudsen_pump/sim.toml
```