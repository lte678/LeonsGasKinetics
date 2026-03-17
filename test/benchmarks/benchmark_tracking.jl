using LeonsGasKinetics

using BenchmarkTools
using StaticArrays
using LinearAlgebra
using Random
using InteractiveUtils
using Profile
using PProf


const N_TESTS = 10_000
Random.seed!(42)

# Cell setup
const normals = SVector(
    SVector(1.0, 0.0, 0.0), SVector(-1.0, 0.0, 0.0),
    SVector(0.0, 1.0, 0.0), SVector(0.0, -1.0, 0.0),
    SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -1.0)
)
const vertices = SVector(
    SVector(0.0, 0.5, 0.5), SVector(1.0, 0.5, 0.5),
    SVector(0.5, 0.0, 0.5), SVector(0.5, 1.0, 0.5),
    SVector(0.5, 0.5, 0.0), SVector(0.5, 0.5, 1.0)
)

# Generate Random Origins (within unit cube) and Directions (normalized)
const origins = [SVector{3, Float64}(rand(3)...) for _ in 1:N_TESTS]
const dirs = [normalize(SVector{3, Float64}(randn(3)...)) for _ in 1:N_TESTS]

# --- Benchmarking ---
println("Benchmarking find_exit_face with $N_TESTS iterations...")



# Benchmark wrapper function to ensure compiler doesn't optimize away the loop
function run_bench(origins, dirs)
    for i in 1:N_TESTS
        find_exit_face(origins[i], dirs[i], normals, vertices)
    end
end

run_bench(origins, dirs)

# Run the benchmark
#Profile.Allocs.clear()
#Profile.Allocs.@profile run_bench(origins, dirs)
#PProf.Allocs.pprof()

suite = @benchmarkable run_bench($origins, $dirs)
results = run(suite, verbose=true)

display(results)