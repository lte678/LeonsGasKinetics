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

# Want to see something crazy? Float32 is slower than Float64 (on my Ryzen 3600X machine that is)
# Remember kids, always benchmark first.
T = Float64

# Cell setup
const normals = SVector{6, SVector{3, T}}(
    SVector(1.0, 0.0, 0.0), SVector(-1.0, 0.0, 0.0),
    SVector(0.0, 1.0, 0.0), SVector(0.0, -1.0, 0.0),
    SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -1.0)
)
const vertices = SVector{6, SVector{3, T}}(
    SVector(0.0, 0.5, 0.5), SVector(1.0, 0.5, 0.5),
    SVector(0.5, 0.0, 0.5), SVector(0.5, 1.0, 0.5),
    SVector(0.5, 0.5, 0.0), SVector(0.5, 0.5, 1.0)
)

# Generate Random Origins (within unit cube) and Directions (normalized)
const origins = [SVector{3, T}(rand(3)...) for _ in 1:N_TESTS]
const dirs = [SVector{3, T}(randn(3)...) for _ in 1:N_TESTS]
const dirs2 = [SVector{3, T}(randn(3)...) for _ in 1:N_TESTS]

# --- Benchmarking ---
println("Benchmarking find_exit_face with $N_TESTS iterations...")



# Benchmark wrapper function to ensure compiler doesn't optimize away the loop
function run_bench(origins, dirs)
    for i in 1:N_TESTS
        find_exit_face(origins[i], dirs[i], normals, vertices)
    end
end

function run_bench2d(origins, dirs)
    total_dist = 0.0
    for i in 1:N_TESTS
        coll_i, coll_dist = find_exit_face_2d(origins[i], dirs[i], normals, vertices)
        total_dist += coll_dist
    end
    return total_dist
end

run_bench(origins, dirs)

# Run the benchmark
#Profile.Allocs.clear()
#Profile.Allocs.@profile run_bench(origins, dirs)
#PProf.Allocs.pprof()

suite = @benchmarkable run_bench2d($origins, $dirs)
GC.gc()
results = run(suite, verbose=true)

display(results)