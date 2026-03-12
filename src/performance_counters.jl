struct PerformanceCounters
    # Cumulative times for all phases (nanoseconds)
    cumulative_times::Dict{Symbol, UInt64}
    
    # Start times for current timing block
    current_starts::Dict{Symbol, UInt64}
    
    function PerformanceCounters()
        return new(
            Dict(), Dict(),
        )
    end
end

"""
Starts alotting time to the execution step `phase`. May be nested. 
"""
function start_timing!(phase::Symbol, counters::PerformanceCounters)
    # Auto-register phase if it doesn't exist
    if !haskey(counters.cumulative_times, phase)
        counters.cumulative_times[phase] = 0
    end
    
    counters.current_starts[phase] = time_ns()
end

"""
Stops alotting time to the execution step `phase`.
"""
function end_timing!(phase::Symbol, counters::PerformanceCounters)
    # Ensure phase exists (should be true if start_timing! was called)
    if !haskey(counters.current_starts, phase)
        error("Phase $phase not registered. Call start_timing! first.")
    end
    
    start_time = counters.current_starts[phase]
    elapsed = time_ns() - start_time
    counters.cumulative_times[phase] += elapsed
end

"""
Time function using do-block syntax.
Example:
```
timed_region(counters, :advection) do
    ...
end
```

May be nested.
"""
function timed_region(f::Function, counters::PerformanceCounters, phase::Symbol)
    start_timing!(phase, counters)
    try
        return f()
    finally
        end_timing!(phase, counters)
    end
end


function print_performance_stats(counters::PerformanceCounters)
    println("=== Performance Statistics ===")
    
    # Sort phases by total execution time.
    phases = collect(keys(counters.cumulative_times))
    times = collect(values(counters.cumulative_times))
    sort_idx = sortperm(times, rev=true)    
    phases = phases[sort_idx]
    times = times[sort_idx]
    
    # Print each phase
    for (phase, time) in zip(phases, times)
        time_s = round(time / 1e9, digits=3)
        println("  $(rpad(string(phase), 20)) $(lpad(time_s, 8)) s")
    end
end
