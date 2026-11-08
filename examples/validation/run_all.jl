#!/usr/bin/env julia
# Run every comparison in this directory in turn and report which ones completed (see
# README.md). Each driver prints its own agreement figures; this only checks that it ran.
#
# Usage: julia examples/validation/run_all.jl [name ...]
#   e.g. julia examples/validation/run_all.jl gmst psd

const COMPARISONS = ["psd", "filtering", "resample", "spectrogram", "qtransform", "coherence",
                     "gmst", "antenna_response", "time_delay", "reference_curves"]

selected = isempty(ARGS) ? COMPARISONS : ARGS
failed = String[]
for name in selected
    name in COMPARISONS || error("unknown comparison '$name'; choose from $(join(COMPARISONS, ", "))")
    println("\n", "="^20, " compare_$name ", "="^20)
    try
        run(`$(Base.julia_cmd()) $(joinpath(@__DIR__, "compare_$name.jl"))`)
    catch err
        err isa ProcessFailedException || rethrow()
        push!(failed, name)
    end
end

println()
if isempty(failed)
    println("All $(length(selected)) comparisons ran.")
else
    println("Failed: ", join(failed, ", "))
    exit(1)
end
