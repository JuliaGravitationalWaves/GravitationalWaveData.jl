#!/usr/bin/env julia
# GravitationalWaveData.jl side of compare_gmst.jl (see README.md): computes `gmst` at each given GPS
# time and writes it to a CSV.
#
# Usage: julia --project=<GravitationalWaveData.jl root> compare_gmst_gwdata_step.jl <out_csv> <gps1> <gps2> ...

using GravitationalWaveData: gmst

out_csv = ARGS[1]
gps_times = parse.(Float64, ARGS[2:end])

open(out_csv, "w") do io
    println(io, "gps,gmst")
    for g in gps_times
        println(io, "$(repr(g)),$(repr(gmst(g)))")
    end
end
