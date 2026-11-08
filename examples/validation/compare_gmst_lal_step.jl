#!/usr/bin/env julia
# Reference side of compare_gmst.jl (see README.md): computes
# lal.GreenwichMeanSiderealTime at each given GPS time and writes it to a CSV. Runs in this
# directory's project and must not load GWData.
#
# Usage: julia --project=examples/validation compare_gmst_lal_step.jl <out_csv> <gps1> <gps2> ...

using PythonCall
const lal = pyimport("lal")

out_csv = ARGS[1]
gps_times = parse.(Float64, ARGS[2:end])

open(out_csv, "w") do io
    println(io, "gps,gmst")
    for g in gps_times
        gps = lal.LIGOTimeGPS(g)
        val = pyconvert(Float64, lal.GreenwichMeanSiderealTime(gps))
        println(io, "$(repr(g)),$(repr(val))")
    end
end
