#!/usr/bin/env julia
# Compare GWData.jl's `gmst` (ported from LALSuite's XLALGreenwichMeanSiderealTime, see
# src/GWData/SiderealTime.jl) against lal.GreenwichMeanSiderealTime itself, across a range of GPS
# times chosen to exercise: the J2000.0 reference epoch, a leap-second boundary (2015-Jul-01,
# GPS 1119744016 -- one second before/at/after), real GW event trigger times spanning O1-O4, and
# a time carrying a fractional second.
#
# Usage (see README.md in this directory):
#   julia examples/validation/compare_gmst.jl

include("common.jl")

const GPS_TIMES = [
    630763213.0,     # J2000.0 epoch (2000-01-01 12:00:00 UTC) -- XLAL_EPOCH_J2000_0_GPS
    1119744015.0,    # 1s before the 2015-Jul-01 leap second
    1119744016.0,    # exactly at the 2015-Jul-01 leap second (TAI-UTC 35->36)
    1119744017.0,    # 1s after
    1126259462.4,    # GW150914 trigger time (fractional second)
    1187008882.4,    # GW170817 trigger time
    1238166018.5,    # GW190425-ish O3 time
    1440000000.0,    # a "recent" O4-era-ish time, no further leap seconds in LALSuite's table
                      # since 2017 -- matches reality, so this exercises the "beyond the last
                      # table entry, hold the last value" branch consistently
]

gwdata_csv = output_path("gmst_gwdata.csv")
lal_csv = output_path("gmst_lal.csv")

gps_args = string.(GPS_TIMES)

run_gwdata_step("compare_gmst_gwdata_step.jl", gwdata_csv, gps_args...)
run_reference_step("compare_gmst_lal_step.jl", lal_csv, gps_args...)

function read_csv(path)
    lines = readlines(path)
    gps = Float64[]
    val = Float64[]
    for line in lines[2:end]
        g, v = split(line, ',')
        push!(gps, parse(Float64, g))
        push!(val, parse(Float64, v))
    end
    return gps, val
end

gwdata_gps, gwdata_gmst = read_csv(gwdata_csv)
lal_gps, lal_gmst = read_csv(lal_csv)

@assert gwdata_gps == lal_gps "GPS time lists don't match"

println("Comparing GWData.jl's gmst() vs lal.GreenwichMeanSiderealTime, $(length(GPS_TIMES)) GPS times")
println()
maxabs = 0.0
maxrel = 0.0
for i in eachindex(gwdata_gps)
    g = gwdata_gps[i]
    a = gwdata_gmst[i]
    b = lal_gmst[i]
    absdiff = abs(a - b)
    reldiff = absdiff / abs(b)
    global maxabs = max(maxabs, absdiff)
    global maxrel = max(maxrel, reldiff)
    println("  gps=$g  GWData=$a  lal=$b  abs_diff=$absdiff  rel_diff=$reldiff")
end
println()
println("max absolute difference (radians): $maxabs")
println("max relative difference: $maxrel")
