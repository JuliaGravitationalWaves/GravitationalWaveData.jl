#!/usr/bin/env julia
# Compare GravitationalWaveData.jl's time_delay_from_geocenter (src/GravitationalWaveData/Detectors.jl, ported from LALSuite's
# XLALTimeDelayFromEarthCenter via its XLALArrivalTimeDiff helper) against
# lal.TimeDelayFromEarthCenter, across a grid of sky positions/detectors/times including
# near-pole and leap-second-boundary edge cases.
#
# Usage (see README.md in this directory):
#   julia examples/validation/compare_time_delay.jl

include("common.jl")

gwdata_csv = output_path("time_delay_gwdata.csv")
lal_csv = output_path("time_delay_lal.csv")

run_gwdata_step("compare_time_delay_gwdata_step.jl", gwdata_csv)
run_reference_step("compare_time_delay_lal_step.jl", lal_csv)

function read_csv(path)
    lines = readlines(path)
    rows = [split(line, ',') for line in lines[2:end]]
    return rows
end

println("=== time_delay_from_geocenter vs lal.TimeDelayFromEarthCenter ===")
gwdata_rows = read_csv(gwdata_csv)
lal_rows = read_csv(lal_csv)
@assert length(gwdata_rows) == length(lal_rows)

maxabs = 0.0
maxrel = 0.0
for (ga, la) in zip(gwdata_rows, lal_rows)
    @assert ga[1] == la[1] && ga[2] == la[2] && ga[3] == la[3] && ga[4] == la[4]
    a = parse(Float64, ga[5])
    b = parse(Float64, la[5])
    absdiff = abs(a - b)
    reldiff = absdiff / max(abs(b), eps())
    global maxabs = max(maxabs, absdiff)
    global maxrel = max(maxrel, reldiff)
    println("  det=$(ga[1]) ra=$(ga[2]) dec=$(ga[3]) gps=$(ga[4])  GravitationalWaveData=$a  lal=$b  abs_diff=$absdiff")
end
println()
println("max absolute time-delay difference: $maxabs")
println("max relative time-delay difference: $maxrel")
