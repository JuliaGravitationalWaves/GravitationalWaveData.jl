#!/usr/bin/env julia
# Compare GravitationalWaveData.jl's H1/L1/V1 detector geometry (src/GravitationalWaveData/Detectors.jl, ported verbatim from
# LALSuite's LALDetectors.h constants) and compute_antenna_response (ported from
# XLALComputeDetAMResponse) against lal.CachedDetectors/lal.ComputeDetAMResponse, across a grid
# of sky positions/polarizations/times including near-pole (near-overhead-or-nadir, detector-
# latitude dependent) and leap-second-boundary edge cases.
#
# Usage (see README.md in this directory):
#   julia examples/validation/compare_antenna_response.jl

include("common.jl")

gwdata_det_csv = output_path("detectors_gwdata.csv")
gwdata_ant_csv = output_path("antenna_gwdata.csv")
lal_det_csv = output_path("detectors_lal.csv")
lal_ant_csv = output_path("antenna_lal.csv")

run_gwdata_step("compare_antenna_response_gwdata_step.jl", gwdata_det_csv, gwdata_ant_csv)
run_reference_step("compare_antenna_response_lal_step.jl", lal_det_csv, lal_ant_csv)

function read_csv(path)
    lines = readlines(path)
    rows = [split(line, ',') for line in lines[2:end]]
    return rows
end

# --- Detector geometry comparison ---
println("=== Detector geometry: location vector + response tensor ===")
gwdata_rows = Dict(r[1] => parse.(Float64, r[2:end]) for r in read_csv(gwdata_det_csv))
lal_rows = Dict(r[1] => parse.(Float64, r[2:end]) for r in read_csv(lal_det_csv))
det_maxabs = 0.0
det_maxrel = 0.0
for name in ("H1", "L1", "V1")
    a = gwdata_rows[name]
    b = lal_rows[name]
    absdiff = abs.(a .- b)
    reldiff = absdiff ./ max.(abs.(b), eps())
    global det_maxabs = max(det_maxabs, maximum(absdiff))
    global det_maxrel = max(det_maxrel, maximum(reldiff))
    println("  $name: max abs diff = $(maximum(absdiff)), max rel diff = $(maximum(reldiff))")
end
println("  overall: max abs = $det_maxabs, max rel = $det_maxrel")
println()

# --- Antenna response comparison ---
println("=== Antenna response F+/Fx across sky positions/times/polarizations ===")
gwdata_ant = read_csv(gwdata_ant_csv)
lal_ant = read_csv(lal_ant_csv)
@assert length(gwdata_ant) == length(lal_ant)

ant_maxabs = 0.0
ant_maxrel = 0.0
for (ga, la) in zip(gwdata_ant, lal_ant)
    @assert ga[1] == la[1] && ga[2] == la[2] && ga[3] == la[3] && ga[4] == la[4] && ga[5] == la[5]
    fp_a, fc_a = parse(Float64, ga[6]), parse(Float64, ga[7])
    fp_b, fc_b = parse(Float64, la[6]), parse(Float64, la[7])
    absdiff = max(abs(fp_a - fp_b), abs(fc_a - fc_b))
    reldiff = max(abs(fp_a - fp_b) / max(abs(fp_b), eps()), abs(fc_a - fc_b) / max(abs(fc_b), eps()))
    global ant_maxabs = max(ant_maxabs, absdiff)
    global ant_maxrel = max(ant_maxrel, reldiff)
    println("  det=$(ga[1]) ra=$(ga[2]) dec=$(ga[3]) psi=$(ga[4]) gps=$(ga[5])  ",
        "GravitationalWaveData=($fp_a,$fc_a)  lal=($fp_b,$fc_b)  abs_diff=$absdiff")
end
println()
println("max absolute F+/Fx difference: $ant_maxabs")
println("max relative F+/Fx difference: $ant_maxrel")
