#!/usr/bin/env julia
# Compare GWData.jl's aligo_design_psd/avirgo_design_psd (src/GWData/ReferenceCurves.jl, ported
# from LALSimulation's XLALSimNoisePSDaLIGOZeroDetHighPower/XLALSimNoisePSDAdvVirgo) against
# lalsimulation.SimNoisePSDaLIGOZeroDetHighPower/lalsimulation.SimNoisePSDAdvVirgo, across a
# frequency grid spanning the band these curves are meant for (10-2048 Hz), plus a couple of
# points below the aLIGO model's documented ~9 Hz validity floor to see how it behaves there too.
#
# Usage (see README.md in this directory):
#   julia examples/validation/compare_reference_curves.jl

include("common.jl")

# 10-2048 Hz (the band these design curves are meant to be used over) on a log-ish grid, plus a
# couple of points below the aLIGO model's documented ~9 Hz validity floor.
const FREQS = vcat([5.0, 7.0], collect(10.0:10.0:200.0), collect(200.0:25.0:1000.0),
                    collect(1000.0:100.0:2048.0))

gwdata_csv = output_path("reference_curves_gwdata.csv")
lal_csv = output_path("reference_curves_lal.csv")

freq_args = string.(FREQS)

run_gwdata_step("compare_reference_curves_gwdata_step.jl", gwdata_csv, freq_args...)
run_reference_step("compare_reference_curves_lal_step.jl", lal_csv, freq_args...)

function read_csv(path)
    lines = readlines(path)
    f = Float64[]
    aligo = Float64[]
    avirgo = Float64[]
    for line in lines[2:end]
        a, b, c = split(line, ',')
        push!(f, parse(Float64, a))
        push!(aligo, parse(Float64, b))
        push!(avirgo, parse(Float64, c))
    end
    return f, aligo, avirgo
end

gwdata_f, gwdata_aligo, gwdata_avirgo = read_csv(gwdata_csv)
lal_f, lal_aligo, lal_avirgo = read_csv(lal_csv)

@assert gwdata_f == lal_f "Frequency lists don't match"

function report(label, f, a, b)
    println("=== $label ===")
    maxabs = 0.0
    maxrel = 0.0
    for i in eachindex(f)
        absdiff = abs(a[i] - b[i])
        reldiff = absdiff / abs(b[i])
        maxabs = max(maxabs, absdiff)
        maxrel = max(maxrel, reldiff)
        println("  f=$(f[i])  GWData=$(a[i])  lal=$(b[i])  rel_diff=$reldiff")
    end
    println("  max absolute difference: $maxabs")
    println("  max relative difference: $maxrel")
    println()
    return maxrel
end

println("Comparing GWData.jl's reference design PSDs vs lalsimulation, $(length(FREQS)) frequencies")
println()
rel_ligo = report("aLIGO zero-detuning high-power PSD", gwdata_f, gwdata_aligo, lal_aligo)
rel_virgo = report("AdvVirgo PSD", gwdata_f, gwdata_avirgo, lal_avirgo)

println("Overall max relative difference: aLIGO=$rel_ligo  AdvVirgo=$rel_virgo")
