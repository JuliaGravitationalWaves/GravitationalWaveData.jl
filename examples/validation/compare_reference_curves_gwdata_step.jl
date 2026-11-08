#!/usr/bin/env julia
# GravitationalWaveData.jl side of compare_reference_curves.jl (see README.md): writes
# `aligo_design_psd`/`avirgo_design_psd` over a frequency grid to a CSV.
#
# Usage: julia --project=<GravitationalWaveData.jl root> compare_reference_curves_gwdata_step.jl <out_csv> <f1> <f2> ...

using GravitationalWaveData: aligo_design_psd, avirgo_design_psd

out_csv = ARGS[1]
freqs = parse.(Float64, ARGS[2:end])

open(out_csv, "w") do io
    println(io, "f,aligo_psd,avirgo_psd")
    for f in freqs
        p_ligo = only(aligo_design_psd([f]))
        p_virgo = only(avirgo_design_psd([f]))
        println(io, "$(repr(f)),$(repr(p_ligo)),$(repr(p_virgo))")
    end
end
