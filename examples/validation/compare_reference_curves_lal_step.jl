#!/usr/bin/env julia
# Reference side of compare_reference_curves.jl (see README.md): writes
# lalsimulation.SimNoisePSDaLIGOZeroDetHighPower and SimNoisePSDAdvVirgo over the same
# frequency grid to a CSV. Runs in this
# directory's project and must not load GWData.
#
# Usage: julia --project=examples/validation compare_reference_curves_lal_step.jl <out_csv> <f1> <f2> ...

using PythonCall
const lalsim = pyimport("lalsimulation")

out_csv = ARGS[1]
freqs = parse.(Float64, ARGS[2:end])

open(out_csv, "w") do io
    println(io, "f,aligo_psd,avirgo_psd")
    for f in freqs
        p_ligo = pyconvert(Float64, lalsim.SimNoisePSDaLIGOZeroDetHighPower(f))
        p_virgo = pyconvert(Float64, lalsim.SimNoisePSDAdvVirgo(f))
        println(io, "$(repr(f)),$(repr(p_ligo)),$(repr(p_virgo))")
    end
end
