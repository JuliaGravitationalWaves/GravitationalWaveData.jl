#!/usr/bin/env julia
# Reference side of compare_coherence.jl (see README.md): computes gwpy's
# `TimeSeries.csd()`/`.coherence()` over the same windows and writes both to a CSV. Runs in this
# directory's project and must not load GWData.
#
# Usage: julia --project=examples/validation compare_coherence_gwpy_step.jl <h1_path> <l1_path> <off_start> <off_end> <duration> <out_csv>

using PythonCall

h1_path, l1_path, off_start_s, off_end_s, duration_s, out_csv = ARGS
off_source_start = parse(Float64, off_start_s)
off_source_end = parse(Float64, off_end_s)
duration = parse(Float64, duration_s)

gwpy_ts_mod = pyimport("gwpy.timeseries")
TimeSeries = gwpy_ts_mod.TimeSeries

py_h1 = TimeSeries.read(h1_path, format="hdf5.gwosc")
py_l1 = TimeSeries.read(l1_path, format="hdf5.gwosc")

py_h1_off = py_h1.crop(off_source_start, off_source_end)
py_l1_off = py_l1.crop(off_source_start, off_source_end)

py_csd = py_h1_off.csd(py_l1_off, fftlength=duration, overlap=duration * 0.5, window="hann")
py_coh = py_h1_off.coherence(py_l1_off, fftlength=duration, overlap=duration * 0.5, window="hann")

freqs = pyconvert(Vector{Float64}, py_csd.frequencies.value)
csd_vals = pyconvert(Vector{ComplexF64}, py_csd.value)
coh_vals = pyconvert(Vector{Float64}, py_coh.value)

open(out_csv, "w") do io
    println(io, "freq,csd_re,csd_im,coherence")
    for (f, c, k) in zip(freqs, csd_vals, coh_vals)
        println(io, "$f,$(real(c)),$(imag(c)),$k")
    end
end
