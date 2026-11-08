#!/usr/bin/env julia
# Reference side of compare_psd.jl (see README.md): computes gwpy's median-Welch
# PSD over the same off-source window and writes it to a CSV. Runs in this
# directory's project and must not load GravitationalWaveData.
#
# Usage: julia --project=examples/validation compare_psd_gwpy_step.jl <h1_path> <off_start> <off_end> <duration> <out_csv>

using PythonCall

h1_path, off_start_s, off_end_s, duration_s, out_csv = ARGS
off_source_start = parse(Float64, off_start_s)
off_source_end = parse(Float64, off_end_s)
duration = parse(Float64, duration_s)

gwpy_ts_mod = pyimport("gwpy.timeseries")
TimeSeries = gwpy_ts_mod.TimeSeries

py_ts = TimeSeries.read(h1_path, format="hdf5.gwosc")
py_off_source = py_ts.crop(off_source_start, off_source_end)
py_psd = py_off_source.psd(fftlength=duration, overlap=duration * 0.5, method="median", window="hann")

freqs = pyconvert(Vector{Float64}, py_psd.frequencies.value)
psd = pyconvert(Vector{Float64}, py_psd.value)

open(out_csv, "w") do io
    println(io, "freq,psd")
    for (f, p) in zip(freqs, psd)
        println(io, "$f,$p")
    end
end
