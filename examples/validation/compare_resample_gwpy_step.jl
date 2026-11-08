#!/usr/bin/env julia
# Reference side of compare_resample.jl (see README.md): applies gwpy's
# `.resample()` to the same segment and writes its plain Welch PSD (`method="welch"`) to a CSV. Runs in this
# directory's project and must not load GravitationalWaveData.
#
# Usage: julia --project=examples/validation compare_resample_gwpy_step.jl <h1_path> <seg_start> <seg_end> <new_rate> <out_csv>

using PythonCall

h1_path, start_s, end_s, new_rate_s, out_csv = ARGS
seg_start = parse(Float64, start_s)
seg_end = parse(Float64, end_s)
new_rate = parse(Float64, new_rate_s)

gwpy_ts_mod = pyimport("gwpy.timeseries")
TimeSeries = gwpy_ts_mod.TimeSeries

py_ts = TimeSeries.read(h1_path, format="hdf5.gwosc")
py_seg = py_ts.crop(seg_start, seg_end)
py_resampled = py_seg.resample(new_rate)
py_psd = py_resampled.psd(fftlength=4.0, overlap=2.0, method="welch", window="hann")

freqs = pyconvert(Vector{Float64}, py_psd.frequencies.value)
psd = pyconvert(Vector{Float64}, py_psd.value)

open(out_csv, "w") do io
    println(io, "freq,psd")
    for (f, p) in zip(freqs, psd)
        println(io, "$f,$p")
    end
end
