#!/usr/bin/env julia
# Reference side of compare_qtransform.jl (see README.md): computes gwpy's raw
# (pre-interpolation) Q-scan peak over the same segment and writes it to a CSV. Runs in this
# directory's project and must not load GravitationalWaveData.
#
# Calls `gwpy.signal.qtransform.q_scan` directly (the function `TimeSeries.q_transform()` calls
# internally) and reads its raw `QGram.peak`, not the interpolated `Spectrogram`'s argmax pixel.
# GravitationalWaveData.jl's `q_transform` stops at the same raw QGram/peak (the spline regridding is
# presentation-only; see src/GravitationalWaveData/QTransform.jl's module docstring), so this is the
# like-for-like target.
#
# Usage: julia --project=examples/validation compare_qtransform_gwpy_step.jl <h1_path> <seg_start> <seg_end> <fmin> <fmax> <out_csv>

using PythonCall

h1_path, seg_start_s, seg_end_s, fmin_s, fmax_s, out_csv = ARGS
seg_start = parse(Float64, seg_start_s)
seg_end = parse(Float64, seg_end_s)
fmin = parse(Float64, fmin_s)
fmax = parse(Float64, fmax_s)

gwpy_ts_mod = pyimport("gwpy.timeseries")
qmod = pyimport("gwpy.signal.qtransform")
TimeSeries = gwpy_ts_mod.TimeSeries

py_ts = TimeSeries.read(h1_path, format="hdf5.gwosc")
segment = py_ts.crop(seg_start, seg_end)

qgram_far = qmod.q_scan(segment; mismatch=0.2, qrange=(4.0, 64.0), frange=(fmin, fmax))
qgram = qgram_far[0]
peak = qgram.peak

time = pyconvert(Float64, peak["time"])
frequency = pyconvert(Float64, peak["frequency"])
energy = pyconvert(Float64, peak["energy"])
snr = pyconvert(Float64, peak["snr"])
q = pyconvert(Float64, qgram.plane.q)

open(out_csv, "w") do io
    println(io, "field,value")
    println(io, "time,$time")
    println(io, "frequency,$frequency")
    println(io, "energy,$energy")
    println(io, "snr,$snr")
    println(io, "q,$q")
end
