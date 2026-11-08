#!/usr/bin/env julia
# Reference side of compare_spectrogram.jl (see README.md): computes gwpy's
# non-averaged power spectrogram (`TimeSeries.spectrogram2`) over the same segment and writes
# its frequency axis, GPS time axis and power matrix to CSVs. Runs in this
# directory's project and must not load GWData.
#
# Usage: julia --project=examples/validation compare_spectrogram_gwpy_step.jl <h1_path> <seg_start> <seg_end> <n> <noverlap> <out_prefix>

using PythonCall

h1_path, seg_start_s, seg_end_s, n_s, noverlap_s, out_prefix = ARGS
seg_start = parse(Float64, seg_start_s)
seg_end = parse(Float64, seg_end_s)
n = parse(Int, n_s)
noverlap = parse(Int, noverlap_s)

gwpy_ts_mod = pyimport("gwpy.timeseries")
TimeSeries = gwpy_ts_mod.TimeSeries

py_ts = TimeSeries.read(h1_path, format="hdf5.gwosc")
py_seg = py_ts.crop(seg_start, seg_end)

fs = pyconvert(Float64, py_seg.sample_rate.to("Hz").value)
fftlength = n / fs
overlap = noverlap / fs

# `detrend=false`: GWData.spectrogram (via DSP.spectrogram) does not remove each segment's mean
# before the FFT, but scipy.signal.periodogram (which spectrogram2 calls per-segment) defaults to
# `detrend="constant"` (mean removal) -- match GWData's convention explicitly rather than
# comparing two different conventions. `window="hann"`: gwpy's TimeSeries.spectrogram2 default;
# scipy.signal.get_window peak-centers this window in the middle of the segment, matching
# GWData.spectrogram's own default (DSP.Windows.hanning(n), i.e. `zerophase=false` -- NOT
# `zerophase=true`, which centers the window's peak at index 1 instead, a different
# zero-phase-filter-kernel convention; see Spectrogram.jl's `spectrogram` docstring for how this
# was found by numerically cross-checking a single segment against a from-scratch reference).
py_spec = py_seg.spectrogram2(fftlength=fftlength, overlap=overlap, window="hann", detrend=false)

freqs = pyconvert(Vector{Float64}, py_spec.frequencies.value)
gwpy_times = pyconvert(Vector{Float64}, py_spec.times.value)
power = pyconvert(Matrix{Float64}, py_spec.value)  # shape (ntime, nfreq)

open(out_prefix * "_freqs.csv", "w") do io
    for f in freqs
        println(io, f)
    end
end

open(out_prefix * "_times.csv", "w") do io
    for tt in gwpy_times
        println(io, tt)
    end
end

open(out_prefix * "_power.csv", "w") do io
    nt, nf = size(power)
    for i in 1:nf
        println(io, join(power[:, i], ","))
    end
end
