#!/usr/bin/env julia
# Reference side of compare_filtering.jl (see README.md): computes, with gwpy: Runs in this
# directory's project and must not load GravitationalWaveData.
#   - an off-source ASD (median-Welch, same convention as compare_psd.jl) -> asd.csv, so the
#     GravitationalWaveData step whitens against the *exact same* numbers and PSD-estimation differences
#     can't be confused with whitening-algorithm differences.
#   - `.bandpass()`/`.highpass()`/`.lowpass()`, forced to `ftype="butter"` so both sides use
#     the same filter *family* (gwpy's default is an auto-order Chebyshev-I via
#     `scipy.signal.iirdesign`; see Filtering.jl's module docstring), on a raw analysis segment.
#   - `.whiten(asd=<the same off-source ASD>)` on a shorter on-source segment.
#
# Usage: julia --project=examples/validation compare_filtering_gwpy_step.jl <h1_path> <out_dir>

using PythonCall

h1_path, out_dir = ARGS

gwpy_ts_mod = pyimport("gwpy.timeseries")
GwpyTimeSeries = gwpy_ts_mod.TimeSeries

py_ts = GwpyTimeSeries.read(h1_path, format="hdf5.gwosc")

const TRIGGER_TIME = 1126259462.4
const ANALYSIS_START = TRIGGER_TIME - 2.0   # matches compare_psd.jl's (duration=4, post_trigger=2) convention
const PSD_SEGMENT = 1024.0
const FFTLENGTH = 4.0

off_source = py_ts.crop(ANALYSIS_START - PSD_SEGMENT, ANALYSIS_START)
py_asd = off_source.asd(fftlength=FFTLENGTH, overlap=FFTLENGTH / 2, method="median", window="hann")

function write_csv(path, freqs, vals)
    open(path, "w") do io
        println(io, "x,y")
        for (f, v) in zip(freqs, vals)
            println(io, "$f,$v")
        end
    end
end

asd_freqs = pyconvert(Vector{Float64}, py_asd.frequencies.value)
asd_vals = pyconvert(Vector{Float64}, py_asd.value)
write_csv(joinpath(out_dir, "asd.csv"), asd_freqs, asd_vals)

# --- bandpass / highpass / lowpass, on a longer raw (unwindowed) segment ------------------------
const FILTER_SEG_START = TRIGGER_TIME - 16.0
const FILTER_SEG_END = TRIGGER_TIME + 16.0
filter_seg = py_ts.crop(FILTER_SEG_START, FILTER_SEG_END)
t0 = pyconvert(Float64, filter_seg.t0.value)
dt = pyconvert(Float64, filter_seg.dt.value)
raw_vals = pyconvert(Vector{Float64}, filter_seg.value)
n = length(raw_vals)
times = t0 .+ dt .* (0:(n - 1))
write_csv(joinpath(out_dir, "raw_segment.csv"), times, raw_vals)

bp = filter_seg.bandpass(35.0, 350.0, ftype="butter")
hp = filter_seg.highpass(35.0, ftype="butter")
lp = filter_seg.lowpass(350.0, ftype="butter")

write_csv(joinpath(out_dir, "gwpy_bandpass.csv"), times, pyconvert(Vector{Float64}, bp.value))
write_csv(joinpath(out_dir, "gwpy_highpass.csv"), times, pyconvert(Vector{Float64}, hp.value))
write_csv(joinpath(out_dir, "gwpy_lowpass.csv"), times, pyconvert(Vector{Float64}, lp.value))

# --- whiten, on a shorter on-source segment, against the off-source ASD above --------------------
const WHITEN_SEG_START = TRIGGER_TIME - 6.0
const WHITEN_SEG_END = TRIGGER_TIME + 2.0
whiten_seg = py_ts.crop(WHITEN_SEG_START, WHITEN_SEG_END)
wt0 = pyconvert(Float64, whiten_seg.t0.value)
wdt = pyconvert(Float64, whiten_seg.dt.value)
wn = pyconvert(Int, whiten_seg.size)
wtimes = wt0 .+ wdt .* (0:(wn - 1))

whitened = whiten_seg.whiten(asd=py_asd, fduration=2.0)
write_csv(joinpath(out_dir, "gwpy_whiten.csv"), wtimes, pyconvert(Vector{Float64}, whitened.value))

println("gwpy step done: wrote asd.csv, raw_segment.csv, gwpy_{bandpass,highpass,lowpass,whiten}.csv to $out_dir")
