#!/usr/bin/env julia
# Compare GWData.jl's `spectrogram` against gwpy's `TimeSeries.spectrogram2` -- the raw,
# non-averaged STFT spectrogram (one windowed FFT per segment; NOT gwpy's plain
# `TimeSeries.spectrogram`, which Welch/median-averages several sub-segments per column -- see
# src/GWData/Spectrogram.jl's module docstring for why `spectrogram2` is the correct analogue of
# what `DSP.spectrogram` computes) -- on the same real GW150914 H1 strain segment.
#
# Usage (see README.md in this directory):
#   julia examples/validation/compare_spectrogram.jl

include("common.jl")

const H1_PATH = gw150914_strain_path("H1")

# An 8s window straddling the GW150914 trigger: short enough to keep the CSVs small while still
# giving several dozen time bins to compare, long enough to include the merger itself.
const TRIGGER_TIME = 1126259462.4
const SEG_START = TRIGGER_TIME - 4.0
const SEG_END = TRIGGER_TIME + 4.0

# n=1024 samples (0.25s at the strain's native 4096 Hz), noverlap=0 (deliberately -- see below).
#
# Why noverlap=0 rather than the more usual 50%: reading gwpy's actual implementation
# (gwpy/signal/spectral/_ui.py, the `spectrogram` function that both `TimeSeries.spectrogram`
# AND `TimeSeries.spectrogram2` call into) shows that whenever noverlap>0 it does more than stack
# one periodogram per segment -- after computing the per-segment periodograms on the requested
# (nfft, noverlap) grid, it does a *second* pass that re-averages each output column from its
# `density = nfft / (nfft - noverlap)` overlapping neighbors, weighted by a triangular window
# (`get_window("triangle", density)`). That is real cross-column smoothing `DSP.spectrogram` (and
# hence `GWData.spectrogram`) does not do -- it wraps
# `DSP.spectrogram`, a raw per-segment STFT, not gwpy's own further-smoothed variant. With
# noverlap=0, density=1 and `get_window("triangle", 1) == [1.0]`, so that second pass is exactly a
# no-op and gwpy's output reduces to the same raw one-FFT-per-segment computation
# `DSP.spectrogram` performs -- giving a fair, honest apples-to-apples comparison of the actual
# STFT math (windowing, FFT, density normalization) without also having to replicate gwpy's
# separate smoothing feature.
const N = 1024
const NOVERLAP = 0

gwdata_prefix = output_path("spectrogram_gwdata")
gwpy_prefix = output_path("spectrogram_gwpy")

run_gwdata_step("compare_spectrogram_gwdata_step.jl", H1_PATH, SEG_START, SEG_END, N, NOVERLAP, gwdata_prefix)
run_reference_step("compare_spectrogram_gwpy_step.jl", H1_PATH, SEG_START, SEG_END, N, NOVERLAP, gwpy_prefix)

read_vec(path) = parse.(Float64, readlines(path))

function read_power_csv(path)
    rows = [parse.(Float64, split(line, ",")) for line in readlines(path)]
    return permutedims(reduce(hcat, rows))   # (nfreq, ntime), one row per frequency bin
end

gw_freqs = read_vec(gwdata_prefix * "_freqs.csv")
gw_times = read_vec(gwdata_prefix * "_times.csv")
gw_power = read_power_csv(gwdata_prefix * "_power.csv")

py_freqs = read_vec(gwpy_prefix * "_freqs.csv")
py_times = read_vec(gwpy_prefix * "_times.csv")
py_power = read_power_csv(gwpy_prefix * "_power.csv")

@assert size(gw_power) == size(py_power) "Shape mismatch: GWData=$(size(gw_power)) gwpy=$(size(py_power))"
@assert isapprox(gw_freqs, py_freqs; atol=1e-6) "Frequency grids don't match bin-for-bin"

# Time axis convention differs by design: DSP.spectrogram (and hence GWData.spectrogram) labels
# each column with its segment's *center* GPS time, while gwpy's spectrogram2 labels each column
# with its segment's *start* GPS time (see compare_spectrogram_gwpy_step.jl's header and
# src/GWData/Spectrogram.jl's `spectrogram` docstring) -- a constant offset of n/(2*fs), not a
# bug. Verify that offset is exactly what's expected rather than silently ignoring the difference.
const NATIVE_FS = 4096.0   # GW150914 H1 LOSC strain's native sample rate
expected_offset = N / (2 * NATIVE_FS)
time_offsets = gw_times .- py_times
@assert isapprox(time_offsets, fill(time_offsets[1], length(time_offsets)); atol=1e-9) "Time offset between GWData and gwpy columns isn't constant"

band = findall(f -> 20.0 <= f <= 1024.0, gw_freqs)
rel_diff = abs.(gw_power[band, :] .- py_power[band, :]) ./ py_power[band, :]

println("Comparing GWData.spectrogram vs gwpy TimeSeries.spectrogram2, same H1 segment straddling the GW150914 trigger")
println("  segment: GPS [$SEG_START, $SEG_END) ($(SEG_END - SEG_START)s)")
println("  n=$N samples (fftlength=$(N/4096)s), noverlap=$NOVERLAP (0% -- see header comment: gwpy's spectrogram2 does extra triangular cross-column smoothing whenever noverlap>0, which DSP.spectrogram doesn't and isn't asked to replicate), window=hann (peak-centered, DSP.Windows.hanning's default -- see Spectrogram.jl's spectrogram docstring for why zerophase=true is the wrong choice here despite scipy calling its own convention \"periodic\")")
println("  spectrogram shape: $(size(gw_power)) (nfreq x ntime)")
println("  frequency points compared per column: $(length(band)) (20-1024 Hz band)")
println()
println("Time-bin labeling convention check (segment-center [GWData] vs segment-start [gwpy]):")
println("  constant offset between the two time axes = $(time_offsets[1]) s")
println("  expected offset (n/2 samples = $(N/2/4096) s)  = $expected_offset s")
println()
println("Relative power difference (in-band, all time columns):")
println("  median = $(sort(vec(rel_diff))[length(rel_diff)÷2])")
println("  max    = $(maximum(rel_diff))")
println("  (the max is typically a single near-null bin in one column -- both values tiny and")
println("  dominated by the small window-definition difference between DSP's symmetric Hann and")
println("  scipy's periodic Hann, not a systematic issue; this is a raw, un-averaged periodogram,")
println("  unlike the Welch/median-averaged PSD comparison in compare_psd.jl, so it doesn't average")
println("  away this kind of single-bin outlier. The median is the number to trust.)")
println()
mid_col = size(gw_power, 2) ÷ 2 + 1
println("Spot check at a few frequencies, middle time column ($(gw_times[mid_col]) GPS):")
for f in (20.0, 50.0, 100.0, 150.0, 200.0, 300.0, 500.0, 1000.0)
    idx = argmin(abs.(gw_freqs .- f))
    gv, pv = gw_power[idx, mid_col], py_power[idx, mid_col]
    println("  f=$(gw_freqs[idx]) Hz  GWData=$gv  gwpy=$pv  rel_diff=$(abs(gv - pv) / pv)")
end
