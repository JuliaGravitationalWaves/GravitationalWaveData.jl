#!/usr/bin/env julia
# Cross-validate GravitationalWaveData.jl's bandpass/highpass/lowpass/whiten against gwpy on the same real
# GW150914 H1 strain segment.
#
# bandpass/highpass/lowpass caveat (see Filtering.jl's module docstring for the full rationale):
# gwpy's own default `.bandpass()`/`.highpass()`/`.lowpass()` design an *auto-order* Chebyshev-I
# filter via `scipy.signal.iirdesign` to hit explicit gpass/gstop specs; this package instead uses
# a *fixed-order* Butterworth filter. The gwpy step forces `ftype="butter"` so both sides are at
# least the same filter family, but the order is still auto-selected on gwpy's side vs. fixed
# (order=4) here, so exact sample-for-sample agreement is not expected or the right thing to
# check. What IS checked: agreement in each filter's unambiguous, order-independent region -- e.g.
# bandpass's *interior passband* should match the original unfiltered signal on both sides.
#
# whiten has no such caveat: this package's algorithm is a direct port of gwpy's own
# `TimeSeries.whiten()`/`.convolve()`/`filter_design.fir_from_transfer` (see Filtering.jl's
# docstrings), and both sides are given the *exact same* ASD (computed once by the gwpy step and
# shared via `asd.csv`) so PSD-estimation differences can't be a confound -- so whiten's agreement
# is expected to be tight, the same order as compare_psd.jl's ~0.1-1% PSD agreement.
#
# Usage (see README.md in this directory):
#   julia examples/validation/compare_filtering.jl

include("common.jl")

const H1_PATH = gw150914_strain_path("H1")

println("Step 1/2: gwpy (asd + bandpass/highpass/lowpass/whiten)...")
run_reference_step("compare_filtering_gwpy_step.jl", H1_PATH, OUT_DIR)

println("Step 2/2: GravitationalWaveData.jl (bandpass/highpass/lowpass/whiten)...")
run_gwdata_step("compare_filtering_gwdata_step.jl", H1_PATH, OUT_DIR)

function read_csv(path)
    lines = readlines(path)
    xs = Float64[]
    ys = Float64[]
    for line in lines[2:end]
        x, y = split(line, ',')
        push!(xs, parse(Float64, x))
        push!(ys, parse(Float64, y))
    end
    return xs, ys
end

"""
    aligned(gt, gy, jt, jy) -> (gy_aligned, jy_aligned)

Align two (time, value) series sampled at the same `dt` but each independently cropped/rounded to
the nearest sample by gwpy vs. this package -- shifts by the (small, usually 0) integer sample
offset between `gt[1]` and `jt[1]` and truncates to the overlapping length, so index-for-index
comparison below is meaningful even if the two crops picked adjacent samples.
"""
function aligned(gt, gy, jt, jy)
    dt = gt[2] - gt[1]
    shift = round(Int, (jt[1] - gt[1]) / dt)   # julia's t0 minus gwpy's t0, in samples
    if shift >= 0
        g = gy[(1 + shift):end]
        j = jy
    else
        g = gy
        j = jy[(1 - shift):end]
    end
    n = min(length(g), length(j))
    return g[1:n], j[1:n]
end

rel_rms(a, b) = sqrt(sum(abs2, a .- b) / length(a)) / sqrt(sum(abs2, b) / length(b))

println()
println("=== bandpass/highpass/lowpass (order=4 Butterworth, both sides; gwpy forced to ftype=\"butter\") ===")
for name in ("bandpass", "highpass", "lowpass")
    gt, gy = read_csv(joinpath(OUT_DIR, "gwpy_$name.csv"))
    jt, jy = read_csv(joinpath(OUT_DIR, "julia_$name.csv"))
    ga, ja = aligned(gt, gy, jt, jy)
    fs = round(Int, 1 / (gt[2] - gt[1]))
    edge = fs   # drop 1s from each end (filtfilt edge-transient region)
    core_g = ga[(edge + 1):(end - edge)]
    core_j = ja[(edge + 1):(end - edge)]
    println("  $name: n=$(length(ga)) (full), $(length(core_g)) (interior, 1s edges dropped)")
    println("    whole-segment relative RMS diff (julia vs gwpy): $(rel_rms(ja, ga))")
    println("    interior      relative RMS diff (julia vs gwpy): $(rel_rms(core_j, core_g))")
end

println()
println("=== whiten (fduration=2.0, same off-source ASD fed to both sides) ===")
gt, gy = read_csv(joinpath(OUT_DIR, "gwpy_whiten.csv"))
jt, jy = read_csv(joinpath(OUT_DIR, "julia_whiten.csv"))
ga, ja = aligned(gt, gy, jt, jy)
fs = round(Int, 1 / (gt[2] - gt[1]))
edge = round(Int, 1.0 * fs)   # gwpy's own documented fduration/2 = 1.0s edge-corruption region
core_g = ga[(edge + 1):(end - edge)]
core_j = ja[(edge + 1):(end - edge)]
abs_diff = abs.(core_j .- core_g)
denom = abs.(core_g)
rel_diff = abs_diff ./ max.(denom, 1e-3 * sqrt(sum(abs2, core_g) / length(core_g)))
println("  n=$(length(ga)) (full), $(length(core_g)) (interior, 1s edges dropped each side)")
println("  interior relative RMS diff (julia vs gwpy): $(rel_rms(core_j, core_g))")
println("  interior |diff| : median=$(sort(rel_diff)[length(rel_diff) ÷ 2])  max=$(maximum(rel_diff))")
println("  interior RMS: julia=$(sqrt(sum(abs2,core_j)/length(core_j)))  gwpy=$(sqrt(sum(abs2,core_g)/length(core_g)))")
