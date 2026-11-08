#!/usr/bin/env julia
# Compare GWData.jl's `csd`/`coherence` against gwpy's `TimeSeries.csd()`/`.coherence()`, using a
# genuine two-channel case: real H1 and L1 GW150914 strain over the same off-source GPS window.
# A synthetic or self-coherence case would be a much weaker check.
#
# Usage (see README.md in this directory):
#   julia examples/validation/compare_coherence.jl

include("common.jl")

const H1_PATH = gw150914_strain_path("H1")
const L1_PATH = gw150914_strain_path("L1")

# Same off-source window compare_psd.jl uses (trigger_time, duration=4.0, post_trigger=2.0):
# the 1024s immediately before the 4s analysis window, well clear of the GW150914 signal itself,
# so H1/L1 are expected to be *incoherent* detector noise here (a genuinely non-trivial, real
# check -- unlike a signal segment where astrophysical + noise coherence would be entangled).
const TRIGGER_TIME = 1126259462.4
const DURATION = 4.0
const POST_TRIGGER = 2.0
const CSD_SEGMENT = 1024.0

analysis_start = TRIGGER_TIME - (DURATION - POST_TRIGGER)
off_source_end = analysis_start
off_source_start = off_source_end - CSD_SEGMENT

gwdata_csv = output_path("coherence_gwdata.csv")
gwpy_csv = output_path("coherence_gwpy.csv")

run_gwdata_step("compare_coherence_gwdata_step.jl", H1_PATH, L1_PATH, off_source_start, off_source_end, DURATION, gwdata_csv)
run_reference_step("compare_coherence_gwpy_step.jl", H1_PATH, L1_PATH, off_source_start, off_source_end, DURATION, gwpy_csv)

function read_csv(path)
    lines = readlines(path)
    freqs = Float64[]
    re = Float64[]
    im = Float64[]
    coh = Float64[]
    for line in lines[2:end]
        f, r, i, c = split(line, ',')
        push!(freqs, parse(Float64, f))
        push!(re, parse(Float64, r))
        push!(im, parse(Float64, i))
        push!(coh, parse(Float64, c))
    end
    return freqs, re, im, coh
end

gwdata_freqs, gwdata_re, gwdata_im, gwdata_coh = read_csv(gwdata_csv)
gwpy_freqs, gwpy_re, gwpy_im, gwpy_coh = read_csv(gwpy_csv)

@assert length(gwpy_freqs) == length(gwdata_freqs) "Frequency grid length mismatch: gwpy=$(length(gwpy_freqs)) GWData=$(length(gwdata_freqs))"
@assert isapprox(gwpy_freqs, gwdata_freqs; atol=1e-6) "Frequency grids don't match bin-for-bin"

band = findall(f -> 20.0 <= f <= 1024.0, gwdata_freqs)

gwdata_csd_mag = sqrt.(gwdata_re[band] .^ 2 .+ gwdata_im[band] .^ 2)
gwpy_csd_mag = sqrt.(gwpy_re[band] .^ 2 .+ gwpy_im[band] .^ 2)
csd_rel_diff = abs.(gwdata_csd_mag .- gwpy_csd_mag) ./ gwpy_csd_mag

coh_abs_diff = abs.(gwdata_coh[band] .- gwpy_coh[band])

println("Comparing GWData.jl's `csd`/`coherence` vs gwpy TimeSeries.csd/coherence, real H1 vs L1")
println("  off-source window: GPS [$off_source_start, $off_source_end) ($(CSD_SEGMENT)s)")
println("  seglen = $DURATION s, overlap = 50%, window = hann, average = mean")
println("  frequency points compared: $(length(band)) (20-1024 Hz band)")
println()
println("CSD magnitude relative difference (in-band):")
println("  median = $(sort(csd_rel_diff)[length(csd_rel_diff)÷2])")
println("  max    = $(maximum(csd_rel_diff))")
println()
println("Coherence absolute difference (in-band, coherence itself is in [0,1]):")
println("  median = $(sort(coh_abs_diff)[length(coh_abs_diff)÷2])")
println("  max    = $(maximum(coh_abs_diff))")
println("  mean gwpy coherence in-band   = $(sum(gwpy_coh[band]) / length(band))")
println("  mean GWData coherence in-band = $(sum(gwdata_coh[band]) / length(band))")
println()
for f in (20.0, 50.0, 100.0, 150.0, 200.0, 300.0, 500.0, 1000.0)
    idx = argmin(abs.(gwdata_freqs .- f))
    println("  f=$(gwdata_freqs[idx]) Hz  " *
            "GWData csd=$(gwdata_re[idx])+$(gwdata_im[idx])im coh=$(gwdata_coh[idx])  " *
            "gwpy csd=$(gwpy_re[idx])+$(gwpy_im[idx])im coh=$(gwpy_coh[idx])")
end
