#!/usr/bin/env julia
# Compare GWData.jl's `psd` against gwpy's median-Welch PSD on the same real GW150914 off-source
# H1 strain segment.
#
# Usage (see README.md in this directory):
#   julia examples/validation/compare_psd.jl

include("common.jl")

const H1_PATH = gw150914_strain_path("H1")

# The 1024 s immediately preceding a 4 s analysis window ending 2 s after the GW150914 trigger
# (the usual off-source PSD window for parameter estimation).
const TRIGGER_TIME = 1126259462.4
const DURATION = 4.0
const POST_TRIGGER = 2.0
const PSD_SEGMENT = 1024.0

analysis_start = TRIGGER_TIME - (DURATION - POST_TRIGGER)
off_source_end = analysis_start
off_source_start = off_source_end - PSD_SEGMENT

gwdata_csv = output_path("psd_gwdata.csv")
gwpy_csv = output_path("psd_gwpy.csv")

run_gwdata_step("compare_psd_gwdata_step.jl", H1_PATH, off_source_start, off_source_end, DURATION, gwdata_csv)
run_reference_step("compare_psd_gwpy_step.jl", H1_PATH, off_source_start, off_source_end, DURATION, gwpy_csv)

function read_psd_csv(path)
    lines = readlines(path)
    freqs = Float64[]
    psd = Float64[]
    for line in lines[2:end]
        f, p = split(line, ',')
        push!(freqs, parse(Float64, f))
        push!(psd, parse(Float64, p))
    end
    return freqs, psd
end

gwdata_freqs, gwdata_psd = read_psd_csv(gwdata_csv)
gwpy_freqs, gwpy_psd = read_psd_csv(gwpy_csv)

@assert length(gwpy_freqs) == length(gwdata_freqs) "Frequency grid length mismatch: gwpy=$(length(gwpy_freqs)) GWData=$(length(gwdata_freqs))"
@assert isapprox(gwpy_freqs, gwdata_freqs; atol=1e-6) "Frequency grids don't match bin-for-bin"

band = findall(f -> 20.0 <= f <= 1024.0, gwdata_freqs)
rel_diff = abs.(gwdata_psd[band] .- gwpy_psd[band]) ./ gwpy_psd[band]

println("Comparing GWData.jl's `psd` vs gwpy TimeSeries.psd, same H1 off-source segment")
println("  off-source window: GPS [$off_source_start, $off_source_end) ($(PSD_SEGMENT)s)")
println("  seglen = $DURATION s, overlap = 50%, method = median, window = hann")
println("  frequency points compared: $(length(band)) (20-1024 Hz band)")
println()
println("Relative PSD difference (in-band):")
println("  median = $(sort(rel_diff)[length(rel_diff)÷2])")
println("  max    = $(maximum(rel_diff))")
println()
for f in (20.0, 50.0, 100.0, 150.0, 200.0, 300.0, 500.0, 1000.0)
    idx = argmin(abs.(gwdata_freqs .- f))
    println("  f=$(gwdata_freqs[idx]) Hz  GWData=$(gwdata_psd[idx])  gwpy=$(gwpy_psd[idx])  rel_diff=$(abs(gwdata_psd[idx]-gwpy_psd[idx])/gwpy_psd[idx])")
end
