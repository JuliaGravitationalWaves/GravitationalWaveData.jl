#!/usr/bin/env julia
# Compare GravitationalWaveData.jl's `QTransform.q_transform` against gwpy's `gwpy.signal.qtransform.q_scan`
# (the function gwpy's own `TimeSeries.q_transform()` calls internally) on the same real GW150914
# H1 strain segment.
#
# Both sides scan the same un-whitened data, so the gwpy side calls the `q_scan` primitive
# directly rather than the `TimeSeries.q_transform()` wrapper, which whitens by default. Both
# sides also read the RAW (pre-interpolation) peak -- GravitationalWaveData.jl's `q_transform` stops there by
# design (see the QTransform.jl module docstring on why the spline regridding is presentation-only).
#
# Usage (see README.md in this directory):
#   julia examples/validation/compare_qtransform.jl

include("common.jl")

const H1_PATH = gw150914_strain_path("H1")

# 4s segment straddling the GW150914 trigger, unwindowed (a plain crop -- NOT
# `condition_segment`'s Tukey-tapered analysis segment; gwpy's `.fft()`/`q_scan` on a cropped
# `TimeSeries` doesn't apply any window either, so this keeps both sides doing the same thing).
# `frange=(20, 500)` restricts the scan to the band GW150914's chirp actually sweeps through,
# matching standard Omega-scan practice of excluding the sub-20 Hz seismic wall and calibration
# lines above the signal band -- both sides use the same restriction.
const TRIGGER_TIME = 1126259462.4
const SEG_START = TRIGGER_TIME - 2.0
const SEG_END = TRIGGER_TIME + 2.0
const FMIN = 20.0
const FMAX = 500.0

gwdata_csv = output_path("qtransform_gwdata.csv")
gwpy_csv = output_path("qtransform_gwpy.csv")

run_gwdata_step("compare_qtransform_gwdata_step.jl", H1_PATH, SEG_START, SEG_END, FMIN, FMAX, gwdata_csv)
run_reference_step("compare_qtransform_gwpy_step.jl", H1_PATH, SEG_START, SEG_END, FMIN, FMAX, gwpy_csv)

function read_result_csv(path)
    out = Dict{String,Float64}()
    for line in readlines(path)[2:end]
        k, v = split(line, ',')
        out[k] = parse(Float64, v)
    end
    return out
end

gwdata = read_result_csv(gwdata_csv)
gwpy = read_result_csv(gwpy_csv)

println("Comparing QTransform.q_transform (GravitationalWaveData.jl) vs gwpy.signal.qtransform.q_scan, same H1 segment")
println("  segment: GPS [$SEG_START, $SEG_END) (4s around GW150914 trigger $TRIGGER_TIME)")
println("  frange = ($FMIN, $FMAX) Hz, qrange = (4, 64), mismatch = 0.2, norm = median, no whitening")
println()
for field in ("time", "frequency", "energy", "snr", "q")
    gd = gwdata[field]
    gp = gwpy[field]
    reldiff = gp == 0 ? abs(gd - gp) : abs(gd - gp) / abs(gp)
    println("  $field: GravitationalWaveData=$gd  gwpy=$gp  " * (field == "time" ? "abs_diff=$(abs(gd-gp))s" : "rel_diff=$reldiff"))
end
