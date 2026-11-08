#!/usr/bin/env julia
# GWData.jl side of compare_qtransform.jl (see README.md): computes the raw
# `q_transform` peak over the given H1 segment and writes it to a CSV.
#
# Usage: julia --project=<GWData.jl root> compare_qtransform_gwdata_step.jl <h1_path> <seg_start> <seg_end> <fmin> <fmax> <out_csv>

using GWData

h1_path, seg_start_s, seg_end_s, fmin_s, fmax_s, out_csv = ARGS
seg_start = parse(Float64, seg_start_s)
seg_end = parse(Float64, seg_end_s)
fmin = parse(Float64, fmin_s)
fmax = parse(Float64, fmax_s)

strain, _dq = load_gwosc_strain(h1_path)
segment = crop(strain, seg_start, seg_end)

_scan, peak = q_transform(segment; frange=(fmin, fmax))

open(out_csv, "w") do io
    println(io, "field,value")
    println(io, "time,$(peak.time)")
    println(io, "frequency,$(peak.frequency)")
    println(io, "energy,$(peak.energy)")
    println(io, "snr,$(peak.snr)")
    println(io, "q,$(peak.q)")
end
