#!/usr/bin/env julia
# GravitationalWaveData.jl side of compare_psd.jl (see README.md): computes `psd` over the given
# off-source window and writes it to a CSV.
#
# Usage: julia --project=<GravitationalWaveData.jl root> compare_psd_gwdata_step.jl <h1_path> <off_start> <off_end> <duration> <out_csv>

using GravitationalWaveData

h1_path, off_start_s, off_end_s, duration_s, out_csv = ARGS
off_source_start = parse(Float64, off_start_s)
off_source_end = parse(Float64, off_end_s)
duration = parse(Float64, duration_s)

strain, _dq_mask = load_gwosc_strain(h1_path)
off_source = crop(strain, off_source_start, off_source_end)
nseg = round(Int, duration / strain.dt)

result = psd(off_source; seglen=nseg, overlap=0.5, method=:median)

open(out_csv, "w") do io
    println(io, "freq,psd")
    for (f, p) in zip(frequencies(result), result.data)
        println(io, "$f,$p")
    end
end
