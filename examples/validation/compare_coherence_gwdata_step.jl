#!/usr/bin/env julia
# GWData.jl side of compare_coherence.jl (see README.md): computes `csd`/`coherence`
# between the given H1/L1 off-source windows and writes both to a CSV.
#
# Usage: julia --project=<GWData.jl root> compare_coherence_gwdata_step.jl <h1_path> <l1_path> <off_start> <off_end> <duration> <out_csv>

using GWData

h1_path, l1_path, off_start_s, off_end_s, duration_s, out_csv = ARGS
off_source_start = parse(Float64, off_start_s)
off_source_end = parse(Float64, off_end_s)
duration = parse(Float64, duration_s)

h1_strain, _h1_dq = load_gwosc_strain(h1_path)
l1_strain, _l1_dq = load_gwosc_strain(l1_path)

h1_off = crop(h1_strain, off_source_start, off_source_end)
l1_off = crop(l1_strain, off_source_start, off_source_end)

nseg = round(Int, duration / h1_strain.dt)

csd_result = csd(h1_off, l1_off; seglen=nseg, overlap=0.5)
coh_result = coherence(h1_off, l1_off; seglen=nseg, overlap=0.5)

open(out_csv, "w") do io
    println(io, "freq,csd_re,csd_im,coherence")
    for (f, c, k) in zip(frequencies(csd_result), csd_result.data, coh_result.data)
        println(io, "$f,$(real(c)),$(imag(c)),$k")
    end
end
