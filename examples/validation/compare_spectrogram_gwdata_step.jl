#!/usr/bin/env julia
# GravitationalWaveData.jl side of compare_spectrogram.jl (see README.md): computes `spectrogram` over
# the given segment and writes its frequency axis, GPS time axis and power matrix to CSVs.
#
# Usage: julia --project=<GravitationalWaveData.jl root> compare_spectrogram_gwdata_step.jl <h1_path> <seg_start> <seg_end> <n> <noverlap> <out_prefix>

using GravitationalWaveData

h1_path, seg_start_s, seg_end_s, n_s, noverlap_s, out_prefix = ARGS
seg_start = parse(Float64, seg_start_s)
seg_end = parse(Float64, seg_end_s)
n = parse(Int, n_s)
noverlap = parse(Int, noverlap_s)

strain, _ = load_gwosc_strain(h1_path)
seg = crop(strain, seg_start, seg_end)

spec = spectrogram(seg; n=n, noverlap=noverlap)

open(out_prefix * "_freqs.csv", "w") do io
    for f in frequencies(spec)
        println(io, f)
    end
end

open(out_prefix * "_times.csv", "w") do io
    for tt in times(spec)
        println(io, tt)
    end
end

open(out_prefix * "_power.csv", "w") do io
    nf, nt = size(spec)
    for i in 1:nf
        println(io, join(spec.data[i, :], ","))
    end
end
