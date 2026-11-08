#!/usr/bin/env julia
# GravitationalWaveData.jl side of compare_resample.jl (see README.md): crops the H1 segment,
# resamples it with `resample`, and writes a plain (mean, not median) Welch PSD of the result,
# computed with `DSP.welch_pgram` to match gwpy's `method="welch"`, to a CSV. The PSD is only a
# lens on resampling fidelity, not a test of `psd` itself.
#
# Usage: julia --project=<GravitationalWaveData.jl root> compare_resample_gwdata_step.jl <h1_path> <seg_start> <seg_end> <new_rate> <out_csv>

using GravitationalWaveData
import DSP   # qualified: DSP also exports `resample`, which would otherwise clash with GravitationalWaveData's

h1_path, start_s, end_s, new_rate_s, out_csv = ARGS
seg_start = parse(Float64, start_s)
seg_end = parse(Float64, end_s)
new_rate = parse(Float64, new_rate_s)

strain, _ = load_gwosc_strain(h1_path)
seg = crop(strain, seg_start, seg_end)
resampled = GravitationalWaveData.resample(seg, new_rate)

n = round(Int, 4.0 * new_rate)     # 4s FFT segments, same convention as compare_psd.jl
noverlap = n ÷ 2                    # 50% overlap
pgram = DSP.welch_pgram(resampled.data, n, noverlap; fs=new_rate, window=DSP.hanning)

open(out_csv, "w") do io
    println(io, "freq,psd")
    for (f, p) in zip(DSP.freq(pgram), DSP.power(pgram))
        println(io, "$f,$p")
    end
end
