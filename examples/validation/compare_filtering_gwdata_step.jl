#!/usr/bin/env julia
# GWData.jl side of compare_filtering.jl (see README.md): applies bandpass/highpass/
# lowpass to the analysis segment and whitens against the ASD the reference step wrote, writing
# each result to a CSV in <out_dir>.
#
# Usage: julia --project=<GWData.jl root> compare_filtering_gwdata_step.jl <h1_path> <out_dir>

using GWData

h1_path, out_dir = ARGS

strain, _dq = load_gwosc_strain(h1_path)

const TRIGGER_TIME = 1126259462.4
const FILTER_SEG_START = TRIGGER_TIME - 16.0
const FILTER_SEG_END = TRIGGER_TIME + 16.0
const WHITEN_SEG_START = TRIGGER_TIME - 6.0
const WHITEN_SEG_END = TRIGGER_TIME + 2.0

function write_csv(path, xs, ys)
    open(path, "w") do io
        println(io, "x,y")
        for (x, y) in zip(xs, ys)
            println(io, "$x,$y")
        end
    end
end

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

# --- bandpass / highpass / lowpass ---------------------------------------------------------------
filter_seg = crop(strain, FILTER_SEG_START, FILTER_SEG_END)
write_csv(joinpath(out_dir, "raw_segment_julia.csv"), times(filter_seg), filter_seg.data)

bp = bandpass(filter_seg, 35.0, 350.0; order=4)
hp = highpass(filter_seg, 35.0; order=4)
lp = lowpass(filter_seg, 350.0; order=4)

write_csv(joinpath(out_dir, "julia_bandpass.csv"), times(bp), bp.data)
write_csv(joinpath(out_dir, "julia_highpass.csv"), times(hp), hp.data)
write_csv(joinpath(out_dir, "julia_lowpass.csv"), times(lp), lp.data)

# --- whiten, against the gwpy-computed off-source ASD written to asd.csv by the gwpy step --------
asd_freqs, asd_vals = read_csv(joinpath(out_dir, "asd.csv"))
asd_df = asd_freqs[2] - asd_freqs[1]
asd = FrequencySeries(asd_vals, asd_freqs[1], asd_df)

whiten_seg = crop(strain, WHITEN_SEG_START, WHITEN_SEG_END)
whitened = whiten(whiten_seg, asd; fduration=2.0)
write_csv(joinpath(out_dir, "julia_whiten.csv"), times(whitened), whitened.data)

println("GWData step done: wrote julia_{bandpass,highpass,lowpass,whiten}.csv to $out_dir")
