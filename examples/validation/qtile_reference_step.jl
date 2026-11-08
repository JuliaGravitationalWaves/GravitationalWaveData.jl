#!/usr/bin/env julia
# Reference values for a single Q-transform tile, used by test/qtransform.jl. Computes a QTile's
# window, data indices, ntiles, windowsize, padding and raw energy with gwpy itself, for the same
# deterministic signal and tile parameters the test checks `QTransform.jl` against. The exact
# indexing/normalization convention (in particular a `+1` offset) is pinned down only by gwpy's
# source, not its docs or Chatterji's thesis, hence this direct check.
#
# Not part of run_all.jl: run it by hand in this directory's project (see README.md) when the
# test's reference values need regenerating:
#   julia --project=examples/validation examples/validation/qtile_reference_step.jl \
#       <signal_csv> <q> <freq> <dur> <fs> <mismatch> <out_csv>

using PythonCall

signal_csv, q_s, freq_s, dur_s, fs_s, mismatch_s, out_csv = ARGS
q = parse(Float64, q_s)
freq = parse(Float64, freq_s)
dur = parse(Float64, dur_s)
fs = parse(Float64, fs_s)
mismatch = parse(Float64, mismatch_s)

np = pyimport("numpy")
gwpy_ts_mod = pyimport("gwpy.timeseries")
qmod = pyimport("gwpy.signal.qtransform")

data = np.loadtxt(signal_csv)
ts = gwpy_ts_mod.TimeSeries(data, sample_rate=fs, epoch=0.0)
fseries = ts.fft()

tile = qmod.QTile(q, freq, dur, fs, mismatch=mismatch)
window = pyconvert(Vector{Float64}, tile.get_window())
indices0 = pyconvert(Vector{Int}, tile.get_data_indices())
padding = pyconvert(Tuple{Int,Int}, tile.padding)
ntiles = pyconvert(Int, tile.ntiles)
windowsize = pyconvert(Int, tile.windowsize)
qprime = pyconvert(Float64, tile.qprime)
energy_ts = tile.transform(fseries, norm=false, epoch=0.0)
energy = pyconvert(Vector{Float64}, energy_ts.value)

open(out_csv, "w") do io
    println(io, "ntiles,$ntiles")
    println(io, "windowsize,$windowsize")
    println(io, "qprime,$qprime")
    println(io, "padding_left,$(padding[1])")
    println(io, "padding_right,$(padding[2])")
    println(io, "data_indices_0based," * join(indices0, ";"))
    println(io, "window," * join(window, ";"))
    println(io, "energy," * join(energy, ";"))
end

println("ntiles=$ntiles windowsize=$windowsize qprime=$qprime padding=$padding")
println("data_indices[1:5] (0-based)=", indices0[1:5])
println("window[1:5]=", window[1:5])
maxidx = argmax(energy)
println("energy max=$(energy[maxidx]) at 0-based index $(maxidx-1)")
