# GravitationalWaveData.jl Bonito notebook: Q-transform explorer
#
# Bonito.jl counterpart of `../pluto/qtransform_explorer.jl`. Interactively explore GravitationalWaveData.jl's
# constant-Q transform (`q_plane_transform`) on the GW150914 H1 segment. Adjust **Q**, the
# frequency search band (**low**/**high**), and whether the segment is whitened first, and watch
# the time-frequency energy map (the "Q-scan"/"chirp plot") update live.
#
# This mirrors `docs/literate/06_qtransform.jl`'s static tutorial, but lets you interactively pick
# a single fixed Q (via `q_plane_transform` directly) instead of the automatic multi-Q search
# `q_transform` performs.
#
# Run standalone and open the printed URL in a browser:
#   julia --project=notebooks/bonito notebooks/bonito/qtransform_explorer.jl
#
# Uses its own environment (this directory) rather than `docs`, since Bonito's HTTP.jl
# requirement (2.x) conflicts with Pluto's (1.x) and the two can't share one environment.

import Pkg
Pkg.activate(@__DIR__)

using GravitationalWaveData, Plots, Bonito, Base64
import Bonito.TailwindDashboard as D

gr()

# --- Fetching real strain data ---
# The same standard GW150914 H1 segment the rest of this project's tutorials use, fetched live
# from GWOSC via `fetch_open_data` (cached locally under `~/.cache/gravitationalwavedata/strain`, so
# re-running this notebook doesn't re-download the file). Cropped to 16 s centered on the
# merger, same as `docs/literate/06_qtransform.jl`.
const gps_start, gps_end = 1126256640, 1126260736
const gps_merger = 1126259462.4   # GW150914 merger GPS time

const qt_path = fetch_open_data("H1", gps_start, gps_end)
const qt_strain, qt_dq_mask = load_gwosc_strain(qt_path)
const qt_seg = crop(qt_strain, gps_merger - 8.0, gps_merger + 8.0)

# A whitened version of the same window is also precomputed once here (not reactive -- whitening
# a ~24 s stretch is the expensive part, so it's done a single time up front rather than on every
# slider move). Padding beyond the 16 s window we care about keeps the whitening FIR filter's
# edge transient (fduration/2 = 1 s each side by default) away from the merger.
const qt_pad = 8.0 + 4.0
const qt_wide = crop(qt_strain, gps_merger - qt_pad, gps_merger + qt_pad)
const qt_asd_est = asd(qt_strain; seglen=round(Int, 4 * sample_rate(qt_strain)))
const qt_white_wide = whiten(qt_wide, qt_asd_est; highpass=20.0)
const qt_white_seg = crop(qt_white_wide, gps_merger - 8.0, gps_merger + 8.0)

# `QScan.rows` holds one energy row per tile frequency, each with its own time resolution (see
# `QScan`'s docstring) -- deliberately not spline-regridded onto a common grid by the package
# itself (that's presentation-only). This short linear-interpolation helper, lifted from
# `docs/literate/06_qtransform.jl`, resamples every row onto one shared time grid for the
# heatmap.
function _lininterp(xs, ys, x)
    x <= xs[1] && return ys[1]
    x >= xs[end] && return ys[end]
    i = searchsortedlast(xs, x)
    x0, x1 = xs[i], xs[i + 1]
    y0, y1 = ys[i], ys[i + 1]
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
end

function qscan_grid(scan, tlo, thi; ntime=400)
    common_times = range(tlo, thi; length=ntime)
    freqs = [row.frequency for row in scan.rows]
    z = Array{Float64}(undef, length(scan.rows), ntime)
    for (i, row) in enumerate(scan.rows)
        n = length(row.energy)
        row_times = row.t0 .+ row.dt .* (0:(n - 1))
        z[i, :] .= [_lininterp(row_times, row.energy, t) for t in common_times]
    end
    return collect(common_times), freqs, z
end

# Render a Plots.jl figure to a data-URI <img>, since Bonito doesn't natively embed Plots.jl
# figures the way it does Makie scenes.
function _img(p)
    io = IOBuffer()
    show(io, MIME"image/png"(), p)
    return DOM.img(src="data:image/png;base64,$(base64encode(take!(io)))",
                    style="max-width:100%;")
end

# Everything here depends on the controls below and re-runs automatically whenever they change.
function qscan_and_plot(q, flo, fhi, whiten_on)
    flo, fhi = min(flo, fhi), max(flo, fhi)
    input = whiten_on ? qt_white_seg : qt_seg
    fseries = normalized_fft(input)
    scan = q_plane_transform(fseries.data, q, (flo, fhi), duration(input), sample_rate(input),
                              input.t0)
    peak = qscan_peak(scan; search=(gps_merger - 0.2, gps_merger + 0.1))

    grid_times, grid_freqs, energy_grid = qscan_grid(scan, gps_merger - 0.15, gps_merger + 0.08)
    rel_times = grid_times .- gps_merger

    p = heatmap(rel_times, grid_freqs, energy_grid; xlabel="Time from merger (s)",
                ylabel="Frequency (Hz)", yscale=:log10, colorbar_title="Normalized energy",
                color=:viridis, title="H1 Q-transform (Q=$(round(q; digits=2)), whiten=$whiten_on)")

    peak_info = DOM.p("Peak tile: t - merger = $(round(peak.time - gps_merger; digits=4)) s, ",
                       "f = $(round(peak.frequency; digits=1)) Hz, ",
                       "SNR = $(round(peak.snr; digits=2))")

    return DOM.div(_img(p), peak_info)
end

app = App() do session
    # - Q: the constant-Q value of the single Q-plane computed (gwpy/LIGO detchar convention;
    #   higher Q means better frequency resolution, worse time resolution).
    # - low/high: the frequency band (Hz) the plane is tiled over.
    # - whiten: whether the segment is whitened (against a full-segment ASD estimate) before the
    #   Q-transform, matching the "whiten before Q-scan" step common in detector
    #   characterization.
    q = D.Slider("Q: ", 4.0:0.5:64.0; value=8.0)
    flo = D.Slider("Frequency low (Hz): ", 10.0:5.0:995.0; value=20.0)
    fhi = D.Slider("Frequency high (Hz): ", 15.0:5.0:1000.0; value=500.0)
    whiten_cb = D.Checkbox("Whiten", false)

    plots = map(qscan_and_plot, q.value, flo.value, fhi.value, whiten_cb.value)

    controls = D.Card(D.FlexCol(D.Title("Q-transform controls"), q, flo, fhi, whiten_cb))
    return D.FlexRow(controls, D.Card(plots))
end

server = Bonito.Server(app, "127.0.0.1", 8082)
println("Q-transform explorer running at http://127.0.0.1:8082")
wait(server)
