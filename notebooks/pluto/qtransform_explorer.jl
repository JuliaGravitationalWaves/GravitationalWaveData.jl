### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 41062102-aad6-11f1-bc7b-5149286d63a4
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", "docs"))
end


# ╔═╡ 41062102-aad6-11f1-ab39-85978b114f29
using GWData, Plots, PlutoUI, Markdown, HypertextLiteral

# ╔═╡ 4105abd0-aad6-11f1-b0c9-99fcc22264bb
md"""
# Q-transform explorer

Interactively explore GWData.jl's constant-Q transform ([`q_plane_transform`](@ref)) on the
GW150914 H1 segment. Adjust **Q**, the frequency search band (**frange**), and whether the
segment is whitened first, and watch the time-frequency energy map (the "Q-scan"/"chirp plot")
update live.

This mirrors `docs/literate/06_qtransform.jl`'s static tutorial, but lets you interactively pick
a single fixed Q (via [`q_plane_transform`](@ref) directly) instead of the automatic multi-Q
search [`q_transform`](@ref) performs.
"""


# ╔═╡ 41062102-aad6-11f1-a0b6-21df7edcf934
gr()

# ╔═╡ 41062102-aad6-11f1-a19c-fd85daf6b6c5
md"""
## Fetching real strain data

The same standard GW150914 H1 segment the rest of this project's tutorials use, fetched live
from GWOSC via [`fetch_open_data`](@ref) (cached locally under `~/.cache/gwdata/strain`, so
re-running this notebook doesn't re-download the file). We crop to 16 s centered on the merger,
same as `docs/literate/06_qtransform.jl`.
"""


# ╔═╡ 41062102-aad6-11f1-b521-b901dc77d6cb
begin
    gps_start, gps_end = 1126256640, 1126260736
    gps_merger = 1126259462.4   # GW150914 merger GPS time

    qt_path = fetch_open_data("H1", gps_start, gps_end)
    qt_strain, qt_dq_mask = load_gwosc_strain(qt_path)
    qt_seg = crop(qt_strain, gps_merger - 8.0, gps_merger + 8.0)
    (duration(qt_seg), sample_rate(qt_seg))
end


# ╔═╡ 41062102-aad6-11f1-bb94-db22b75a0560
md"""
A whitened version of the same window is also precomputed once here (not reactive -- whitening
a ~24 s stretch is the expensive part, so it's done a single time up front rather than on every
slider move). Padding beyond the 16 s window we care about keeps the whitening FIR filter's
edge transient (`fduration/2` = 1 s each side by default) away from the merger.
"""


# ╔═╡ 41062102-aad6-11f1-9d3f-5d6f241b9c17
begin
    qt_pad = 8.0 + 4.0
    qt_wide = crop(qt_strain, gps_merger - qt_pad, gps_merger + qt_pad)
    qt_asd_est = asd(qt_strain; seglen=round(Int, 4 * sample_rate(qt_strain)))
    qt_white_wide = whiten(qt_wide, qt_asd_est; highpass=20.0)
    qt_white_seg = crop(qt_white_wide, gps_merger - 8.0, gps_merger + 8.0)
end;


# ╔═╡ 41064810-aad6-11f1-b109-b53670354214
md"""
## Controls

- **Q**: the constant-Q value of the single Q-plane computed (LIGO detchar convention;
  higher Q means better frequency resolution, worse time resolution).
- **frange**: the (low, high) Hz band the plane is tiled over.
- **whiten**: whether the segment is whitened (against a full-segment ASD estimate) before the
  Q-transform, matching the "whiten before Q-scan" step common in detector characterization.
"""


# ╔═╡ 41064810-aad6-11f1-89f3-151418985e08
@bind qt_controls PlutoUI.combine() do Child
    @htl("""
    <div style="position:fixed; top:80px; right:20px; z-index:500; background:white;
                color:#222; border:1px solid #999; border-radius:10px; padding:12px 16px;
                box-shadow:0 2px 12px rgba(0,0,0,0.3); width:230px; font-size:14px;">
      <div class="qt-handle" style="cursor:move; font-weight:bold; user-select:none;
                                     padding-bottom:6px; margin-bottom:6px; border-bottom:1px solid #ddd;">
        ⠿ Q-transform controls
      </div>
      <div style="margin-top:8px;">Q<br>$(Child("q", Slider(4.0:0.5:64.0; default=8.0, show_value=true)))</div>
      <div style="margin-top:8px;">Frequency range (Hz)<br>$(Child("frange", RangeSlider(10:5:1000; default=20:5:500)))</div>
      <div style="margin-top:8px;">Whiten $(Child("whiten", CheckBox(default=false)))</div>
    </div>
    <script>
    const panel = currentScript.previousElementSibling;
    const handle = panel.querySelector('.qt-handle');
    let dragging = false, startX = 0, startY = 0, origX = 0, origY = 0;
    handle.addEventListener('mousedown', (e) => {
        dragging = true;
        startX = e.clientX; startY = e.clientY;
        const r = panel.getBoundingClientRect();
        origX = r.left; origY = r.top;
        e.preventDefault();
    });
    window.addEventListener('mousemove', (e) => {
        if (!dragging) return;
        panel.style.left = (origX + e.clientX - startX) + 'px';
        panel.style.top = (origY + e.clientY - startY) + 'px';
        panel.style.right = 'auto';
    });
    window.addEventListener('mouseup', () => { dragging = false; });
    </script>
    """)
end

# ╔═╡ 41064810-aad6-11f1-bf04-99f7b3430e2c
qt_q, qt_frange, qt_whiten = qt_controls.q, qt_controls.frange, qt_controls.whiten;

# ╔═╡ 41064810-aad6-11f1-bf8f-df28573ab677
md"""
## Regridding helper

`QScan.rows` holds one energy row per tile frequency, each with its own time resolution (see
[`QScan`](@ref)'s docstring) -- deliberately not spline-regridded onto a common grid by the
package itself (that's presentation-only). This short linear-interpolation helper, lifted from
`docs/literate/06_qtransform.jl`, resamples every row onto one shared time grid for the heatmap.
"""


# ╔═╡ 41064810-aad6-11f1-b2d4-e106276f0f63
function _lininterp(xs, ys, x)
    x <= xs[1] && return ys[1]
    x >= xs[end] && return ys[end]
    i = searchsortedlast(xs, x)
    x0, x1 = xs[i], xs[i + 1]
    y0, y1 = ys[i], ys[i + 1]
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
end


# ╔═╡ 41064810-aad6-11f1-921e-839647347ded
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


# ╔═╡ 41064810-aad6-11f1-86ca-95000ef5d488
md"""
## Live Q-scan

Everything below depends on the sliders above and re-runs automatically whenever they change.
"""


# ╔═╡ 41064810-aad6-11f1-b55a-6f7ac41a295c
begin
    qt_flo, qt_fhi = Float64(first(qt_frange)), Float64(last(qt_frange))
    qt_input = qt_whiten ? qt_white_seg : qt_seg
    qt_fseries = normalized_fft(qt_input)
    qt_scan = q_plane_transform(qt_fseries.data, qt_q, (qt_flo, qt_fhi),
                                 duration(qt_input), sample_rate(qt_input), qt_input.t0)
    qt_peak = qscan_peak(qt_scan; search=(gps_merger - 0.2, gps_merger + 0.1))
    qt_scan, qt_peak
end;


# ╔═╡ 41064810-aad6-11f1-8728-5ba3f0b40193
begin
    qt_grid_times, qt_grid_freqs, qt_energy_grid = qscan_grid(qt_scan, gps_merger - 0.15,
                                                               gps_merger + 0.08)
    qt_rel_times = qt_grid_times .- gps_merger
end;


# ╔═╡ 41064810-aad6-11f1-bd1c-031fcde54d04
heatmap(qt_rel_times, qt_grid_freqs, qt_energy_grid;
        xlabel="Time from merger (s)", ylabel="Frequency (Hz)", yscale=:log10,
        colorbar_title="Normalized energy", color=:viridis,
        title="H1 Q-transform (Q=$(round(qt_q; digits=2)), whiten=$(qt_whiten))")


# ╔═╡ 41064810-aad6-11f1-9502-05f20ad9473a
md"""
**Peak tile:** t - merger = $(round(qt_peak.time - gps_merger; digits=4)) s,
f = $(round(qt_peak.frequency; digits=1)) Hz, SNR = $(round(qt_peak.snr; digits=2))
"""


# ╔═╡ Cell order:
# ╠═4105abd0-aad6-11f1-b0c9-99fcc22264bb
# ╠═41062102-aad6-11f1-bc7b-5149286d63a4
# ╠═41062102-aad6-11f1-ab39-85978b114f29
# ╠═41062102-aad6-11f1-a0b6-21df7edcf934
# ╠═41064810-aad6-11f1-b109-b53670354214
# ╠═41064810-aad6-11f1-89f3-151418985e08
# ╠═41064810-aad6-11f1-bf04-99f7b3430e2c
# ╠═41064810-aad6-11f1-86ca-95000ef5d488
# ╠═41064810-aad6-11f1-bd1c-031fcde54d04
# ╠═41064810-aad6-11f1-9502-05f20ad9473a
# ╠═41062102-aad6-11f1-a19c-fd85daf6b6c5
# ╠═41062102-aad6-11f1-b521-b901dc77d6cb
# ╠═41062102-aad6-11f1-bb94-db22b75a0560
# ╠═41062102-aad6-11f1-9d3f-5d6f241b9c17
# ╠═41064810-aad6-11f1-bf8f-df28573ab677
# ╠═41064810-aad6-11f1-b2d4-e106276f0f63
# ╠═41064810-aad6-11f1-921e-839647347ded
# ╠═41064810-aad6-11f1-b55a-6f7ac41a295c
# ╠═41064810-aad6-11f1-8728-5ba3f0b40193
