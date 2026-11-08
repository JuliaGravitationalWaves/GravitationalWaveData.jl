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

# ╔═╡ 4495fe80-aad6-11f1-8949-4f7f54d8dc4f
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, "..", "..", "docs"))
end


# ╔═╡ 4495fe80-aad6-11f1-84fb-617462881561
using GravitationalWaveData, Plots, PlutoUI, Markdown, HypertextLiteral

# ╔═╡ 4495fe80-aad6-11f1-9443-f1fb60600b92
md"""
# Filter explorer

Interactively design a zero-phase Butterworth filter ([`bandpass`](@ref)/[`highpass`](@ref)/
[`lowpass`](@ref)/[`notch`](@ref)) on the GW150914 H1 segment and watch the before/after ASD
update live. Mirrors `docs/literate/03_filtering_whitening.jl`'s static filter-comparison plots.
"""


# ╔═╡ 4495fe80-aad6-11f1-b632-ff9ac4d64d8e
gr()

# ╔═╡ 4495fe80-aad6-11f1-b302-a74bf24002d5
md"""
## Fetching real strain data

Same standard GW150914 H1 segment as the rest of this project's tutorials, fetched via
[`fetch_open_data`](@ref) (cached locally under `~/.cache/gravitationalwavedata/strain`). Cropped to 16 s
around the merger, same as `docs/literate/03_filtering_whitening.jl`.
"""


# ╔═╡ 4495fe80-aad6-11f1-a285-57af32948ffe
begin
    gps_start, gps_end = 1126256640, 1126260736
    gps_merger = 1126259462.4   # GW150914 merger GPS time

    fx_path = fetch_open_data("H1", gps_start, gps_end)
    fx_strain, fx_dq_mask = load_gwosc_strain(fx_path)
    fx_seg = crop(fx_strain, gps_merger - 8.0, gps_merger + 8.0)
    fx_seglen = round(Int, 4 * sample_rate(fx_seg))   # 4 s Welch segments for the ASD estimates
    (duration(fx_seg), sample_rate(fx_seg))
end


# ╔═╡ 4495fe80-aad6-11f1-9932-777651f79356
begin
    fx_seg_asd = crop(asd(fx_seg; seglen=fx_seglen), 8.0, nothing)   # drop f=0 bin (log axes)
end;


# ╔═╡ 4495fe80-aad6-11f1-b9b5-efd34dc683f2
md"""
## Controls

- **filter type**: bandpass / highpass / lowpass / notch.
- **low / cutoff frequency**: used as the low corner (bandpass), the cutoff (highpass), or the
  notch center frequency `f0` (notch). Ignored for lowpass.
- **high frequency**: used as the high corner (bandpass) or the cutoff (lowpass). Ignored for
  highpass/notch.
- **notch bandwidth**: only used when filter type is notch.
- **order**: Butterworth filter order (applied zero-phase via `filtfilt`, so the *effective*
  order is `2*order`).
"""


# ╔═╡ 4495fe80-aad6-11f1-88fe-ef7df79e74c5
@bind fx_controls PlutoUI.combine() do Child
    @htl("""
    <div style="position:fixed; top:80px; right:20px; z-index:500; background:white;
                color:#222; border:1px solid #999; border-radius:10px; padding:12px 16px;
                box-shadow:0 2px 12px rgba(0,0,0,0.3); width:230px; font-size:14px;">
      <div class="fx-handle" style="cursor:move; font-weight:bold; user-select:none;
                                     padding-bottom:6px; margin-bottom:6px; border-bottom:1px solid #ddd;">
        ⠿ Filter controls
      </div>
      <div style="margin-top:8px;">Type<br>$(Child("filter_type", Select(["bandpass", "highpass", "lowpass", "notch"]; default="bandpass")))</div>
      <div style="margin-top:8px;">Low / cutoff (Hz)<br>$(Child("flow", Slider(1.0:1.0:1000.0; default=50.0, show_value=true)))</div>
      <div style="margin-top:8px;">High (Hz)<br>$(Child("fhigh", Slider(1.0:1.0:1000.0; default=250.0, show_value=true)))</div>
      <div style="margin-top:8px;">Notch bandwidth<br>$(Child("notch_bw", Slider(0.5:0.5:20.0; default=2.0, show_value=true)))</div>
      <div style="margin-top:8px;">Order<br>$(Child("order", Slider(1:1:8; default=4, show_value=true)))</div>
    </div>
    <script>
    const panel = currentScript.previousElementSibling;
    const handle = panel.querySelector('.fx-handle');
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

# ╔═╡ 4495fe80-aad6-11f1-b90b-0fc20b9a2137
fx_filter_type, fx_flow, fx_fhigh, fx_notch_bw, fx_order = fx_controls.filter_type,
    fx_controls.flow, fx_controls.fhigh, fx_controls.notch_bw, fx_controls.order;

# ╔═╡ 4495fe80-aad6-11f1-9c12-89cfd86be48b
md"""
## Live filter design + before/after ASD

Everything below depends on the controls above and re-runs automatically whenever they change.
"""


# ╔═╡ 4495fe80-aad6-11f1-b108-617d9496a7c1
begin
    fx_filtered = if fx_filter_type == "bandpass"
        fx_lo, fx_hi = min(fx_flow, fx_fhigh), max(fx_flow, fx_fhigh)
        bandpass(fx_seg, fx_lo, fx_hi; order=fx_order)
    elseif fx_filter_type == "highpass"
        highpass(fx_seg, fx_flow; order=fx_order)
    elseif fx_filter_type == "lowpass"
        lowpass(fx_seg, fx_fhigh; order=fx_order)
    elseif fx_filter_type == "notch"
        notch(fx_seg, fx_flow; bandwidth=fx_notch_bw, order=fx_order)
    else
        error("unknown filter type: $fx_filter_type")
    end
    fx_filtered
end;


# ╔═╡ 4495fe80-aad6-11f1-b272-0198a1c2eddf
fx_filt_asd = crop(asd(fx_filtered; seglen=fx_seglen), 8.0, nothing);


# ╔═╡ 4495fe80-aad6-11f1-a3ba-d30c313b190b
begin
    fx_p = plot(fx_seg_asd; label="raw", xlim=(8, 2000))
    plot!(fx_p, fx_filt_asd; label="$(fx_filter_type)")
    fx_p
end


# ╔═╡ 4495fe80-aad6-11f1-8981-5d908e8ce33a
md"""
Below: a short time-domain zoom around the merger, raw vs. filtered, so you can see the effect
in the time domain too (not just the ASD).
"""


# ╔═╡ 4495fe80-aad6-11f1-98d3-17491096c80a
begin
    fx_zoom_lo, fx_zoom_hi = gps_merger - 0.3, gps_merger + 0.3
    fx_raw_zoom = crop(fx_seg, fx_zoom_lo, fx_zoom_hi)
    fx_filt_zoom = crop(fx_filtered, fx_zoom_lo, fx_zoom_hi)
    fx_zoom_rel = collect(times(fx_raw_zoom)) .- gps_merger
    fx_tp = plot(fx_zoom_rel, fx_raw_zoom.data; xlabel="Time from merger (s)",
                 ylabel="Amplitude (strain)", label="raw")
    plot!(fx_tp, fx_zoom_rel, fx_filt_zoom.data; label="$(fx_filter_type)")
    fx_tp
end


# ╔═╡ Cell order:
# ╠═4495fe80-aad6-11f1-9443-f1fb60600b92
# ╠═4495fe80-aad6-11f1-8949-4f7f54d8dc4f
# ╠═4495fe80-aad6-11f1-84fb-617462881561
# ╠═4495fe80-aad6-11f1-b632-ff9ac4d64d8e
# ╠═4495fe80-aad6-11f1-b9b5-efd34dc683f2
# ╠═4495fe80-aad6-11f1-88fe-ef7df79e74c5
# ╠═4495fe80-aad6-11f1-b90b-0fc20b9a2137
# ╠═4495fe80-aad6-11f1-9c12-89cfd86be48b
# ╠═4495fe80-aad6-11f1-a3ba-d30c313b190b
# ╠═4495fe80-aad6-11f1-8981-5d908e8ce33a
# ╠═4495fe80-aad6-11f1-98d3-17491096c80a
# ╠═4495fe80-aad6-11f1-b302-a74bf24002d5
# ╠═4495fe80-aad6-11f1-a285-57af32948ffe
# ╠═4495fe80-aad6-11f1-9932-777651f79356
# ╠═4495fe80-aad6-11f1-b108-617d9496a7c1
# ╠═4495fe80-aad6-11f1-b272-0198a1c2eddf
