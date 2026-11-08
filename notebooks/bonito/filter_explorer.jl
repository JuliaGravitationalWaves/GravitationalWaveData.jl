# GWData.jl Bonito notebook: Filter explorer
#
# Bonito.jl counterpart of `../pluto/filter_explorer.jl`. Interactively design a zero-phase
# Butterworth filter (`bandpass`/`highpass`/`lowpass`/`notch`) on the GW150914 H1 segment and
# watch the before/after ASD (and a time-domain zoom) update live. Mirrors
# `docs/literate/03_filtering_whitening.jl`'s static filter-comparison plots.
#
# Run standalone and open the printed URL in a browser:
#   julia --project=notebooks/bonito notebooks/bonito/filter_explorer.jl
#
# Uses its own environment (this directory) rather than `docs`, since Bonito's HTTP.jl
# requirement (2.x) conflicts with Pluto's (1.x) and the two can't share one environment.

import Pkg
Pkg.activate(@__DIR__)

using GWData, Plots, Bonito, Base64
import Bonito.TailwindDashboard as D

gr()

# --- Fetching real strain data ---
# Same standard GW150914 H1 segment as the rest of this project's tutorials, fetched via
# `fetch_open_data` (cached locally under `~/.cache/gwdata/strain`). Cropped to 16 s around
# the merger, same as `docs/literate/03_filtering_whitening.jl`.
const gps_start, gps_end = 1126256640, 1126260736
const gps_merger = 1126259462.4   # GW150914 merger GPS time

const fx_path = fetch_open_data("H1", gps_start, gps_end)
const fx_strain, fx_dq_mask = load_gwosc_strain(fx_path)
const fx_seg = crop(fx_strain, gps_merger - 8.0, gps_merger + 8.0)
const fx_seglen = round(Int, 4 * sample_rate(fx_seg))   # 4 s Welch segments for the ASD estimates
const fx_seg_asd = crop(asd(fx_seg; seglen=fx_seglen), 8.0, nothing)   # drop f=0 bin (log axes)

# Render a Plots.jl figure to a data-URI <img>, since Bonito doesn't natively embed Plots.jl
# figures the way it does Makie scenes.
function _img(p)
    io = IOBuffer()
    show(io, MIME"image/png"(), p)
    return DOM.img(src="data:image/png;base64,$(base64encode(take!(io)))",
                    style="max-width:100%;")
end

# Everything here depends on the controls below and re-runs automatically whenever they change.
function filter_and_plot(filter_type, flow, fhigh, notch_bw, order)
    filtered = if filter_type == "bandpass"
        lo, hi = min(flow, fhigh), max(flow, fhigh)
        bandpass(fx_seg, lo, hi; order=order)
    elseif filter_type == "highpass"
        highpass(fx_seg, flow; order=order)
    elseif filter_type == "lowpass"
        lowpass(fx_seg, fhigh; order=order)
    elseif filter_type == "notch"
        notch(fx_seg, flow; bandwidth=notch_bw, order=order)
    else
        error("unknown filter type: $filter_type")
    end

    filt_asd = crop(asd(filtered; seglen=fx_seglen), 8.0, nothing)
    p_asd = plot(fx_seg_asd; label="raw", xlim=(8, 2000))
    plot!(p_asd, filt_asd; label=filter_type)

    # Short time-domain zoom around the merger, raw vs. filtered, so the effect is visible in
    # the time domain too (not just the ASD).
    zoom_lo, zoom_hi = gps_merger - 0.3, gps_merger + 0.3
    raw_zoom = crop(fx_seg, zoom_lo, zoom_hi)
    filt_zoom = crop(filtered, zoom_lo, zoom_hi)
    zoom_rel = collect(times(raw_zoom)) .- gps_merger
    p_time = plot(zoom_rel, raw_zoom.data; xlabel="Time from merger (s)",
                   ylabel="Amplitude (strain)", label="raw")
    plot!(p_time, zoom_rel, filt_zoom.data; label=filter_type)

    return DOM.div(_img(p_asd), _img(p_time))
end

app = App() do session
    # - filter type: bandpass / highpass / lowpass / notch.
    # - low / cutoff frequency: the low corner (bandpass), the cutoff (highpass), or the notch
    #   center frequency f0 (notch). Ignored for lowpass.
    # - high frequency: the high corner (bandpass) or the cutoff (lowpass). Ignored for
    #   highpass/notch.
    # - notch bandwidth: only used when filter type is notch.
    # - order: Butterworth filter order (applied zero-phase via filtfilt, so the *effective*
    #   order is 2*order).
    filter_type = D.Dropdown("Type: ", ["bandpass", "highpass", "lowpass", "notch"])
    flow = D.Slider("Low / cutoff (Hz): ", 1.0:1.0:1000.0; value=50.0)
    fhigh = D.Slider("High (Hz): ", 1.0:1.0:1000.0; value=250.0)
    notch_bw = D.Slider("Notch bandwidth: ", 0.5:0.5:20.0; value=2.0)
    order = D.Slider("Order: ", 1:1:8; value=4)

    plots = map(filter_and_plot, filter_type.value, flow.value, fhigh.value, notch_bw.value,
                order.value)

    controls = D.Card(D.FlexCol(D.Title("Filter controls"), filter_type, flow, fhigh, notch_bw,
                                 order))
    return D.FlexRow(controls, D.Card(plots))
end

server = Bonito.Server(app, "127.0.0.1", 8081)
println("Filter explorer running at http://127.0.0.1:8081")
wait(server)
