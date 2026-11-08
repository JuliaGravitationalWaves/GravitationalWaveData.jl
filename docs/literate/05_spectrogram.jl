# # Spectrogram

# A spectrogram shows how a signal's frequency content evolves over time: a sequence of
# short-time Fourier transforms (STFTs), one per time bin, stacked side by side into a
# time-frequency image. [`spectrogram`](@ref) (`src/GravitationalWaveData/Spectrogram.jl`) is a thin wrapper
# around `DSP.spectrogram`: a raw, un-averaged STFT power spectrogram (no Welch/median averaging of
# sub-segments within each column). This tutorial uses a short, heavily-overlapped FFT window to
# resolve a short-duration transient, then shows why a transient this short can be invisible at the
# wrong resolution.

using GravitationalWaveData
using Plots
gr()

# ## Fetching and conditioning real strain data
#
# We reuse the same standard 4096 s H1 segment covering GW150914 that the other tutorials use
# (cached locally after the first fetch), and the same whitening approach as the Filtering &
# Whitening tutorial: estimate an ASD from the whole 4096 s segment (long enough that GW150914's
# own sub-second signal can't meaningfully bias it) and whiten against it, hard-zeroing the
# response below 20 Hz via `whiten`'s own `highpass` kwarg. Whitening flattens the spectrum so the
# STFT below isn't dominated by low-frequency seismic noise -- without it, a spectrogram of raw
# strain shows essentially nothing but that noise floor.

gps_start, gps_end = 1126256640, 1126260736
gps_merger = 1126259462.4   # GW150914 merger GPS time

path = fetch_open_data("H1", gps_start, gps_end)
strain, dq_mask = load_gwosc_strain(path)

asd_est = asd(strain; seglen=round(Int, 4 * sample_rate(strain)))
white = whiten(strain, asd_est; highpass=20.0);

# [`crop`](@ref) down to a 32 s window centered on the merger -- "a few tens of seconds", per this
# tutorial's brief, rather than the full 4096 s segment (unnecessary and slow to render as a
# spectrogram). This window sits deep in the interior of the 4096 s segment we whitened, well clear
# of the `fduration/2`-second edge transient `whiten` leaves at the very start/end of its output
# (see `docs/CONVENTIONS.md`'s "Whitening" section).

seg = crop(white, gps_merger - 16, gps_merger + 16)

# ## An over-dense STFT spectrogram
#
# We use a very short, heavily-overlapped FFT window -- 1/16 s segments overlapping by 15/16 -- to
# over-sample the frequency content and give a short transient a chance to stand out. In samples at
# H1's 4096 Hz sample rate, that's an `n=256`-sample segment with `noverlap=240` (a ~94% overlap
# between consecutive segments, well above `spectrogram`'s own 50%-overlap default):

n = 256          # 1/16 s at 4096 Hz
noverlap = 240   # 15/16 = ~94% overlap

spec = spectrogram(seg; n=n, noverlap=noverlap)
size(spec)   # (nfreq bins, ntime bins)

# `spectrogram` returns power; we take its square root to get a more familiar amplitude spectral
# density scale before plotting, building a new [`Spectrogram`](@ref) from the square-rooted data
# (the axis metadata -- `t0`, `dt`, `f0`, `df` -- is unchanged by this, only the values):

spec_asd = Spectrogram(sqrt.(spec.data), spec.t0, spec.dt, spec.f0, spec.df;
                        name="Whitened strain ASD")

# ## A first, full-window look
#
# The [`Spectrogram`](@ref) plot recipe defaults to a heatmap with a log-scaled colorbar --
# appropriate here since STFT power spans many orders of magnitude across frequency. We restrict
# the frequency axis to `20-1024 Hz` (below 20 Hz is already hard-zeroed by whitening's own
# highpass; above ~1 kHz there's nothing of astrophysical interest for GW150914):

plot(spec_asd; title="H1, whitened, over-dense spectrogram (32 s around GW150914)",
     ylim=(20, 1024))
#-

# GW150914's chirp lasts well under a second, out of the 32 s plotted above -- so before drawing
# any conclusion from squinting at the image, it's worth checking numerically whether there's
# actually anything to see. There is: the single loudest pixel in this entire 32 s, `size(spec, 2)`
# = 8177-time-bin spectrogram sits almost exactly at the merger GPS time:

global_max, global_idx = findmax(spec_asd.data)
peak_freq = frequencies(spec_asd)[global_idx[1]]
peak_time_from_merger = times(spec_asd)[global_idx[2]] - gps_merger
global_max, peak_freq, peak_time_from_merger

# That confirms the transient really is there, and roughly where physics says it should be
# (GW150914's merger frequency is close to that ~64 Hz bin). But at the full 32 s width plotted
# above, that one bright column is a needle in a haystack of ~8000 columns -- easy to miss entirely
# by eye: at this resolution and dynamic range, a signal this short is practically invisible. The
# Q-transform tutorial (`06_qtransform`) is specifically the better tool for finding a sub-second
# transient like this -- it adapts its time/frequency resolution per frequency rather than using one
# fixed STFT window throughout. Here, we can still get a fair look at it with a plain spectrogram by
# zooming in.

# ## Zooming in around the merger
#
# GravitationalWaveData.jl doesn't yet provide a dedicated `crop` for `Spectrogram` the way it does for
# `TimeSeries`/`FrequencySeries` (`docs/CONVENTIONS.md`), so we restrict both axes by hand, the
# same half-open-window idea applied to a matrix instead of a vector:

function crop_spectrogram(s::Spectrogram, t_range, f_range)
    tvals = collect(times(s))
    fvals = collect(frequencies(s))
    tmask = t_range[1] .<= tvals .<= t_range[2]
    fmask = f_range[1] .<= fvals .<= f_range[2]
    return Spectrogram(s.data[fmask, tmask], first(tvals[tmask]), s.dt, first(fvals[fmask]), s.df;
                        name=s.name)
end

zoomed = crop_spectrogram(spec_asd, (gps_merger - 0.15, gps_merger + 0.05), (20.0, 400.0))
size(zoomed)

# Restricting the data itself (not just the plotted axis limits) also tightens the color scale to
# the range that actually matters here, rather than one set by the 32 s window's overall extremes.
# The [`Spectrogram`](@ref) recipe always plots absolute GPS time; for a ~0.2 s window we'd rather
# read the axis relative to the merger, so we call `heatmap` directly instead of going through the
# recipe (the same `(freq, time)` row/column order `Spectrogram` stores, and the same
# `colorbar_scale`/`colorbar_title` attributes the recipe sets, just spelled out by hand):

rel_times = collect(times(zoomed)) .- gps_merger
heatmap(rel_times, collect(frequencies(zoomed)), zoomed.data;
        xlabel="Time from merger (s)", ylabel="Frequency (Hz)",
        colorbar_scale=:log10, colorbar_title="Whitened strain ASD",
        title="H1 spectrogram around GW150914 (zoomed)", color=:viridis)
#-

# ## What this actually shows
#
# Zoomed in like this, a rising streak is visible by eye: power climbs from around 32 Hz roughly
# 0.1 s before the merger, up through the 60-150 Hz band right around the merger GPS time, then
# fades within another ~0.03 s -- the classic inspiral-merger-ringdown shape of a binary black hole
# coalescence, and broadly consistent with GW150914's actual reported peak frequency (around 150
# Hz). We can quantify the contrast against the background rather than just eyeballing it: comparing
# the same global peak found above to the whole spectrogram's typical (median) pixel value,

background_level = sort(vec(spec_asd.data))[length(spec_asd.data) ÷ 2]
global_max, background_level, global_max / background_level

# the peak is roughly an order of magnitude above the typical background pixel -- a real, if modest,
# excursion, not an artifact of the color scale. That said, this is a much less dramatic picture
# than the sharp, high-contrast chirp track a matched-filter-style time-frequency method produces: a
# fixed-window STFT (even a heavily-overlapped one like this) trades time resolution against
# frequency resolution the same way at every frequency, whereas GW150914's chirp sweeps through
# several octaves in a fraction of a second -- exactly the regime where a fixed-resolution STFT
# struggles and a multi-resolution method doesn't. See the Q-transform tutorial (`06_qtransform`)
# for that sharper, purpose-built view of the same GW150914 chirp.
