module Coherence

export csd, coherence

using FFTW
using ..Types: TimeSeries, FrequencySeries
using ..PSD: psd, _welch_starts, _welch_window

# ================================================================================================
# Cross-spectral density (CSD) and (magnitude-squared) coherence between two channels, following
# `scipy.signal.csd`/`scipy.signal.coherence`:
#   - Always mean-averaged (`average="mean"`). A median of *complex* cross-spectral values isn't a
#     well-defined robust statistic the way a median of real, non-negative PSD values is.
#   - Segmentation and windowing are `PSD.jl`'s own (`_welch_starts`/`_welch_window`: Hann window,
#     50% overlap by default), reused directly rather than re-derived.
#   - Normalization (from `scipy.signal._spectral_py._spectral_helper`): cross-periodogram
#     `Pxy = conj(rfft(x)) * rfft(y) * scale` with `scale = dt / (seglen * mean(win.^2))`,
#     one-sided bins above DC/Nyquist doubled (DC and, for even `seglen`, Nyquist are left
#     alone) -- exactly `psd`'s normalization generalized from `abs2(F)` to a cross term.
#   - No detrending. scipy subtracts the per-segment mean by default (`detrend="constant"`), but
#     `psd` does not, and detrending only the cross term would make `coherence`'s numerator and
#     denominator inconsistent. The effect on real (near-zero-mean) strain data is negligible; see
#     `examples/validation/compare_coherence.jl` for the measured agreement.
# ================================================================================================

function _check_csd_compatible(a::TimeSeries, b::TimeSeries)
    a.dt == b.dt ||
        throw(ArgumentError("TimeSeries sampling rates differ: dt=$(a.dt)s ($(1/a.dt) Hz) " *
                             "vs dt=$(b.dt)s ($(1/b.dt) Hz)"))
    length(a) == length(b) ||
        throw(ArgumentError("TimeSeries lengths differ: $(length(a)) vs $(length(b))"))
    return nothing
end

"""
    csd(ts1::TimeSeries, ts2::TimeSeries; seglen::Int, overlap::Real=0.5) -> FrequencySeries

Welch-method one-sided cross-spectral density estimate between two `TimeSeries` sampled on the same
grid (same `dt`; `t0` need not match -- only the bin-by-bin sample correspondence used for each FFT
segment matters, matching `scipy.signal.csd`'s own lack of an epoch concept). See this module's
header comment for the normalization/averaging conventions.

Segments as in `psd` (same `_welch_starts`/`_welch_window` helpers: `seglen`-sample segments,
fractional `overlap`, Hann window), mean-averaged (not median). Complex-valued: `Pxy(f) =
mean_over_segments(conj(rfft(x)) .* rfft(y)) * scale`, one-sided density normalization
(`scipy.signal.welch`'s `scaling="density"`).

Returns a `FrequencySeries{ComplexF64}` with `f0=0.0`, `df=1/(seglen*dt)`.
"""
function csd(ts1::TimeSeries{T1}, ts2::TimeSeries{T2}; seglen::Int,
             overlap::Real=0.5) where {T1<:Real,T2<:Real}
    _check_csd_compatible(ts1, ts2)
    x = ts1.data
    y = ts2.data
    dt = ts1.dt
    n = length(x)
    starts = _welch_starts(n, seglen, overlap)
    win, win_power = _welch_window(seglen)

    freqs = FFTW.rfftfreq(seglen, 1 / dt)
    nf = length(freqs)
    cross = Matrix{ComplexF64}(undef, nf, length(starts))
    for (k, s) in enumerate(starts)
        xseg = @view(x[s:(s + seglen - 1)]) .* win
        yseg = @view(y[s:(s + seglen - 1)]) .* win
        Fx = FFTW.rfft(xseg)
        Fy = FFTW.rfft(yseg)
        @. cross[:, k] = (2 * dt / seglen) * conj(Fx) * Fy / win_power
    end
    # Same DC/Nyquist bin-halving rationale as `psd` (see PSD.jl): the blanket "times 2" above
    # double-counts the two bins that have no negative-frequency mirror to fold in.
    cross[1, :] ./= 2
    iseven(seglen) && (cross[end, :] ./= 2)

    csd_data = vec(sum(cross; dims=2)) ./ length(starts)

    name = (ts1.name === nothing && ts2.name === nothing) ? nothing :
           "$(something(ts1.name, "?"))---$(something(ts2.name, "?")) CSD"
    df = 1 / (seglen * dt)
    return FrequencySeries(csd_data, 0.0, df; name=name)
end

"""
    coherence(ts1::TimeSeries, ts2::TimeSeries; seglen::Int, overlap::Real=0.5) -> FrequencySeries

Magnitude-squared (Welch-averaged) coherence between two `TimeSeries`: `|Pxy(f)|^2 / (Pxx(f) *
Pyy(f))`, real-valued, in `[0, 1]`. Follows `scipy.signal.coherence`; see this
module's header comment for the normalization/averaging conventions.

Requires averaging over multiple segments to be meaningful: a single-segment estimate is
trivially 1.0 at every frequency (`Pxy`, `Pxx`, `Pyy` all computed from the same one FFT), which is
exactly why `csd`/`coherence` always mean-average over `seglen`-sample,
`overlap`-fraction-overlapping segments rather than accepting `seglen == length(ts)`.

`Pxx`/`Pyy` are computed via `psd(ts; seglen, overlap, method=:mean)` (mean-averaged, matching
`scipy.signal.coherence`'s own `average="mean"` default and keeping `Pxy`'s numerator and `Pxx`/
`Pyy`'s denominator on the exact same averaging convention); `Pxy` via `csd`. Requires matching
`dt`/length, same as `csd`.

Returns a `FrequencySeries{Float64}` with `f0=0.0`, `df=1/(seglen*dt)`.
"""
function coherence(ts1::TimeSeries{T1}, ts2::TimeSeries{T2}; seglen::Int,
                    overlap::Real=0.5) where {T1<:Real,T2<:Real}
    _check_csd_compatible(ts1, ts2)
    pxy = csd(ts1, ts2; seglen=seglen, overlap=overlap)
    pxx = psd(ts1; seglen=seglen, overlap=overlap, method=:mean)
    pyy = psd(ts2; seglen=seglen, overlap=overlap, method=:mean)

    coh_data = abs2.(pxy.data) ./ (pxx.data .* pyy.data)

    name = (ts1.name === nothing && ts2.name === nothing) ? nothing :
           "Coherence between $(something(ts1.name, "?")) and $(something(ts2.name, "?"))"
    return FrequencySeries(coh_data, pxy.f0, pxy.df; name=name)
end

end
