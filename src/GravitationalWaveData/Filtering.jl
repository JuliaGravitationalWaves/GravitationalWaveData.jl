module Filtering

export bandpass, highpass, lowpass, notch, whiten

using DSP
using FFTW: irfft
using ..Types: TimeSeries, FrequencySeries, sample_rate, duration, frequencies

# ================================================================================================
# Bandpass/highpass/lowpass/notch: fixed-order Butterworth filter design (`DSP.Filters.Butterworth`)
# applied zero-phase via `DSP.filtfilt`. The filter order is explicit rather than auto-selected
# from passband-ripple/stopband-attenuation specs (as `scipy.signal.iirdesign`-based tools do), so
# results are predictable and comparable across sample rates. `examples/validation/
# compare_filtering.jl` therefore checks functional behavior (passband transparency, stopband
# attenuation) against auto-order designs rather than sample-for-sample agreement.
# ================================================================================================

"""
    bandpass(ts::TimeSeries, flow::Real, fhigh::Real; order::Integer=4) -> TimeSeries

Zero-phase Butterworth band-pass filter, passing frequencies in `(flow, fhigh)` Hz. Filter design
via `DSP.digitalfilter(DSP.Bandpass(flow, fhigh), DSP.Butterworth(order); fs=sample_rate(ts))`,
applied with `DSP.filtfilt` (forward-backward, so the *effective* order is `2*order` and there is
no phase distortion).
"""
function bandpass(ts::TimeSeries, flow::Real, fhigh::Real; order::Integer=4)
    filt = digitalfilter(Bandpass(flow, fhigh), Butterworth(order); fs=sample_rate(ts))
    return TimeSeries(filtfilt(filt, ts.data), ts.t0, ts.dt; name=ts.name)
end

"""
    highpass(ts::TimeSeries, fcut::Real; order::Integer=4) -> TimeSeries

Zero-phase Butterworth high-pass filter with corner frequency `fcut` Hz. See [`bandpass`](@ref)
for the design/application convention.
"""
function highpass(ts::TimeSeries, fcut::Real; order::Integer=4)
    filt = digitalfilter(Highpass(fcut), Butterworth(order); fs=sample_rate(ts))
    return TimeSeries(filtfilt(filt, ts.data), ts.t0, ts.dt; name=ts.name)
end

"""
    lowpass(ts::TimeSeries, fcut::Real; order::Integer=4) -> TimeSeries

Zero-phase Butterworth low-pass filter with corner frequency `fcut` Hz. See [`bandpass`](@ref) for
the design/application convention.
"""
function lowpass(ts::TimeSeries, fcut::Real; order::Integer=4)
    filt = digitalfilter(Lowpass(fcut), Butterworth(order); fs=sample_rate(ts))
    return TimeSeries(filtfilt(filt, ts.data), ts.t0, ts.dt; name=ts.name)
end

"""
    notch(ts::TimeSeries, f0::Real; bandwidth::Real=1.0, order::Integer=4) -> TimeSeries

Zero-phase Butterworth band-*stop* (notch) filter, rejecting a narrow band `(f0 - bandwidth/2,
f0 + bandwidth/2)` Hz -- e.g. mains-hum removal at `f0=60` (US) or `f0=50` (EU). Same
Butterworth+filtfilt convention as [`bandpass`](@ref), using `DSP.Bandstop` + `DSP.Butterworth`.
"""
function notch(ts::TimeSeries, f0::Real; bandwidth::Real=1.0, order::Integer=4)
    filt = digitalfilter(Bandstop(f0 - bandwidth / 2, f0 + bandwidth / 2), Butterworth(order);
                          fs=sample_rate(ts))
    return TimeSeries(filtfilt(filt, ts.data), ts.t0, ts.dt; name=ts.name)
end

# ================================================================================================
# Whitening by inverse spectrum truncation: build a *time-domain FIR filter* approximating
# `1/ASD(f)` (truncated to `fduration` seconds and tapered so it doesn't have unbounded support),
# then convolve that FIR filter with the (detrended, edge-tapered) data using the overlap-save
# method -- never a bare frequency-domain divide, which rings badly at segment edges. Background
# theory: arXiv:gr-qc/0509116.
#
# The helpers below follow gwpy's implementation function-for-function
# (`gwpy/signal/filter_design.py::fir_from_transfer`/`_truncate_transfer`/`_truncate_impulse`,
# `gwpy/timeseries/timeseries.py::TimeSeries.whiten`/`.convolve`,
# `FrequencySeries.interpolate`), so window conventions and normalization match it exactly.
# ================================================================================================

"""
    _periodic_hann(n) -> Vector{Float64}

The periodic/DFT-even Hann window `scipy.signal.get_window("hann", n)` returns by default
(`fftbins=true`) -- distinct from `DSP.Windows.hanning(n)`, which defaults to the *symmetric*
convention (`fftbins=false`). Computed directly as `0.5*(1-cos(2*pi*k/n))` for `k in 0:n-1`, which
is exactly a length-`(n+1)` symmetric Hann window with its last (duplicate-of-first) sample
dropped -- the standard periodic-window trick. gwpy relies on this specific convention (via
`scipy.signal.get_window`'s default) in `_truncate_impulse` and `TimeSeries.convolve`'s edge
tapers, so it's reproduced exactly here rather than approximated with `DSP.Windows.hanning`.
"""
_periodic_hann(n::Integer) = [0.5 * (1 - cos(2pi * k / n)) for k in 0:(n - 1)]

function _get_window(window::Symbol, n::Integer)
    window === :hann && return _periodic_hann(n)
    throw(ArgumentError("whiten: unsupported window $(repr(window)) (only :hann is implemented)"))
end

"""
    _truncate_transfer(transfer, ncorner) -> Vector{Float64}

Direct port of gwpy's `filter_design._truncate_transfer`: zero the first `ncorner` bins (hard
highpass) and smoothly taper 5 samples at each edge of the remaining array with a symmetric
(`sym=true`-equivalent) Tukey(alpha=1) == Hann(10) window, so the transfer function doesn't have a
sharp discontinuity at its own array boundaries before being inverse-FFT'd.
"""
function _truncate_transfer(transfer::Vector{Float64}, ncorner::Integer)
    nsamp = length(transfer)
    out = copy(transfer)
    ncorner > 0 && (out[1:ncorner] .= 0.0)
    ntaper = nsamp - ncorner
    if ntaper > 0
        nsamples = 5
        w = tukey(2 * nsamples, 1.0)   # == symmetric Hann(10); DSP.Windows default (zerophase=false)
                                        # matches scipy.signal.windows.tukey's default sym=true, which
                                        # is what gwpy calls here (not the periodic get_window default).
        out[(ncorner + 1):(ncorner + nsamples)] .*= w[1:nsamples]
        out[(nsamp - nsamples + 1):nsamp] .*= w[(end - nsamples + 1):end]
    end
    return out
end

"""
    _truncate_impulse(impulse, ntaps, window) -> Vector{Float64}

Direct port of gwpy's `filter_design._truncate_impulse`: window and zero the impulse response so
only `ntaps` total non-tapered-to-zero samples remain (split between the two ends, wrapped around
what will become the FIR's center once rolled in [`_fir_from_transfer`](@ref)).
"""
function _truncate_impulse(impulse::Vector{Float64}, ntaps::Integer, window::Symbol)
    out = copy(impulse)
    trunc_start = ntaps ÷ 2
    trunc_stop = length(out) - trunc_start
    w = _get_window(window, ntaps)
    out[1:trunc_start] .*= w[(trunc_start + 1):ntaps]
    out[(trunc_stop + 1):end] .*= w[1:trunc_start]
    out[(trunc_start + 1):trunc_stop] .= 0.0
    return out
end

"""
    _fir_from_transfer(transfer, ntaps; window=:hann, ncorner=0) -> Vector{Float64}

Direct port of gwpy's `filter_design.fir_from_transfer`: design a length-`ntaps` (should be even)
FIR filter realizing the given one-sided (`rfft`-style, DC-first) frequency-domain `transfer`
function (here always `1 ./ asd`, see [`whiten`](@ref)), via truncate -> inverse-FFT -> truncate ->
circularly-shift-to-center.
"""
function _fir_from_transfer(transfer::Vector{Float64}, ntaps::Integer; window::Symbol=:hann,
                             ncorner::Integer=0)
    t = _truncate_transfer(transfer, ncorner)
    impulse = irfft(Complex{Float64}.(t), 2 * (length(t) - 1))
    impulse = _truncate_impulse(impulse, ntaps, window)
    return circshift(impulse, ntaps ÷ 2 - 1)[1:ntaps]
end

"""
    _numpy_interp(xnew, x, y) -> Vector{Float64}

1-D linear interpolation matching `numpy.interp`'s convention: `x` must be sorted ascending;
`xnew` values outside `[x[1], x[end]]` are clamped to the corresponding boundary `y` value (no
extrapolation), matching what `FrequencySeries.interpolate` relies on.
"""
function _numpy_interp(xnew::AbstractVector{<:Real}, x::AbstractVector{<:Real},
                        y::AbstractVector{<:Real})
    n = length(x)
    out = Vector{Float64}(undef, length(xnew))
    for (i, xv) in enumerate(xnew)
        if xv <= x[1]
            out[i] = y[1]
        elseif xv >= x[n]
            out[i] = y[n]
        else
            j = clamp(searchsortedlast(x, xv), 1, n - 1)
            x0, x1 = x[j], x[j + 1]
            y0, y1 = y[j], y[j + 1]
            out[i] = x1 == x0 ? y0 : y0 + (y1 - y0) * (xv - x0) / (x1 - x0)
        end
    end
    return out
end

"""
    _interpolate(fs::FrequencySeries, new_df::Real) -> FrequencySeries

Direct port of gwpy's `FrequencySeries.interpolate`: linear interpolation (`numpy.interp`) onto a
new frequency grid with spacing `new_df`, same `f0`, length `round((N-1)*(df/new_df)) + 1`.
"""
function _interpolate(fs::FrequencySeries, new_df::Real)
    n_samples = round(Int, (length(fs) - 1) * (fs.df / new_df)) + 1
    new_freqs = fs.f0 .+ new_df .* (0:(n_samples - 1))
    new_vals = _numpy_interp(new_freqs, collect(frequencies(fs)), fs.data)
    return FrequencySeries(new_vals, fs.f0, Float64(new_df); name=fs.name)
end

"""
    _convolve_same(a, b) -> Vector{Float64}

Full linear convolution of `a` and `b`, cropped to `scipy.signal.fftconvolve(a, b,
mode="same")`'s convention: the `length(a)`-sample slice of the full (`length(a)+length(b)-1`)
convolution, centered (`start = (full_len - length(a)) ÷ 2`, floor division, 0-based in scipy's own
formula).
"""
function _convolve_same(a::Vector{Float64}, b::Vector{Float64})
    full = conv(a, b)
    n = length(a)
    total = length(full)
    start = (total - n) ÷ 2
    return full[(start + 1):(start + n)]
end

"""
    _overlap_save_convolve(data, fir, window) -> Vector{Float64}

Direct port of gwpy's `TimeSeries.convolve` (overlap-save FFT convolution with edge tapering to
suppress spectral leakage from the filter's settle-in transient at the two boundaries). For the
segment lengths this package's whitening is normally used on (a handful of seconds to a few
hundred seconds, `fir` a couple of seconds long), the `nfft >= n/2` branch (single `fftconvolve`
call) is what actually executes; the chunked overlap-save branch is exercised by
`test/runtests.jl`'s dedicated long-signal test (forcing `n` well past `8*length(fir)`) checked
against the single-shot path, since real whitening usage on this package's typical segment lengths
never reaches it otherwise.
"""
function _overlap_save_convolve(data::Vector{Float64}, fir::Vector{Float64}, window::Symbol)
    n = length(data)
    m = length(fir)
    pad = cld(m, 2)   # ceil(m/2)
    nfft = min(8 * m, n)
    w = _get_window(window, m)
    in_ = copy(data)
    in_[1:pad] .*= w[1:pad]
    in_[(end - pad + 1):end] .*= w[(end - pad + 1):end]

    if nfft >= n / 2
        return _convolve_same(in_, fir)
    end

    nstep = nfft - 2 * pad
    conv_out = zeros(Float64, n)
    # first chunk
    first_yk = _convolve_same(in_[1:nfft], fir)
    conv_out[1:(nfft - pad)] .= first_yk[1:(nfft - pad)]
    # middle chunks
    k = nfft - pad
    while k < n - nfft + pad
        yk = _convolve_same(in_[(k - pad + 1):(k + nstep + pad)], fir)
        ylen = length(yk)
        conv_out[(k + 1):(k + ylen - 2 * pad)] .= yk[(pad + 1):(ylen - pad)]
        k += nstep
    end
    # last chunk
    last_yk = _convolve_same(in_[(n - nfft + 1):n], fir)
    conv_out[(n - nfft + pad + 1):n] .= last_yk[(pad + 1):nfft]
    return conv_out
end

"""
    whiten(ts::TimeSeries, psd_or_asd::FrequencySeries; kind=:asd, fduration=2.0,
           highpass=nothing, window=:hann) -> TimeSeries

Whiten `ts` using inverse spectrum truncation: design a length-`fduration`
(seconds) time-domain FIR filter approximating `1/ASD(f)` (via [`_fir_from_transfer`](@ref)) and
convolve it with the detrended data using the overlap-save method (via
[`_overlap_save_convolve`](@ref)) -- *not* a naive FFT/divide/inverse-FFT, which would ring badly
at the segment edges.

`psd_or_asd` supplies the spectrum to whiten against, typically from
[`asd`](@ref GravitationalWaveData.PSD.asd) or [`psd`](@ref GravitationalWaveData.PSD.psd) on an off-source stretch; pass
`kind=:psd` if you're passing a power (not amplitude) spectral density. It is linearly interpolated onto the frequency grid with spacing
`1/duration(ts)` before use, so it need not already be sampled at that resolution (e.g. it's fine to
reuse the same off-source-estimated ASD for whitening several different on-source segments).

`highpass` (Hz), if given, hard-zeros the FIR's frequency response below that corner in addition
to whatever `1/ASD` already does there.

On stationary Gaussian input, the output has zero mean and unit variance. A segment of length
`fduration/2` at each end of the output is corrupted by filter settle-in (tapered instead of
exact) -- crop it off if that matters for downstream use.
"""
function whiten(ts::TimeSeries, psd_or_asd::FrequencySeries; kind::Symbol=:asd,
                 fduration::Real=2.0, highpass::Union{Real,Nothing}=nothing,
                 window::Symbol=:hann)
    asd = if kind === :asd
        psd_or_asd
    elseif kind === :psd
        FrequencySeries(sqrt.(psd_or_asd.data), psd_or_asd.f0, psd_or_asd.df; name=psd_or_asd.name)
    else
        throw(ArgumentError("whiten: kind must be :asd or :psd, got $(repr(kind))"))
    end

    target_df = 1.0 / duration(ts)
    asd_i = _interpolate(asd, target_df)

    ncorner = highpass === nothing ? 0 : trunc(Int, highpass / asd_i.df)
    ntaps = trunc(Int, fduration * sample_rate(ts))
    ntaps > 0 || throw(ArgumentError("whiten: fduration*sample_rate(ts) must be positive, got $ntaps"))

    tdw = _fir_from_transfer(1.0 ./ Vector{Float64}(asd_i.data), ntaps; window=window,
                              ncorner=ncorner)

    x = Vector{Float64}(ts.data)
    detrended = x .- (sum(x) / length(x))   # detrend="constant": subtract the mean
    conv_out = _overlap_save_convolve(detrended, tdw, window)
    out = conv_out .* sqrt(2 * ts.dt)
    return TimeSeries(out, ts.t0, ts.dt; name=ts.name)
end

end
