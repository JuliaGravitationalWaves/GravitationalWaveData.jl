module PSD

export psd, asd

using FFTW
using DSP: DSP
using Statistics: median
using ..Types: TimeSeries, FrequencySeries

# ================================================================================================
# Median-Welch PSD/ASD estimation: Hann window, DC/Nyquist bin-halving, `ln(2)` median debiasing,
# standard density normalization -- see `psd`'s docstring below and `docs/CONVENTIONS.md`.
#
# Window: `DSP.Windows.hanning(N)` (default `padding=0, zerophase=false`), the "symmetric"
# (non-periodic) Hann convention `0.5 * (1 - cos(2*pi*n/(N-1)))`, endpoints exactly zero.
# ================================================================================================

"""
    _welch_starts(n::Integer, seglen::Integer, overlap::Real) -> AbstractRange{Int}

1-based starting indices of Welch segments of length `seglen` with fractional `overlap`, covering
a series of length `n`. Factored out of `psd` so `Coherence.jl`'s `csd`/`coherence` use the exact
same segmentation convention (same step size, same "shorter than one segment" error) rather than a
hand-mirrored copy that could silently drift from this one.
"""
function _welch_starts(n::Integer, seglen::Integer, overlap::Real)
    step = max(1, round(Int, seglen * (1 - overlap)))
    starts = 1:step:(n - seglen + 1)
    isempty(starts) &&
        throw(ArgumentError("series (length $n) is shorter than one segment (length $seglen)"))
    return starts
end

"""
    _welch_window(seglen::Integer) -> (win, win_power)

The Hann window (`DSP.Windows.hanning`) used for every Welch-segmented estimate in this package
(`psd`/`asd` here, `csd`/`coherence` in `Coherence.jl`), plus its power-loss correction factor
`mean(win.^2)`.
"""
function _welch_window(seglen::Integer)
    win = DSP.Windows.hanning(seglen)
    win_power = sum(abs2, win) / seglen   # mean(w.^2), corrects for the window's power loss
    return win, win_power
end

"""
    psd(ts::TimeSeries; seglen::Int, overlap::Real=0.5, method::Symbol=:median) -> FrequencySeries

Welch-method one-sided power spectral density estimate over a (typically long, off-source)
`TimeSeries`: split into `seglen`-sample segments with fractional `overlap`, Hann-window each
(`DSP.Windows.hanning`, also `scipy.signal.welch`'s default). This is deliberately a different
window from `condition_segment`'s Tukey taper on the actual analysis segment -- that one exists to
suppress edge discontinuities in a single non-periodic segment being FFT'd once, not to run a Welch
average.

`method=:median` (default, the standard LIGO convention) uses the per-frequency-bin median across
segments rather than the mean, which is far more robust to loud transient glitches contaminating a
handful of segments; debiased by the exponential distribution's `ln(2)` median correction
(asymptotically exact for the segment counts here -- O(500) segments for a 1024s off-source
stretch at 4s/50% overlap -- see e.g. Allen et al. 2005 for the finite-sample correction this
approximates). `method=:mean` skips both the median and the bias correction.

Normalization is the standard density convention (`scipy.signal.welch`'s `scaling="density"`):
`∫psd(f) df` over `[0, Nyquist]` approximates the time series' variance.

Segments containing any non-finite sample are left out of the average, with a warning, since GWOSC
files mark gaps with `NaN` and one such segment would otherwise make every bin `NaN`. An
`ArgumentError` is thrown if no segment is finite.

Returns a `FrequencySeries` with `f0=0.0`, `df=1/(seglen*dt)`, real-valued data.
"""
function psd(ts::TimeSeries{T}; seglen::Int, overlap::Real=0.5,
             method::Symbol=:median) where {T<:Real}
    data = ts.data
    dt = ts.dt
    n = length(data)
    all_starts = _welch_starts(n, seglen, overlap)
    starts = [s for s in all_starts if all(isfinite, @view(data[s:(s + seglen - 1)]))]
    isempty(starts) && throw(ArgumentError("psd: all $(length(all_starts)) segments contain non-finite samples"))
    if length(starts) < length(all_starts)
        @warn "psd: left out $(length(all_starts) - length(starts)) of $(length(all_starts)) segments containing non-finite samples"
    end
    win, win_power = _welch_window(seglen)

    freqs = FFTW.rfftfreq(seglen, 1 / dt)
    nf = length(freqs)
    periodograms = Matrix{Float64}(undef, nf, length(starts))
    for (k, s) in enumerate(starts)
        seg = @view(data[s:(s + seglen - 1)]) .* win
        F = FFTW.rfft(seg)
        @. periodograms[:, k] = (2 * dt / seglen) * abs2(F) / win_power
    end
    # Halve the DC and (for even seglen) Nyquist bins: rfft's one-sided output already carries
    # the full two-sided power at those two bins (they have no negative-frequency mirror to fold
    # in), so the blanket "times 2" above double-counts them specifically.
    periodograms[1, :] ./= 2
    iseven(seglen) && (periodograms[end, :] ./= 2)

    if method === :median
        psd_data = vec(median(periodograms; dims=2)) ./ log(2)
    elseif method === :mean
        psd_data = vec(sum(periodograms; dims=2)) ./ length(starts)
    else
        throw(ArgumentError("Unknown method $method, expected :median or :mean"))
    end

    df = 1 / (seglen * dt)
    return FrequencySeries(psd_data, 0.0, df; name=ts.name === nothing ? nothing : "$(ts.name) PSD")
end

"""
    asd(ts::TimeSeries; kwargs...) -> FrequencySeries

Amplitude spectral density, `sqrt.(psd(ts; kwargs...))` -- same Welch/median convention as `psd`
(see its docstring), just square-rooted.
"""
function asd(ts::TimeSeries; kwargs...)
    p = psd(ts; kwargs...)
    return FrequencySeries(sqrt.(p.data), p.f0, p.df;
                            name=p.name === nothing ? nothing : replace(p.name, " PSD" => " ASD"))
end

end
