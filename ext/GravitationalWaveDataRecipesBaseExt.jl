module GravitationalWaveDataRecipesBaseExt

using GravitationalWaveData: TimeSeries, FrequencySeries, Spectrogram, times, frequencies
using RecipesBase: RecipesBase, @recipe
using Statistics: quantile

# ================================================================================================
# RecipesBase plotting recipes for `TimeSeries`, `FrequencySeries`, and `Spectrogram`. This file is a Julia
# *package extension* (see the `[weakdeps]`/`[extensions]` sections of `Project.toml`): it is only
# loaded -- and only needs to compile -- once both `GravitationalWaveData` and `RecipesBase` are loaded in the
# same session (in practice, whenever a caller also loads a RecipesBase-backed plotting package
# such as `Plots.jl` or `UnicodePlots.jl`, since those depend on `RecipesBase`). `GravitationalWaveData` itself
# carries no hard dependency on `RecipesBase` or any plotting package.
#
# Once loaded, `plot(ts)`, `plot(fs)`, and `plot(spec)` "just work" for `GravitationalWaveData`'s own types with
# no further glue code.
#
# Convention notes (only the load-bearing choices are explained here; see each recipe's own
# docstring for the rest):
#   - Time axes:       x stays absolute GPS seconds, so `vline!(p, [gps])` and friends work, but
#     ticks are labelled relative to a round epoch named in the axis label: ten-digit GPS ticks
#     collide with each other at any ordinary figure width. A caller's `epoch` keyword instead
#     shifts the x data, which GR needs for sub-second windows (see the Spectrogram docstring).
#     Shifting is not the default because series with different windows would then silently
#     misalign on one axis.
#   - TimeSeries:      linear-linear line plot.
#   - FrequencySeries: log-log line plot -- the near-universal convention for PSD/ASD curves in
#     GW data analysis. Set with `-->` (not `:=`) so a caller with non-PSD-shaped data (e.g. a
#     linearly-scaled coherence curve) can override, e.g. `plot(fs; xscale=:identity)`.
#   - Spectrogram:     heatmap, x=time, y=frequency, color=power in dB on a linear colorbar, with
#     color limits from percentiles. Plots.jl's GR backend draws its colorbar ticks itself, ignoring
#     `colorbar_ticks` and `colorbar_formatter`, and puts the colorbar title at a fixed distance from
#     the plot. A log-scaled colorbar over raw strain power therefore gets one crowded label per
#     decade running into the title, whereas dB labels stay short.
#
# Every attribute below is set with `-->` (recipe-default precedence: the *caller's* explicit
# keyword wins), never `:=`, since none of these choices is something a caller must be prevented
# from overriding -- unlike, say, a fixed data transform, there is no GravitationalWaveData-internal invariant
# that would break if a user picks a different colormap, scale, or label.
# ================================================================================================

# GPS time as a plain decimal to the microsecond (a Float64 would print as 1.126259462e9).
function _gps_string(t::Real)
    whole, micro = fldmod(round(Int, t * 1e6), 10^6)
    micro == 0 && return string(whole)
    return string(whole, ".", rstrip(lpad(micro, 6, '0'), '0'))
end

_time_label(epoch::Real) = "Time (s) from GPS $(_gps_string(epoch))"

# Round GPS epoch at or before `t0`, at the power of ten just below `duration`. Tick positions are
# round in absolute GPS seconds, so subtracting an epoch this round leaves round labels, which the
# start of the data would not.
function _round_epoch(t0::Real, duration::Real)
    k = duration > 0 ? max(floor(Int, log10(duration)), -6) : 0
    return floor(t0 / 10.0^k) * 10.0^k
end

# Tick formatter that labels absolute GPS ticks relative to `epoch`.
function _relative_formatter(epoch::Real)
    return function (t)
        value = round(t - epoch; digits=6)
        return isinteger(value) ? string(Int(value)) : string(value)
    end
end

# Power in dB, with zero power (for example, bins a highpass has hard-zeroed) as NaN, which draws
# as a gap and does not count towards the color limits.
_power_db(power::AbstractMatrix) = map(p -> p > 0 ? 10 * log10(p) : NaN, power)

# 1st and 99.9th percentiles of the finite values, so a lone outlier cannot set the color range.
# `nothing` if there is no range to use.
function _color_limits(z::AbstractMatrix)
    finite_values = filter(isfinite, vec(z))
    isempty(finite_values) && return nothing
    low, high = quantile(finite_values, (0.01, 0.999))
    return low < high ? (low, high) : nothing
end

"""
    plot(ts::TimeSeries)

Plot recipe for [`TimeSeries`](@ref) (via `RecipesBase`/`@recipe`): a line plot of sample value
against GPS time (`times(ts)`). Complex-valued data (rare for a time-domain series, but permitted
by the type) is plotted as its magnitude. The x data is absolute GPS time, so `vline!(p, [gps])`
works, while tick labels count from a round epoch named in `xlabel`. `plot(ts; epoch=gps)` instead
plots time relative to `gps`, as for [`Spectrogram`](@ref). `label` uses `ts.name` when set. All
defaults are set with `-->` and so can be overridden by the caller, e.g.
`plot(ts; label="H1", ylabel="Strain")`.

Requires a `RecipesBase`-backed plotting package (e.g. `Plots.jl`) to be loaded to actually render
a plot; `GravitationalWaveData` itself has no hard dependency on one.
"""
@recipe function plot(ts::TimeSeries)
    seriestype --> :path
    epoch = pop!(plotattributes, :epoch, nothing)
    if epoch === nothing
        round_epoch = _round_epoch(ts.t0, ts.dt * (length(ts) - 1))
        xlabel --> _time_label(round_epoch)
        xformatter --> _relative_formatter(round_epoch)
    else
        xlabel --> _time_label(epoch)
    end
    ylabel --> "Amplitude"
    label --> something(ts.name, "")
    y = eltype(ts) <: Complex ? abs.(ts.data) : ts.data
    return collect(times(ts)) .- something(epoch, 0), y
end

"""
    plot(fs::FrequencySeries)

Plot recipe for [`FrequencySeries`](@ref) (via `RecipesBase`/`@recipe`): a line plot of sample
value against frequency (`frequencies(fs)`), defaulting to log-log axes (`xscale`/`yscale` -->
`:log10`) -- the standard convention for PSD/ASD-style curves in GW data analysis. Complex-valued
data (e.g. a raw `FFTW.rfft` output rather than a PSD/ASD) is plotted as its magnitude. `label`
uses `fs.name` when set. All defaults are set with `-->`, so a caller with non-PSD-shaped data can
override, e.g. `plot(fs; xscale=:identity, yscale=:identity)`.

Requires a `RecipesBase`-backed plotting package (e.g. `Plots.jl`) to be loaded to actually render
a plot; `GravitationalWaveData` itself has no hard dependency on one.
"""
@recipe function plot(fs::FrequencySeries)
    seriestype --> :path
    xlabel --> "Frequency (Hz)"
    ylabel --> "Amplitude"
    xscale --> :log10
    yscale --> :log10
    label --> something(fs.name, "")
    y = eltype(fs) <: Complex ? abs.(fs.data) : fs.data
    return collect(frequencies(fs)), y
end

"""
    plot(spec::Spectrogram)

Plot recipe for [`Spectrogram`](@ref) (via `RecipesBase`/`@recipe`): a heatmap with x=GPS time
(`times(spec)`), y=frequency (`frequencies(spec)`), and color=power in dB (`10 log10(spec.data)`,
in the `(freq, time)` row/column order a heatmap recipe expects: `z[i, j]` renders at `(x[j], y[i])`
-- see `Spectrogram`'s own docstring on this axis-order convention). Color limits default to the
1st and 99.9th percentiles of the dB values, so that a few extreme pixels do not wash out the rest;
override with, e.g., `clims=(-60, -20)`. The x data is absolute GPS time, so `vline!(p, [gps])`
works, while tick labels count from a round epoch named in `xlabel`.

`plot(spec; epoch=gps)` instead plots time relative to `gps`, so later series on the same axes
must use the same origin. Use it for windows under a few seconds: the GR backend lacks the precision
to draw a sub-second window at absolute GPS coordinates, and loses the colorbar's color gradient.

Passing `colorbar_scale` plots raw power on that scale instead of dB, with the full data range as
the color limits, e.g. `plot(spec; colorbar_scale=:identity)` for a linear color scale.

Requires a `RecipesBase`-backed plotting package (e.g. `Plots.jl`) to be loaded to actually render
a plot; `GravitationalWaveData` itself has no hard dependency on one.
"""
@recipe function plot(spec::Spectrogram)
    seriestype --> :heatmap
    epoch = pop!(plotattributes, :epoch, nothing)
    if epoch === nothing
        round_epoch = _round_epoch(spec.t0, spec.dt * (size(spec, 2) - 1))
        xlabel --> _time_label(round_epoch)
        xformatter --> _relative_formatter(round_epoch)
    else
        xlabel --> _time_label(epoch)
    end
    ylabel --> "Frequency (Hz)"
    seriescolor --> :viridis
    quantity = something(spec.name, "Power")
    if haskey(plotattributes, :colorbar_scale)
        colorbar_title --> quantity
        z = spec.data
    else
        colorbar_title --> "$quantity (dB)"
        z = _power_db(spec.data)
        limits = _color_limits(z)
        limits === nothing || (clims --> limits)
    end
    return collect(times(spec)) .- something(epoch, 0), collect(frequencies(spec)), z
end

end
