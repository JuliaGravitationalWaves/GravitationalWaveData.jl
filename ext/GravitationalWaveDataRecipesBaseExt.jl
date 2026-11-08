module GravitationalWaveDataRecipesBaseExt

using GravitationalWaveData: TimeSeries, FrequencySeries, Spectrogram, times, frequencies
using RecipesBase: RecipesBase, @recipe

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
#   - TimeSeries:      linear-linear line plot.
#   - FrequencySeries: log-log line plot -- the near-universal convention for PSD/ASD curves in
#     GW data analysis. Set with `-->` (not `:=`) so a caller with non-PSD-shaped data (e.g. a
#     linearly-scaled coherence curve) can override, e.g. `plot(fs; xscale=:identity)`.
#   - Spectrogram:     heatmap, x=time, y=frequency, color=power, with a log-scaled colorbar, since raw STFT
#     power spans many orders of magnitude across frequency. Also set with `-->`, not `:=`.
#
# Every attribute below is set with `-->` (recipe-default precedence: the *caller's* explicit
# keyword wins), never `:=`, since none of these choices is something a caller must be prevented
# from overriding -- unlike, say, a fixed data transform, there is no GravitationalWaveData-internal invariant
# that would break if a user picks a different colormap, scale, or label.
# ================================================================================================

"""
    plot(ts::TimeSeries)

Plot recipe for [`TimeSeries`](@ref) (via `RecipesBase`/`@recipe`): a line plot of sample value
against GPS time (`times(ts)`). Complex-valued data (rare for a time-domain series, but permitted
by the type) is plotted as its magnitude. Default `xlabel`/`ylabel` are generic; `label` uses
`ts.name` when set. All defaults are set with `-->` and so can be overridden by the caller, e.g.
`plot(ts; label="H1", ylabel="Strain")`.

Requires a `RecipesBase`-backed plotting package (e.g. `Plots.jl`) to be loaded to actually render
a plot; `GravitationalWaveData` itself has no hard dependency on one.
"""
@recipe function plot(ts::TimeSeries)
    seriestype --> :path
    xlabel --> "GPS Time (s)"
    ylabel --> "Amplitude"
    label --> something(ts.name, "")
    y = eltype(ts) <: Complex ? abs.(ts.data) : ts.data
    return collect(times(ts)), y
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
(`times(spec)`), y=frequency (`frequencies(spec)`), and color=power (`spec.data`, already stored in
the `(freq, time)` row/column order a heatmap recipe expects: `z[i, j]` renders at `(x[j], y[i])` --
see `Spectrogram`'s own docstring on this axis-order convention). Defaults to a log-scaled colorbar
(`colorbar_scale --> :log10`), since spectrogram power spans many orders of magnitude across
frequency; overridable, e.g. `plot(spec; colorbar_scale=:identity)` for a linear color scale.

Requires a `RecipesBase`-backed plotting package (e.g. `Plots.jl`) to be loaded to actually render
a plot; `GravitationalWaveData` itself has no hard dependency on one.
"""
@recipe function plot(spec::Spectrogram)
    seriestype --> :heatmap
    xlabel --> "GPS Time (s)"
    ylabel --> "Frequency (Hz)"
    colorbar_scale --> :log10
    colorbar_title --> something(spec.name, "Power")
    return collect(times(spec)), collect(frequencies(spec)), spec.data
end

end
