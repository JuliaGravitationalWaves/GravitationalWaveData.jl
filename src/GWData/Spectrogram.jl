module Spectrograms

import DSP
import ..Types: TimeSeries, FrequencySeries, sample_rate, times, frequencies

export Spectrogram, spectrogram

# ================================================================================================
# Spectrogram: a uniformly-sampled time-frequency series -- one FrequencySeries-shaped column of
# power per time bin.
#
# Design choice: a dedicated `Spectrogram{T,M}` type, not a bare
# `Matrix` + separate axis vectors. This mirrors `TimeSeries`/`FrequencySeries` from Types.jl
# exactly: those store `data` alongside `(t0, dt)` / `(f0, df)` rather than materializing a full
# axis vector, and derive `times`/`frequencies` on demand via the exported generic functions of
# those names. A `Spectrogram` is naturally the same idea in two dimensions at once -- it has a
# GPS-anchored time axis (t0, dt) *and* a frequency axis (f0, df), so it stores both axis-defining
# pairs and reuses (by adding methods to) the same generic `times`/`frequencies` functions Types.jl
# already exports, rather than inventing new accessor names. This keeps one API idiom across the
# whole package instead of a special-cased tuple-return for spectrograms only.
# ================================================================================================

"""
    Spectrogram{T,M<:AbstractMatrix{T}}

A uniformly-sampled time-frequency series: `data` is an `nfreq x ntime` matrix (row `i`, column
`j` is the power in frequency bin `i` during time bin `j` -- the same `(freq, time)` axis order
`DSP.spectrogram` uses for its `power` field). Row `i` is at frequency
`f0 + (i-1) * df` (Hz); column `j` is centered at GPS time `t0 + (j-1) * dt` (seconds) -- i.e. `t0`
here is the GPS time of the *center* of the first STFT segment, not the start of the input
`TimeSeries`, matching how `DSP.spectrogram`'s own `time` field is defined (segment-center
times).

    Spectrogram(data, t0, dt, f0, df; name=nothing)

`name` is optional free-text metadata, as with `TimeSeries`/`FrequencySeries`.

Use the existing `times`/`frequencies` generic functions (this module adds methods to the ones
Types.jl exports) to get the GPS time axis / frequency axis as ranges, exactly as for
`TimeSeries`/`FrequencySeries`.
"""
struct Spectrogram{T,M<:AbstractMatrix{T}}
    data::M
    t0::Float64
    dt::Float64
    f0::Float64
    df::Float64
    name::Union{String,Nothing}
end

function Spectrogram(data::AbstractMatrix, t0::Real, dt::Real, f0::Real, df::Real;
                      name::Union{AbstractString,Nothing}=nothing)
    dt > 0 || throw(ArgumentError("dt must be positive, got $dt"))
    df > 0 || throw(ArgumentError("df must be positive, got $df"))
    isfinite(t0) || throw(ArgumentError("t0 must be finite, got $t0"))
    f0 >= 0 || throw(ArgumentError("f0 must be non-negative, got $f0"))
    return Spectrogram{eltype(data),typeof(data)}(data, Float64(t0), Float64(dt), Float64(f0),
                                                   Float64(df),
                                                   name === nothing ? nothing : String(name))
end

Base.size(s::Spectrogram) = size(s.data)
Base.size(s::Spectrogram, d::Integer) = size(s.data, d)
Base.eltype(::Spectrogram{T}) where {T} = T

"""
    times(s::Spectrogram) -> AbstractRange{Float64}

GPS time of the center of every time bin (column): `t0 .+ dt .* (0:ntime-1)`.
"""
times(s::Spectrogram) = s.t0 .+ s.dt .* (0:(size(s, 2) - 1))

"""
    frequencies(s::Spectrogram) -> AbstractRange{Float64}

Frequency (Hz) of every frequency bin (row): `f0 .+ df .* (0:nfreq-1)`.
"""
frequencies(s::Spectrogram) = s.f0 .+ s.df .* (0:(size(s, 1) - 1))

Base.getindex(s::Spectrogram, i::Integer, j::Integer) = s.data[i, j]

"""
    s[:, j] -> FrequencySeries

The `j`th time bin's power spectrum.
"""
function Base.getindex(s::Spectrogram, ::Colon, j::Integer)
    return FrequencySeries(s.data[:, j], s.f0, s.df; name=s.name)
end

"""
    s[i, :] -> TimeSeries

The `i`th frequency bin's power evolution over time (a "time series of one frequency bin",
uniformly sampled at the spectrogram's time-bin spacing `dt`).
"""
function Base.getindex(s::Spectrogram, i::Integer, ::Colon)
    return TimeSeries(s.data[i, :], s.t0, s.dt; name=s.name)
end

function Base.show(io::IO, s::Spectrogram{T}) where {T}
    print(io, "Spectrogram{$T}($(size(s, 1)) freq bins x $(size(s, 2)) time bins, ",
          "t0=$(s.t0), dt=$(s.dt), f0=$(s.f0), df=$(s.df)")
    s.name === nothing || print(io, ", name=\"$(s.name)\"")
    print(io, ")")
end

function Base.:(==)(a::Spectrogram, b::Spectrogram)
    return a.t0 == b.t0 && a.dt == b.dt && a.f0 == b.f0 && a.df == b.df && a.name == b.name &&
           a.data == b.data
end

# ================================================================================================
# `spectrogram(ts::TimeSeries; ...)`: raw (un-averaged) STFT power spectrogram, wrapping
# `DSP.spectrogram` -- one windowed FFT per segment. Welch/median averaging of several
# sub-segments within each column is a separate feature, not provided here.
# ================================================================================================

"""
    spectrogram(ts::TimeSeries; n=length(ts) >> 3, noverlap=n >> 1, nfft=DSP.nextfastfft(n),
                window=DSP.Windows.hanning, onesided=eltype(ts) <: Real,
                name=ts.name) -> Spectrogram

Short-time Fourier transform (STFT) power spectrogram of `ts`: `n`-sample segments, `noverlap`
samples of overlap between consecutive segments, each segment windowed and FFT'd independently
(no averaging across segments within a column -- see the module-level note above for how this
differs from a Welch/median-averaged spectrogram). Thin wrapper around `DSP.spectrogram`.

Defaults match `DSP.spectrogram`'s own (`n = length(ts) >> 3`, 50% overlap, `nfft =
DSP.nextfastfft(n)`, `window = DSP.Windows.hanning`, i.e. `zerophase=false`). This is the
standard spectral-analysis convention (window peak centered in the middle of the segment, as
`scipy.signal.spectrogram` also uses). It is tempting to reach for
`DSP.Windows.hanning(n; zerophase=true)` for this, since scipy's own docs call its own convention
"periodic" -- but DSP.jl's `zerophase=true` means something else entirely: per that function's own
docstring, it centers the window's *peak* at index 1 with the negative half wrapped to the end
(a zero-phase-filter-kernel arrangement in the time domain, not a windowing-for-spectral-analysis
arrangement) -- using it here silently redistributes power across frequency bins by up to several
orders of magnitude at some frequencies.

The returned `Spectrogram`'s GPS time axis is derived from `ts.t0`: `DSP.spectrogram`'s own `time`
field gives segment-center times in seconds *relative to the start of the input signal*, so
`Spectrogram.t0 = ts.t0 + (that first segment-center offset)`.
"""
function spectrogram(ts::TimeSeries; n::Integer=max(length(ts) >> 3, 1),
                      noverlap::Integer=n >> 1,
                      nfft::Integer=DSP.nextfastfft(n),
                      window=DSP.Windows.hanning,
                      onesided::Bool=eltype(ts) <: Real,
                      name::Union{AbstractString,Nothing}=ts.name)
    fs = sample_rate(ts)
    spec = DSP.spectrogram(ts.data, n, noverlap; fs=fs, nfft=nfft, onesided=onesided,
                            window=window)
    freqs = DSP.freq(spec)
    tvals = DSP.time(spec)

    f0 = Float64(first(freqs))
    df = Float64(step(freqs))
    t0 = ts.t0 + Float64(first(tvals))
    dt = Float64(step(tvals))

    return Spectrogram(DSP.power(spec), t0, dt, f0, df; name=name)
end

end
