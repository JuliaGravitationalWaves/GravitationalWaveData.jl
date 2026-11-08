module Types

export TimeSeries, FrequencySeries, times, frequencies, duration, sample_rate, crop

# ================================================================================================
# TimeSeries: a uniformly-sampled time-domain series anchored to a GPS epoch.
#
# Accessors: `t0`, `dt`, `sample_rate`, `times`, `duration`, `crop`.
# ================================================================================================

"""
    TimeSeries{T,V<:AbstractVector{T}}

A uniformly-sampled time-domain series: `data` (length `N`, real or complex `T`) sampled every
`dt` seconds starting at GPS epoch `t0`. Sample `k` (1-based) occurs at GPS time
`t0 + (k-1) * dt`.

`V` is kept as a free type parameter (rather than hardcoding `Vector{T}`) so views, ranges, and
(eventually) GPU arrays can be stored without an extra copy.

Construct either with `dt` or `sample_rate`:

    TimeSeries(data, t0, dt; name=nothing)
    TimeSeries(data, t0; sample_rate=..., name=nothing)

`name` is optional free-text metadata (e.g. a channel name like `"H1:STRAIN"`); it participates
in `show` but not in arithmetic/equality-of-sampling checks.
"""
struct TimeSeries{T,V<:AbstractVector{T}}
    data::V
    t0::Float64
    dt::Float64
    name::Union{String,Nothing}
end

function TimeSeries(data::AbstractVector, t0::Real, dt::Real;
                     name::Union{AbstractString,Nothing}=nothing)
    dt > 0 || throw(ArgumentError("dt must be positive, got $dt"))
    isfinite(t0) || throw(ArgumentError("t0 must be finite, got $t0"))
    return TimeSeries{eltype(data),typeof(data)}(data, Float64(t0), Float64(dt),
                                                  name === nothing ? nothing : String(name))
end

function TimeSeries(data::AbstractVector, t0::Real; sample_rate::Real,
                     name::Union{AbstractString,Nothing}=nothing)
    sample_rate > 0 || throw(ArgumentError("sample_rate must be positive, got $sample_rate"))
    return TimeSeries(data, t0, 1 / sample_rate; name=name)
end

Base.length(ts::TimeSeries) = length(ts.data)
Base.eltype(::TimeSeries{T}) where {T} = T
Base.size(ts::TimeSeries) = size(ts.data)
Base.firstindex(ts::TimeSeries) = firstindex(ts.data)
Base.lastindex(ts::TimeSeries) = lastindex(ts.data)
Base.iterate(ts::TimeSeries, state...) = iterate(ts.data, state...)

"""
    duration(ts::TimeSeries) -> Float64

Span of the series in seconds, `length(ts) * ts.dt` (i.e. `t0` to `t0 + duration`, half-open --
matches `length(ts) == round(Int, duration(ts) / ts.dt)`).
"""
duration(ts::TimeSeries) = length(ts) * ts.dt

"""
    sample_rate(ts::TimeSeries) -> Float64

`1 / ts.dt`, in Hz.
"""
sample_rate(ts::TimeSeries) = 1 / ts.dt

"""
    times(ts::TimeSeries) -> AbstractRange{Float64}

GPS time of every sample: `t0 .+ dt .* (0:length(ts)-1)`.
"""
times(ts::TimeSeries) = ts.t0 .+ ts.dt .* (0:(length(ts) - 1))

Base.getindex(ts::TimeSeries, i::Integer) = ts.data[i]
function Base.getindex(ts::TimeSeries, r::AbstractUnitRange{<:Integer})
    return TimeSeries(ts.data[r], ts.t0 + ts.dt * (first(r) - 1), ts.dt; name=ts.name)
end

"""
    crop(ts::TimeSeries, t_start=nothing, t_end=nothing) -> TimeSeries

Restrict `ts` to the half-open GPS time window `[t_start, t_end)`, snapping each bound to the
nearest sample (the same rounding `condition_segment` uses). Either bound may be omitted (`nothing`,
the default) to mean "start of the series" / "end of the series" respectively. Throws
`ArgumentError` if the resulting range is empty or falls outside `ts`'s coverage.
"""
function crop(ts::TimeSeries, t_start::Union{Real,Nothing}=nothing,
              t_end::Union{Real,Nothing}=nothing)
    lo = t_start === nothing ? ts.t0 : Float64(t_start)
    hi = t_end === nothing ? ts.t0 + duration(ts) : Float64(t_end)
    lo < hi || throw(ArgumentError("crop requires t_start < t_end, got t_start=$lo, t_end=$hi"))
    i0 = round(Int, (lo - ts.t0) / ts.dt) + 1
    i1 = round(Int, (hi - ts.t0) / ts.dt)
    (i0 >= 1 && i1 <= length(ts) && i1 >= i0) ||
        throw(ArgumentError("crop range [$lo, $hi) falls outside this TimeSeries' coverage " *
                             "[$(ts.t0), $(ts.t0 + duration(ts)))"))
    return ts[i0:i1]
end

function Base.show(io::IO, ts::TimeSeries{T}) where {T}
    print(io, "TimeSeries{$T}($(length(ts)) samples, t0=$(ts.t0), dt=$(ts.dt), ",
          "duration=$(duration(ts))s")
    ts.name === nothing || print(io, ", name=\"$(ts.name)\"")
    print(io, ")")
end

function Base.:(==)(a::TimeSeries, b::TimeSeries)
    return a.t0 == b.t0 && a.dt == b.dt && a.name == b.name && a.data == b.data
end

# ================================================================================================
# FrequencySeries: a uniformly-sampled frequency-domain series.
# ================================================================================================

"""
    FrequencySeries{T,V<:AbstractVector{T}}

A uniformly-sampled frequency-domain series: `data` (length `N`, real or complex `T`, e.g. the
output of `FFTW.rfft`) sampled every `df` Hz starting at frequency `f0` (usually `0.0`). Bin `k`
(1-based) is at frequency `f0 + (k-1) * df`.

    FrequencySeries(data, f0, df; name=nothing)
    FrequencySeries(data, df; name=nothing)   # f0 defaults to 0.0
"""
struct FrequencySeries{T,V<:AbstractVector{T}}
    data::V
    f0::Float64
    df::Float64
    name::Union{String,Nothing}
end

function FrequencySeries(data::AbstractVector, f0::Real, df::Real;
                          name::Union{AbstractString,Nothing}=nothing)
    df > 0 || throw(ArgumentError("df must be positive, got $df"))
    f0 >= 0 || throw(ArgumentError("f0 must be non-negative, got $f0"))
    return FrequencySeries{eltype(data),typeof(data)}(data, Float64(f0), Float64(df),
                                                       name === nothing ? nothing : String(name))
end

function FrequencySeries(data::AbstractVector, df::Real;
                          name::Union{AbstractString,Nothing}=nothing)
    return FrequencySeries(data, 0.0, df; name=name)
end

Base.length(fs::FrequencySeries) = length(fs.data)
Base.eltype(::FrequencySeries{T}) where {T} = T
Base.size(fs::FrequencySeries) = size(fs.data)
Base.firstindex(fs::FrequencySeries) = firstindex(fs.data)
Base.lastindex(fs::FrequencySeries) = lastindex(fs.data)
Base.iterate(fs::FrequencySeries, state...) = iterate(fs.data, state...)

"""
    frequencies(fs::FrequencySeries) -> AbstractRange{Float64}

Frequency (Hz) of every bin: `f0 .+ df .* (0:length(fs)-1)`.
"""
frequencies(fs::FrequencySeries) = fs.f0 .+ fs.df .* (0:(length(fs) - 1))

Base.getindex(fs::FrequencySeries, i::Integer) = fs.data[i]
function Base.getindex(fs::FrequencySeries, r::AbstractUnitRange{<:Integer})
    return FrequencySeries(fs.data[r], fs.f0 + fs.df * (first(r) - 1), fs.df; name=fs.name)
end

"""
    crop(fs::FrequencySeries, f_start=nothing, f_end=nothing) -> FrequencySeries

Restrict `fs` to the half-open frequency window `[f_start, f_end)` Hz, snapping each bound to the
nearest bin. Either bound may be omitted (`nothing`, the default) to mean "first bin" / "last
bin" respectively. Throws `ArgumentError` if the resulting range is empty or falls outside `fs`'s
coverage.
"""
function crop(fs::FrequencySeries, f_start::Union{Real,Nothing}=nothing,
              f_end::Union{Real,Nothing}=nothing)
    lo = f_start === nothing ? fs.f0 : Float64(f_start)
    hi = f_end === nothing ? fs.f0 + fs.df * length(fs) : Float64(f_end)
    lo < hi || throw(ArgumentError("crop requires f_start < f_end, got f_start=$lo, f_end=$hi"))
    i0 = round(Int, (lo - fs.f0) / fs.df) + 1
    i1 = round(Int, (hi - fs.f0) / fs.df)
    (i0 >= 1 && i1 <= length(fs) && i1 >= i0) ||
        throw(ArgumentError("crop range [$lo, $hi) falls outside this FrequencySeries' coverage " *
                             "[$(fs.f0), $(fs.f0 + fs.df * length(fs)))"))
    return fs[i0:i1]
end

function Base.show(io::IO, fs::FrequencySeries{T}) where {T}
    print(io, "FrequencySeries{$T}($(length(fs)) bins, f0=$(fs.f0), df=$(fs.df)")
    fs.name === nothing || print(io, ", name=\"$(fs.name)\"")
    print(io, ")")
end

function Base.:(==)(a::FrequencySeries, b::FrequencySeries)
    return a.f0 == b.f0 && a.df == b.df && a.name == b.name && a.data == b.data
end

# ================================================================================================
# Elementwise arithmetic. Two series can be combined only if they share the same sample grid
# (same `dt`/`df`, same `t0`/`f0`, same length) -- otherwise the bin-by-bin correspondence assumed
# by `+`/`-`/`*`/`/` is meaningless, so we raise a clear `ArgumentError` naming exactly what
# mismatched rather than either erroring opaquely (bare `DimensionMismatch` on length only) or,
# worse, silently truncating/broadcasting past a sampling-rate mismatch.
# ================================================================================================

function _check_compatible(a::TimeSeries, b::TimeSeries)
    a.dt == b.dt ||
        throw(ArgumentError("TimeSeries sampling rates differ: dt=$(a.dt)s ($(sample_rate(a)) Hz) " *
                             "vs dt=$(b.dt)s ($(sample_rate(b)) Hz)"))
    a.t0 == b.t0 ||
        throw(ArgumentError("TimeSeries epochs differ: t0=$(a.t0) vs t0=$(b.t0)"))
    length(a) == length(b) ||
        throw(ArgumentError("TimeSeries lengths differ: $(length(a)) vs $(length(b))"))
    return nothing
end

function _check_compatible(a::FrequencySeries, b::FrequencySeries)
    a.df == b.df ||
        throw(ArgumentError("FrequencySeries bin spacings differ: df=$(a.df) Hz vs df=$(b.df) Hz"))
    a.f0 == b.f0 ||
        throw(ArgumentError("FrequencySeries start frequencies differ: f0=$(a.f0) vs f0=$(b.f0)"))
    length(a) == length(b) ||
        throw(ArgumentError("FrequencySeries lengths differ: $(length(a)) vs $(length(b))"))
    return nothing
end

for op in (:+, :-, :*, :/)
    @eval function Base.$op(a::TimeSeries, b::TimeSeries)
        _check_compatible(a, b)
        return TimeSeries(Base.broadcast($op, a.data, b.data), a.t0, a.dt)
    end
    @eval function Base.$op(a::TimeSeries, b::Number)
        return TimeSeries(Base.broadcast($op, a.data, b), a.t0, a.dt; name=a.name)
    end
    @eval function Base.$op(a::Number, b::TimeSeries)
        return TimeSeries(Base.broadcast($op, a, b.data), b.t0, b.dt; name=b.name)
    end

    @eval function Base.$op(a::FrequencySeries, b::FrequencySeries)
        _check_compatible(a, b)
        return FrequencySeries(Base.broadcast($op, a.data, b.data), a.f0, a.df)
    end
    @eval function Base.$op(a::FrequencySeries, b::Number)
        return FrequencySeries(Base.broadcast($op, a.data, b), a.f0, a.df; name=a.name)
    end
    @eval function Base.$op(a::Number, b::FrequencySeries)
        return FrequencySeries(Base.broadcast($op, a, b.data), b.f0, b.df; name=b.name)
    end
end

end
