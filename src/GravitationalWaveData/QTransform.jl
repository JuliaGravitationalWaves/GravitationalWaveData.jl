module QTransform

using FFTW
using Statistics: median, mean
using ..Types: TimeSeries, FrequencySeries, duration, sample_rate

export normalized_fft, QTileParams, qtile_window, qtile_data_indices, qtile_padding,
       qtile_energy, q_values, tiling_frange, plane_frange, plane_frequencies,
       QRow, QScan, QPeak, energy_series, q_plane_transform, q_scan, q_transform, qscan_peak

# ================================================================================================
# Constant-Q transform, following the Omega-pipeline / Chatterji-thesis algorithm. The tiling,
# windowing and indexing conventions are ported from gwpy's implementation
# (`gwpy.signal.qtransform`: `QTiling`/`QPlane`/`QTile`), because the thesis alone does not pin
# them down -- in particular a `+1` index offset inherited from the original MATLAB Omega code
# and the exact "bi-square" window normalization. Intermediate quantities (window values, data
# indices, a single tile's energy series) are checked against reference values in
# `test/qtransform.jl`.
#
# Contents: the single-Q tile/plane computation (`QTileParams`, `qtile_window`,
# `qtile_data_indices`, `qtile_padding`, `qtile_energy`, `q_plane_transform`), a multi-Q search
# (`q_scan`/`q_transform`: loop over log-spaced Q values, keep the plane containing the loudest
# tile) and peak read-out (`qscan_peak`).
#
# Not implemented: spline re-gridding of the irregular per-row time axes onto a common
# time/frequency grid for pixel plots. That step is presentation-only (it does not change which
# tile is loudest or its energy/time/frequency); a `QScan`'s `rows` expose the raw grid so it can
# be added later without recomputing anything here. Also not implemented: the false-alarm-rate
# estimate.
#
# Input convention: the core functions consume an already-Fourier-transformed spectrum (via
# `normalized_fft` -- NOT whitened). To get a whitened Q-transform, `whiten` the `TimeSeries`
# first and pass the result to `q_transform`.
# ================================================================================================

const DEFAULT_QRANGE = (4.0, 64.0)
const DEFAULT_FRANGE = (0.0, Inf)
const DEFAULT_MISMATCH = 0.2

_deltam(mismatch::Real) = 2 * sqrt(mismatch / 3)

# ------------------------------------------------------------------------------------------------
# Step 0: normalized one-sided FFT, matching gwpy's `TimeSeries.fft()` (NOT the same normalization
# as a plain `rfft` -- divide by N, then double every bin except DC so the result is a "one-sided
# amplitude spectrum": frequency bin k (1-based here, k=1 is DC) is at frequency `(k-1)/duration`).
# ------------------------------------------------------------------------------------------------

"""
    normalized_fft(ts::TimeSeries) -> FrequencySeries{ComplexF64}

One-sided, amplitude-normalized FFT: `X = rfft(data) / N`, then every bin except DC is doubled. This
(not a bare `rfft`) is the `FrequencySeries` convention the rest of this module's tile computations
assume; the per-tile normalization (`qtile_window`) is calibrated against it.
"""
function normalized_fft(ts::TimeSeries)
    n = length(ts)
    X = ComplexF64.(rfft(ts.data)) ./ n
    if length(X) > 1
        X[2:end] .*= 2
    end
    df = 1 / duration(ts)
    return FrequencySeries(X, 0.0, df)
end

# ------------------------------------------------------------------------------------------------
# Step 1: Q values for a (qlow, qhigh) range -- gwpy's `QTiling._iter_qs`.
# ------------------------------------------------------------------------------------------------

"""
    q_values(qrange, dur, mismatch) -> Vector{Float64}

Logarithmically-spaced Q values spanning `qrange = (qlow, qhigh)`, dense enough that neighbouring
planes' tiles mismatch by no more than `mismatch`.
"""
function q_values(qrange::Tuple{<:Real,<:Real}, dur::Real, mismatch::Real)
    dm = _deltam(mismatch)
    cumum = log(qrange[2] / qrange[1]) / sqrt(2)
    nplanes = max(ceil(Int, cumum / dm), 1)
    delta_q = cumum / nplanes
    return [qrange[1] * exp(sqrt(2) * delta_q * (i + 0.5)) for i in 0:(nplanes - 1)]
end

# ------------------------------------------------------------------------------------------------
# Step 2: overall (fmin, fmax) resolution across the whole Q range -- gwpy's `QTiling.__init__`.
# A zero low bound / infinite high bound get resolved to concrete Hz values ONCE, using the
# extreme Q in the range, and that same (fmin, fmax) is then reused for every plane (each `QPlane`
# would only re-derive fmin/fmax from a literal 0/Inf sentinel, which by this point has already
# been replaced).
# ------------------------------------------------------------------------------------------------

"""
    tiling_frange(qrange, frange, dur, fs, mismatch) -> (fmin, fmax)

Resolve a `(fmin, fmax)` frequency range for a multi-Q scan: `fmin == 0` is replaced by
`50*max(q)/(2*pi*dur)`, `fmax == Inf` (or too high for the Q range) is replaced/clamped to
`fs/2/(1+sqrt(11)/min(q))`.
"""
function tiling_frange(qrange::Tuple{<:Real,<:Real}, frange::Tuple{<:Real,<:Real}, dur::Real,
                        fs::Real, mismatch::Real)
    qlist = q_values(qrange, dur, mismatch)
    fmin, fmax = float(frange[1]), float(frange[2])
    if fmin == 0
        fmin = 50 * maximum(qlist) / (2 * pi * dur)
    end
    maxf = fs / 2 / (1 + sqrt(11) / minimum(qlist))
    if isinf(fmax)
        fmax = maxf
    elseif fmax > maxf
        @warn "upper frequency of $fmax Hz is too high for the given Q range, resetting to $maxf Hz"
        fmax = maxf
    end
    return (fmin, fmax)
end

"""
    plane_frange(q, dur, frange, fs) -> (fmin, fmax)

Single-`Q`-plane frange resolution (only relevant when a plane is used standalone with an
unresolved `frange`, e.g. `(0, Inf)` -- `q_scan`/`q_transform` instead resolve `frange` once via
`tiling_frange` and reuse it across every `Q`).
"""
function plane_frange(q::Real, dur::Real, frange::Tuple{<:Real,<:Real}, fs::Real)
    fmin, fmax = float(frange[1]), float(frange[2])
    fmin == 0 && (fmin = 50 * q / (2 * pi * dur))
    isinf(fmax) && (fmax = fs / 2 / (1 + 1 / (q / sqrt(11))))
    return (fmin, fmax)
end

# ------------------------------------------------------------------------------------------------
# Step 3: central frequencies of the tiles in one Q-plane -- gwpy's `QPlane._iter_frequencies`.
# ------------------------------------------------------------------------------------------------

"""
    plane_frequencies(q, frange, dur, mismatch) -> Vector{Float64}

Central frequencies of the tile rows in a single `Q`-plane spanning `frange = (fmin, fmax)`,
log-spaced densely enough for `mismatch`, each snapped down to the nearest multiple of `1/dur`
(the frequency-bin spacing of a `dur`-second segment).
"""
function plane_frequencies(q::Real, frange::Tuple{<:Real,<:Real}, dur::Real, mismatch::Real)
    minf, maxf = frange
    dm = _deltam(mismatch)
    fcum_mismatch = log(maxf / minf) * sqrt(2 + q^2) / 2
    nfreq = max(1, ceil(Int, fcum_mismatch / dm))
    fstep = fcum_mismatch / nfreq
    fstepmin = 1 / dur
    freqs = Float64[]
    last = nothing
    for i in 0:(nfreq - 1)
        this = floor(minf * exp(2 / sqrt(2 + q^2) * (i + 0.5) * fstep) / fstepmin) * fstepmin
        if last === nothing || this != last
            push!(freqs, this)
        end
        last = this
    end
    return freqs
end

# ------------------------------------------------------------------------------------------------
# Step 4: a single (Q, frequency) tile -- gwpy's `QTile`. This is the numerically riskiest part:
# `qtile_data_indices`'s `+1` offset (beyond the usual 0-based->1-based conversion) is a literal
# quirk inherited from gwpy's Python source (itself presumably inherited from the original MATLAB
# Omega code) -- ported verbatim rather than "corrected", and cross-checked directly against
# gwpy's own `QTile.get_data_indices()`/`get_window()` output in `test/qtransform.jl`.
# ------------------------------------------------------------------------------------------------

"""
    QTileParams(q, frequency, duration, sampling, mismatch)

Parameters of a single Q-transform tile row: a fixed `Q`, central `frequency` (Hz), analysis segment
`duration` (s) and `sampling` rate (Hz), and tiling `mismatch`. The tile computations are the free
functions below (`qtile_window`, `qtile_data_indices`, `qtile_padding`, `qtile_energy`).
"""
struct QTileParams
    q::Float64
    frequency::Float64
    duration::Float64
    sampling::Float64
    mismatch::Float64
end

qprime(t::QTileParams) = t.q / sqrt(11)

function _ntiles(t::QTileParams)
    tcum_mismatch = t.duration * 2 * pi * t.frequency / t.q
    x = tcum_mismatch / _deltam(t.mismatch)
    return Int(2^ceil(log(x) / log(2)))
end

_windowsize(t::QTileParams) = 2 * floor(Int, t.frequency / qprime(t) * t.duration) + 1

_window_offsets(t::QTileParams) = (-div(_windowsize(t) - 1, 2)):(div(_windowsize(t) - 1, 2))

"""
    qtile_window(t::QTileParams) -> Vector{Float64}

The "bi-square" frequency-domain window `(1 - x^2)^2`, normalized so the tile's IFFT energy has
the right absolute scale.
"""
function qtile_window(t::QTileParams)
    offsets = collect(_window_offsets(t))
    wfrequencies = offsets ./ t.duration
    xfrequencies = wfrequencies .* (qprime(t) / t.frequency)
    nrm = _ntiles(t) / (t.duration * t.sampling) * sqrt(315 * qprime(t) / (128 * t.frequency))
    return (1.0 .- xfrequencies .^ 2) .^ 2 .* nrm
end

"""
    qtile_data_indices(t::QTileParams) -> UnitRange-like Vector{Int}

1-based Julia indices into a `normalized_fft`-style one-sided `FrequencySeries` selecting the
frequency bins this tile's window applies to: `round(offset + 1 + frequency*duration) + 1`. The
final `+1` converts 0-based to 1-based indexing; the inner `+1` is an offset inherited from the
original Omega pipeline, kept as-is so results match other Omega-derived implementations.
"""
function qtile_data_indices(t::QTileParams)
    offsets = collect(_window_offsets(t))
    py0based = round.(Int, offsets .+ 1 .+ t.frequency * t.duration)
    return py0based .+ 1
end

"""
    qtile_padding(t::QTileParams) -> (left, right)

Zero-padding (in bins) applied to a tile's windowed frequency-domain data before the IFFT, so the
padded length equals `_ntiles(t)` (a power of 2).
"""
function qtile_padding(t::QTileParams)
    pad = _ntiles(t) - _windowsize(t)
    return (div(pad - 1, 2), div(pad + 1, 2))
end

"""
    qtile_energy(fseries_data::AbstractVector{<:Complex}, t::QTileParams) -> (energy, dt)

Un-normalized-by-median/mean tile energy: window+select the relevant bins of `fseries_data` (a
`normalized_fft`-convention one-sided spectrum), zero-pad, `ifftshift`, inverse-FFT, and take
`abs2`. Returns the energy vector and its sample spacing `dt = duration/ntiles`. The median/mean
normalization is applied afterwards by `q_plane_transform`.
"""
function qtile_energy(fseries_data::AbstractVector{<:Complex}, t::QTileParams)
    idx = qtile_data_indices(t)
    (minimum(idx) >= 1 && maximum(idx) <= length(fseries_data)) ||
        throw(ArgumentError("Q-tile (q=$(t.q), f=$(t.frequency) Hz) needs frequency-series " *
                             "indices outside the available range [1, $(length(fseries_data))] " *
                             "-- narrow `frange` or increase the segment duration/sample rate"))
    windowed = fseries_data[idx] .* qtile_window(t)
    n = _ntiles(t)
    left, _ = qtile_padding(t)
    padded = zeros(ComplexF64, n)
    padded[(left + 1):(left + length(windowed))] .= windowed
    td = ifft(ifftshift(padded))
    energy = abs2.(td)
    dt = t.duration / n
    return energy, dt
end

_normalize_energy(energy::Vector{Float64}, norm::Symbol) =
    norm === :none  ? energy :
    norm === :median ? energy ./ median(energy) :
    norm === :mean   ? energy ./ mean(energy) :
    throw(ArgumentError("norm must be :median, :mean, or :none, got $norm"))

# ------------------------------------------------------------------------------------------------
# Step 5: a full Q-plane (all frequency rows at one Q) and a multi-Q scan.
# ------------------------------------------------------------------------------------------------

"""
    QRow(frequency, energy, t0, dt)

One frequency row's (normalized) energy time series from a Q-plane transform: `energy[k]` is the
energy at GPS time `t0 + (k-1)*dt`.
"""
struct QRow
    frequency::Float64
    energy::Vector{Float64}
    t0::Float64
    dt::Float64
end

"""
    energy_series(row::QRow) -> TimeSeries

Wrap a `QRow`'s energy vector as a `TimeSeries` for `times`/`crop`/etc.
"""
energy_series(row::QRow) = TimeSeries(row.energy, row.t0, row.dt)

"""
    QScan(q, frange, rows)

The full energy grid for one Q-plane: `rows` (one per tile frequency), each with its own time
resolution (rows near higher frequencies get finer time resolution). No re-gridding onto a common
time/frequency axis is done; interpolate the rows yourself for display.
"""
struct QScan
    q::Float64
    frange::Tuple{Float64,Float64}
    rows::Vector{QRow}
end

"""
    q_plane_transform(fseries_data, q, frange, dur, fs, epoch; mismatch=0.2, norm=:median) -> QScan

Compute the full energy grid for a single `Q`-plane: one `qtile_energy` row per frequency from
`plane_frequencies`, normalized by `norm` (`:median`, `:mean`, or `:none`).
"""
function q_plane_transform(fseries_data::AbstractVector{<:Complex}, q::Real,
                            frange::Tuple{<:Real,<:Real}, dur::Real, fs::Real, epoch::Real;
                            mismatch::Real=DEFAULT_MISMATCH, norm::Symbol=:median)
    rows = QRow[]
    for freq in plane_frequencies(q, frange, dur, mismatch)
        tile = QTileParams(q, freq, dur, fs, mismatch)
        energy, dt = qtile_energy(fseries_data, tile)
        energy = _normalize_energy(energy, norm)
        push!(rows, QRow(freq, energy, epoch, dt))
    end
    return QScan(q, frange, rows)
end

"""
    QPeak(time, frequency, energy, snr, q)

The loudest tile of a `QScan` (or the loudest across a multi-Q search): its GPS `time`,
`frequency` (Hz), `energy` (normalized, e.g. by median if `norm=:median`), `snr = sqrt(2*energy)`,
and the `q` of the plane it came from.
"""
struct QPeak
    time::Float64
    frequency::Float64
    energy::Float64
    snr::Float64
    q::Float64
end

"""
    qscan_peak(scan::QScan; search=nothing) -> QPeak

Find the loudest tile in `scan` (the maximum-energy sample across every row), optionally
restricted to a `search = (t_lo, t_hi)` GPS time window.
"""
function qscan_peak(scan::QScan; search::Union{Nothing,Tuple{<:Real,<:Real}}=nothing)
    best = QPeak(0.0, 0.0, 0.0, 0.0, scan.q)
    for row in scan.rows
        n = length(row.energy)
        idxs = 1:n
        if search !== nothing
            lo, hi = search
            idxs = filter(i -> lo <= row.t0 + (i - 1) * row.dt < hi, idxs)
            isempty(idxs) && continue
        end
        maxidx = idxs[argmax(view(row.energy, idxs))]
        e = row.energy[maxidx]
        if e > best.energy
            t = row.t0 + (maxidx - 1) * row.dt
            best = QPeak(t, row.frequency, e, sqrt(2 * e), scan.q)
        end
    end
    return best
end

"""
    q_scan(fseries_data, dur, fs, epoch; qrange=(4,64), frange=(0,Inf), mismatch=0.2, norm=:median,
           search=nothing) -> (QScan, QPeak)

Multi-`Q` search: compute a `q_plane_transform` at each of `q_values(qrange, dur, mismatch)`, and
return the plane containing the single loudest tile (and that tile's `QPeak`). The
most significant Q is chosen this way rather than by merging energies across Q.
"""
function q_scan(fseries_data::AbstractVector{<:Complex}, dur::Real, fs::Real, epoch::Real;
                 qrange::Tuple{<:Real,<:Real}=DEFAULT_QRANGE,
                 frange::Tuple{<:Real,<:Real}=DEFAULT_FRANGE,
                 mismatch::Real=DEFAULT_MISMATCH, norm::Symbol=:median,
                 search::Union{Nothing,Tuple{<:Real,<:Real}}=nothing)
    resolved_frange = tiling_frange(qrange, frange, dur, fs, mismatch)
    best_scan = nothing
    best_peak = QPeak(0.0, 0.0, -Inf, 0.0, 0.0)
    for q in q_values(qrange, dur, mismatch)
        scan = q_plane_transform(fseries_data, q, resolved_frange, dur, fs, epoch;
                                  mismatch=mismatch, norm=norm)
        peak = qscan_peak(scan; search=search)
        if peak.energy > best_peak.energy
            best_scan = scan
            best_peak = peak
        end
    end
    return best_scan, best_peak
end

"""
    q_transform(ts::TimeSeries; qrange=(4,64), frange=(0,Inf), mismatch=0.2, norm=:median,
                search=nothing) -> (QScan, QPeak)

Compute the constant-Q transform of `ts` (via `normalized_fft`; the input is not whitened, so call
`whiten` first if needed) and return the loudest Q-plane's full energy grid plus its peak tile.
"""
function q_transform(ts::TimeSeries; qrange::Tuple{<:Real,<:Real}=DEFAULT_QRANGE,
                      frange::Tuple{<:Real,<:Real}=DEFAULT_FRANGE,
                      mismatch::Real=DEFAULT_MISMATCH, norm::Symbol=:median,
                      search::Union{Nothing,Tuple{<:Real,<:Real}}=nothing)
    fseries = normalized_fft(ts)
    return q_scan(fseries.data, duration(ts), sample_rate(ts), ts.t0;
                  qrange=qrange, frange=frange, mismatch=mismatch, norm=norm, search=search)
end

end
