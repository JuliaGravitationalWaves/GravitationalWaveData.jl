# GWData.jl conventions

This page collects the non-obvious numerical and physical conventions GWData.jl commits to. Each
is also documented, with its exact provenance, at its point of use in `src/GWData/`; when the two
disagree, the source comment is authoritative.

The conventions follow established GW-analysis practice (`scipy.signal`, `gwpy`, LALSuite) rather
than inventing new ones, so results are directly comparable with other pipelines. The measured
agreement with those reference implementations is summarized under [Validation](#Validation).

## PSD/ASD normalization

`psd`/`asd` (`src/GWData/PSD.jl`) estimate a one-sided power/amplitude spectral density from a
`TimeSeries` using median-Welch averaging.

- **Segmentation**: the series is split into `seglen`-sample segments with fractional `overlap`
  (default `0.5`), starting at 1-based indices `1, 1+step, 1+2*step, ...` where
  `step = round(seglen * (1 - overlap))`. Too-short input raises `ArgumentError` rather than
  silently returning zero segments. `csd`/`coherence` (`src/GWData/Coherence.jl`) reuse this exact
  segmentation via the shared internal helper `_welch_starts`, so the two never drift apart.
- **Window**: each segment is multiplied by a Hann window (`DSP.Windows.hanning`, the *symmetric*
  convention — endpoints exactly zero, not the *periodic*/DFT-even convention some FFT-adjacent
  code uses), i.e. `0.5 * (1 - cos(2*pi*n/(N-1)))`. This is also `scipy.signal.welch`'s default.
- **Scaling**: each windowed segment's periodogram is `(2 * dt / seglen) * abs2(rfft(segment)) /
  win_power`, where `win_power = mean(window.^2)` corrects for the window's power loss. This is the
  standard *density* normalization (`scipy.signal.welch`'s `scaling="density"`): `∫ psd(f) df` over
  `[0, Nyquist]` approximates the time series' variance.
- **DC/Nyquist bin-halving**: the blanket `2 *` factor above assumes every one-sided bin has a
  mirrored negative-frequency counterpart to fold in. The DC bin (always) and the Nyquist bin (only
  when `seglen` is even) have no such mirror, so they are divided by 2 immediately after — otherwise
  their power would be double-counted. This is a common, easy-to-miss detail in one-sided PSD
  estimators; get it wrong and low-frequency/Nyquist-adjacent bins are off by a factor of 2.
- **Averaging across segments**: `method=:median` (the default, the standard LIGO convention) takes the per-frequency-bin *median* periodogram across segments — far more robust to
  a handful of segments being contaminated by loud transient glitches than a mean would be — and
  then divides by `ln(2)` to debias it (the median of a chi-squared/exponential-distributed
  periodogram systematically underestimates the mean by very close to that factor, asymptotically
  exact for the segment counts this package sees in practice, e.g. O(500) segments for a 1024s
  off-source stretch at 4s/50% overlap). `method=:mean` takes the plain per-bin mean and skips the
  `ln(2)` correction entirely (it isn't needed for a mean estimator).
- **`asd`** is exactly `sqrt.(psd(...))`.
- **Units**: PSD is strain²/Hz; ASD is strain/√Hz (or the equivalent unit² /Hz, unit/√Hz for a
  non-strain input series).

### CSD / coherence normalization

`csd`/`coherence` (`src/GWData/Coherence.jl`) generalize the same Welch machinery to a pair of
channels, following `scipy.signal.csd`/`scipy.signal.coherence`:

- Same segmentation/window as `psd` (reusing `_welch_starts`/`_welch_window` directly), same
  DC/Nyquist bin-halving, same density scaling generalized to a cross term:
  `Pxy(f) = mean_over_segments((2*dt/seglen) * conj(rfft(x)) .* rfft(y) / win_power)`.
- **Always mean-averaged, never median**, because a median of *complex* values isn't the same
  well-defined robust statistic a median of real, non-negative PSD values is.
- **Coherence** is `|Pxy(f)|^2 / (Pxx(f) * Pyy(f))`, real-valued, in `[0, 1]`, with `Pxx`/`Pyy`
  computed via `psd(...; method=:mean)` (mean, not median — so the denominator uses the same
  averaging convention as the numerator's `Pxy`) — matching `scipy.signal.coherence`'s own
  `average="mean"` default.
- **Detrending gap, deliberately not closed**: `scipy.signal.csd`/`coherence` detrend (subtract the
  per-segment mean) by default; this package's `psd` does not detrend, and `csd`/`coherence`
  follow suit so numerator and denominator use the same convention. The effect on real data is
  negligible (see [Validation](#Validation)).

## Whitening

`whiten` (`src/GWData/Filtering.jl`) whitens a `TimeSeries` against a supplied PSD/ASD using
**inverse spectrum truncation** (the algorithm used by `gwpy`'s `TimeSeries.whiten()`), *not* the
naive "FFT, divide by ASD, inverse FFT" approach, which rings badly at segment edges.

The algorithm, in order:

1. **Interpolate** the supplied PSD/ASD (linear interpolation) onto the frequency grid with spacing
   `1/duration(ts)` — so a single off-source ASD estimate can be reused to whiten several
   different on-source segments without needing to already match their exact resolution.
2. **Design a time-domain FIR filter** approximating `1/ASD(f)`: transform to the time domain,
   truncate to `fduration` seconds (default 2s) of taps, and taper the truncation with a window
   (default Hann) so the filter doesn't have unbounded support or a sharp cutoff discontinuity.
   If `highpass` is given, the filter's response below that frequency is hard-zeroed first.
3. **Convolve** that FIR filter with the mean-detrended input data using the **overlap-save**
   method, with the input's own edges tapered — never a bare frequency-domain divide.
4. **Scale** the result by `sqrt(2 * dt)`: on stationary Gaussian input, whitened output has zero
   mean and unit variance.

**Edge effects**: a segment of length `fduration/2` at each end of the *output* is corrupted by the
FIR filter's settle-in transient (tapered rather than exact). Crop it off if that matters
downstream.

Background theory: arXiv:gr-qc/0509116.

## Filtering

`bandpass`/`highpass`/`lowpass`/`notch` (`src/GWData/Filtering.jl`) use a fixed-order (default 4)
Butterworth design applied zero-phase via `DSP.filtfilt`. Filter order is explicit, not
auto-selected from a ripple/attenuation specification, so results differ from tools that design
the order automatically (for example `gwpy`'s `iirdesign`-based `.bandpass()`).

## GPS time

GW data is timestamped in **GPS seconds**: a continuous SI-second count from the GPS epoch
(1980-01-06T00:00:00 UTC) with **no leap seconds** inserted after that epoch — GPS time never
stops, skips, or repeats a second. This is why GW pipelines universally prefer it over UTC as an
internal timestamp: it's monotonic and gap-free.

`gps_to_utc`/`utc_to_gps` (`src/GWData/GPSTime.jl`) convert between GPS seconds and UTC
`DateTime`s, **leap-second aware** via a historical table (`LEAP_SECOND_TABLE`) sourced from IERS
Bulletin C, current through Bulletin C 71 (July 2025; most recent actual leap second inserted
2016-12-31, i.e. `GPS - UTC = 18s` from 2017-01-01 onward). This is a real table lookup, not a
single hardcoded offset constant — a hardcoded constant silently goes stale the next time IERS
announces a new leap second, and there is no way to detect that staleness at runtime. If IERS
announces a new leap second (Bulletin C is issued twice yearly), extend `LEAP_SECOND_TABLE` with
the new entry; both conversion functions derive everything from that one table, so nothing else
needs to change.

A related but independently-sourced leap-second table lives in `src/GWData/SiderealTime.jl`
(`LEAPS`), ported verbatim from LALSuite's own `XLALLeapSeconds.h` rather than IERS directly — this
is deliberate: `gmst` (Greenwich Mean Sidereal Time) must match LALSuite's exact numerical
convention bit-for-bit for interoperability with LALSuite-based pipelines, which means reproducing LALSuite's *specific* leap-second bookkeeping rather than a more
astronomically "correct" modern one. See that module's header comment for the full rationale. The
two tables agree on every leap second either one is aware of; they are kept as two separate tables
because they serve two different faithfulness goals (general UTC correctness vs. LALSuite
bit-for-bit reproduction), not out of oversight.

Caveat shared by essentially every GPS/UTC converter: `Dates.DateTime` cannot represent the leap
second itself (`23:59:60`), since Julia's `Dates` module has no notion of leap seconds — only
ordinary instants strictly before or after an insertion are meaningful inputs to `utc_to_gps`.

### Segment / SegmentList: half-open intervals

`Segment`/`SegmentList` (`src/GWData/Segments.jl`) represent GPS-time intervals for data-quality
vetoing, following the standard LIGO segment-list convention:

- Every interval is **half-open**, `[start, stop)` — `start` is included, `stop` is not. This is
  the same half-open convention `TimeSeries.crop`/`FrequencySeries.crop` use (see below), kept
  consistent across the package so a boundary sample is never double-counted or dropped depending
  on which API touched it.
- A `SegmentList` is always **sorted and coalesced** (overlapping-or-touching segments merged) as a
  maintained class invariant — every `SegmentList` a caller holds is already in canonical form, no
  separate "did you remember to coalesce?" step needed before, say, calling `union`/`intersect`.
- **Touching segments merge, but don't "overlap"**: `overlaps(a, b)` is `false` when `a.stop ==
  b.start` (no shared point in a half-open interval), but `coalesce_segments` (and `SegmentList`
  construction) still merges such touching segments — the more permissive "adjacent enough to
  combine" notion segment-list libraries conventionally use for canonicalization.
- `gaps(sl)` returns the complement of `sl`; without an explicit `bound::Segment`, only *internal*
  gaps (strictly between the first and last segment) are meaningful, since there's no well-defined
  gap before/after the list's own extent otherwise. Pass `bound` to compute the complement within a
  specific analysis window, including any leading/trailing gap against that window's edges.

## TimeSeries / FrequencySeries / Spectrogram: axes and arithmetic

- **TimeSeries**: sample `k` (1-based) is at GPS time `t0 + (k-1)*dt`. **FrequencySeries**: bin `k`
  (1-based) is at frequency `f0 + (k-1)*df` (`f0` defaults to `0.0`). **Spectrogram**: stores both
  axis pairs at once — `data` is an `nfreq x ntime` matrix (row = frequency bin, column = time
  bin), the same `(freq, time)` order `DSP.spectrogram`'s `power` field uses. (Some other tools
  store `(time, freq)`; transpose when exchanging data with them.) A `Spectrogram`'s `t0` is the GPS
  time of the *center* of its first STFT segment (not the start of the input `TimeSeries`),
  matching how `DSP.spectrogram` labels its time bins.
- **`crop`** (on `TimeSeries`/`FrequencySeries`) restricts to a half-open `[start, end)` window,
  snapping each bound to the nearest sample/bin — the same half-open convention `Segment` uses.
  Either bound may be omitted (`nothing`) to mean "start of the series"/"end of the series".
- **Arithmetic (`+ - * /`) between two series** is only defined when they share the *exact* same
  sample grid: same `dt`/`t0`/length for `TimeSeries`, same `df`/`f0`/length for `FrequencySeries`.
  A mismatch raises a specific `ArgumentError` naming exactly which of the three checks
  (rate/epoch/length) failed, rather than either an opaque bare `DimensionMismatch` (length-only)
  or, worse, silently broadcasting/truncating across a sampling-rate mismatch that would make the
  bin-by-bin correspondence meaningless. A series combined with a plain `Number` is unaffected by
  this check (broadcasts elementwise as usual).

## Q-transform

`q_transform`/`q_scan` (`src/GWData/QTransform.jl`) implement the Omega-pipeline constant-Q
transform (Chatterji 2005), with the tiling, bi-square window normalization and index conventions
of `gwpy`'s `gwpy.signal.qtransform`:

- **Input is not whitened.** `q_transform` FFTs its input with `normalized_fft`
  (`rfft(data) / N`, one-sided) and tiles that directly. For a whitened Q-scan, call `whiten`
  first and pass the result.
- **Multi-Q search keeps one plane.** Q values are log-spaced across `qrange` at the given
  `mismatch`, and the result is the single plane containing the loudest tile, not a merge
  across planes.
- **Rows have their own time grids.** Each `QRow` in a `QScan` has its own `t0`/`dt` (higher
  frequencies get finer time resolution). No spline re-gridding onto a common time/frequency grid
  is done; interpolate yourself for display (see the Q-transform tutorial).
- **Energy normalization**: each row's energy is divided by its median (`norm=:median`, default)
  or mean (`norm=:mean`), and `snr = sqrt(2 * energy)`.
- **Not implemented**: false-alarm-rate estimation.

## Reference design curves (aLIGO / AdvVirgo)

`aligo_design_psd`/`aligo_design_asd`/`avirgo_design_psd`/`avirgo_design_asd`
(`src/GWData/ReferenceCurves.jl`) are ported verbatim from LALSimulation's `LALSimNoisePSD.c`, not
re-derived. One naming trap worth flagging for anyone extending this module: the well-known simple
closed-form "Ajith fit" some other codebases call `aLIGOZeroDetHighPower` is **not**
what LALSimulation's function of the same name actually computes — the real function composes a
full Buonanno & Chen (2001) quantum-noise model with a suspension/coating thermal-noise fit. This
module ports that actual composition, since that's what it's cross-validated against. In the
`FrequencySeries`-returning form, the DC bin (`f=0`) is set to exactly `0.0` rather than evaluated
(both underlying formulas are singular there) — matching LALSimulation's own `XLALSimNoisePSD`
convention; the plain-vector form has no such guard and will return whatever IEEE754 NaN/Inf the
raw formula produces at `f=0`, same as calling the LALSimulation function directly would.

## Further reading

Every convention above is documented in more depth, with exact formulas and full source citations,
at its point of use: `src/GWData/PSD.jl`, `src/GWData/Coherence.jl`, `src/GWData/Filtering.jl`,
`src/GWData/GPSTime.jl`, `src/GWData/SiderealTime.jl`, `src/GWData/Segments.jl`,
`src/GWData/Types.jl`, `src/GWData/Spectrogram.jl`, `src/GWData/QTransform.jl`, and
`src/GWData/ReferenceCurves.jl`.

## Validation

Each module was cross-checked numerically against an established reference implementation on the
same inputs, using real GW150914 data wherever relevant. The comparison scripts live in
`examples/validation/` (see its README for how to run them; they need Linux or macOS, since
LALSuite has no Windows build). The spot values they produced are kept as regression tests under
`test/`.

| Function(s) | Reference | Agreement |
|:--|:--|:--|
| `psd`/`asd` | `gwpy` `.psd()` (median-Welch), 1024 s H1 off-source | median rel. diff 0.14%, max 1.14% |
| `csd`/`coherence` | `scipy.signal.csd`/`coherence` via `gwpy`, H1/L1 off-source | CSD median rel. diff 0.0093%, max 0.70%; coherence abs. diff max 2.3e-5 |
| `whiten` | `gwpy` `.whiten()` | interior rel. RMS diff 4.0e-12 |
| `lowpass`/`highpass` | `gwpy` Butterworth filters | rel. RMS diff 0.03% / ~2% |
| `bandpass` | `gwpy` `.bandpass()` | ~24% rel. RMS diff, from `gwpy` auto-selecting a higher filter order (see [Filtering](#Filtering)) |
| `resample` | `gwpy` `.resample()`, 4096 → 2048 Hz | median rel. PSD diff 0.34%, p95 0.56% |
| `spectrogram` | `gwpy` `.spectrogram2()` | median rel. diff 0.35% |
| `q_transform`/`qscan_peak` | `gwpy` `q_transform`, 4 s H1 around GW150914 | peak time/frequency/Q exact; energy/SNR rel. diff ~2e-8 (`gwpy` stores these as Float32) |
| Q-tile internals | `gwpy` `QTile` | window values, indices and single-tile energy bit-for-bit |
| `load_gwf_strain` | `load_gwosc_strain` on the matching HDF5 file | bit-for-bit |
| `gmst` | `lal.GreenwichMeanSiderealTime` | bit-for-bit |
| `compute_antenna_response` | `lal.ComputeDetAMResponse` | max rel. diff ~1.2e-6 |
| `time_delay_from_geocenter` | `lal.TimeDelayFromEarthCenter` | bit-for-bit |
| `aligo_design_*`/`avirgo_design_*` | `lalsimulation.SimNoisePSD*` | max rel. diff ~3e-15 |
| `gps_to_utc` | `gwpy.time.tconvert` | exact |
