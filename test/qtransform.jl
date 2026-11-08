using Random
using GWData.QTransform: normalized_fft, QTileParams, qtile_window, qtile_data_indices,
                          qtile_padding, qtile_energy, q_values, tiling_frange,
                          plane_frange, plane_frequencies, q_plane_transform, qscan_peak,
                          _ntiles, _windowsize, qprime

# ================================================================================================
# Algebraic checks: q_values / tiling_frange / plane_frequencies / QTileParams derived quantities
# are direct, deterministic ports of gwpy's formulas (see src/GWData/QTransform.jl). These
# reference numbers were computed by this package's own functions and are pinned here as
# regression values; the numbers most at risk of a porting mistake (the tile window, its data
# indices, and a full tile's IFFT energy) are separately cross-checked bit-for-bit against gwpy's
# own `QTile` in the "gwpy cross-validation" testset below, which is the check that actually
# resolves the algorithm's ambiguity (see that testset for how/why).
# ================================================================================================

@testset "q_values" begin
    # gwpy DEFAULT_QRANGE=(4,64), DEFAULT_MISMATCH=0.2, duration=4.0s -- QTiling._iter_qs formula.
    qs = q_values((4.0, 64.0), 4.0, 0.2)
    @test length(qs) == 4
    @test isapprox(qs, [5.65685424949238, 11.31370849898476, 22.627416997969522,
                        45.25483399593904]; rtol=1e-12)
    # monotonically increasing, geometric-ish spacing
    @test issorted(qs)
    @test all(qs .> 0)
end

@testset "tiling_frange" begin
    fmin, fmax = tiling_frange((4.0, 64.0), (0.0, Inf), 4.0, 4096.0, 0.2)
    @test isapprox(fmin, 90.0316316157106; rtol=1e-10)
    @test isapprox(fmax, 1291.0530521679268; rtol=1e-10)
    # an explicit, already-valid frange passes through untouched
    @test tiling_frange((4.0, 64.0), (20.0, 500.0), 4.0, 4096.0, 0.2) == (20.0, 500.0)
    # fmax clamped with a warning if above what the Q range supports
    @test (@test_logs (:warn,) match_mode=:any tiling_frange((4.0, 64.0), (20.0, 1e6), 4.0, 4096.0,
                                                              0.2)) == (20.0, 1291.0530521679268)
end

@testset "plane_frange (standalone)" begin
    # QPlane's own frange resolution matches QTiling's for a single Q (used when a plane is
    # constructed directly with an unresolved (0, Inf) range, not through a multi-Q q_scan).
    fmin, fmax = plane_frange(45.25483399593904, 4.0, (0.0, Inf), 4096.0)
    @test isapprox(fmin, 90.0316316157106; rtol=1e-10)
end

@testset "plane_frequencies" begin
    freqs = plane_frequencies(8.0, (20.0, 500.0), 4.0, 0.2)
    @test length(freqs) == 26
    @test isapprox(freqs[1], 21.25; rtol=1e-12)
    @test isapprox(freqs[end], 469.75; rtol=1e-12)
    @test issorted(freqs)
    # every frequency is a multiple of 1/duration = 0.25 Hz (the tiling snaps to the FFT bin grid)
    @test all(f -> isapprox(f / 0.25, round(f / 0.25); atol=1e-9), freqs)
end

@testset "QTileParams derived quantities" begin
    t = QTileParams(8.0, 100.0, 4.0, 4096.0, 0.2)
    @test qprime(t) ≈ 8.0 / sqrt(11)
    @test _ntiles(t) == 1024
    @test _windowsize(t) == 331
    @test qtile_padding(t) == (346, 347)
    idx = qtile_data_indices(t)
    @test length(idx) == 331
    @test issorted(idx)
    win = qtile_window(t)
    @test length(win) == 331
    # bi-square window: zero at the edges, maximal (and positive) at the center
    @test win[1] > 0 && win[end] > 0  # edges are near but not exactly zero (offset != +-windowsize/2)
    @test win[(end + 1) ÷ 2] == maximum(win)
    @test all(win .>= 0)  # (1-x^2)^2 * positive-norm is never negative
end

@testset "normalized_fft matches gwpy's TimeSeries.fft() convention" begin
    # A pure sinusoid at an exact FFT bin frequency should have (twice the analytic amplitude
    # /2, i.e. the amplitude itself) at that bin and ~0 elsewhere, per gwpy's `dft[1:] *= 2`
    # one-sided-amplitude convention.
    fs = 1024.0
    dur = 2.0
    N = Int(fs * dur)
    dt = 1 / fs
    t = (0:(N - 1)) .* dt
    f0 = 50.0  # exact bin: f0*dur = 100, integer
    A = 3.0
    data = A .* sin.(2pi * f0 .* t)
    ts = TimeSeries(data, 0.0, dt)
    fseries = normalized_fft(ts)
    @test fseries.df ≈ 1 / dur
    bin = round(Int, f0 / fseries.df) + 1
    # sin has no cosine component -> the real part is ~0, but |X| at that bin should be ~A
    @test isapprox(abs(fseries.data[bin]), A; atol=1e-9)
    # everywhere else (away from the bin) should be ~0 for this exact-bin-frequency sinusoid
    other = abs.(fseries.data[[1:(bin - 2); (bin + 2):end]])
    @test maximum(other) < 1e-9
end

@testset "qtile_energy: raw energy shape and scale sanity" begin
    # Deterministic sum-of-sinusoids test signal (same construction as the cross-validation
    # reference case), all components at exact FFT-bin frequencies. A tile centered on 100 Hz (no
    # signal there) should see only leakage-floor energy -- this is actually a stringent indexing
    # check on its own: any off-by-one in `qtile_data_indices` breaks the near-exact cancellation
    # and inflates the energy by many orders of magnitude.
    fs = 4096.0
    dur = 4.0
    N = Int(fs * dur)
    dt = 1 / fs
    t = (0:(N - 1)) .* dt
    data = sin.(2pi * 35 .* t) .+ 0.5 .* sin.(2pi * 150 .* t .+ 0.3) .+ 0.1 .* cos.(2pi * 7 .* t)
    ts = TimeSeries(data, 0.0, dt)
    fseries = normalized_fft(ts)

    off_resonance = QTileParams(8.0, 100.0, dur, fs, 0.2)
    e_off, dt_off = qtile_energy(fseries.data, off_resonance)
    @test dt_off ≈ dur / 1024
    @test length(e_off) == 1024
    @test maximum(e_off) < 1e-30  # leakage floor, not a real signal

    on_resonance = QTileParams(8.0, 150.0, dur, fs, 0.2)
    e_on, _ = qtile_energy(fseries.data, on_resonance)
    # A stationary sinusoid excites a single frequency-domain bin, so its IFFT is a
    # constant-*magnitude* complex exponential (only phase rotates) -- not a localized burst --
    # scaled by 1/ntiles (the 150 Hz component lands one bin away from this tile's center due to
    # the same `+1` index-offset quirk noted on `qtile_data_indices`, i.e. at offset -1 not 0).
    # gwpy cross-checked (via examples/validation/qtile_reference_step.jl, same signal.csv,
    # `QTile(8.0, 150.0, 4.0, 4096.0, mismatch=0.2).transform(fseries, norm=False)`): gwpy's energy
    # max was 3.685318469980373e-11 -- bit-for-bit identical (to every digit printed) to this
    # package's own value below, the strongest confirmation obtained for this task's biggest risk
    # (the tiling/indexing convention, see src/GWData/QTransform.jl's module docstring).
    @test all(isapprox.(e_on, e_on[1]; rtol=1e-9))  # flat (constant-magnitude), as expected
    @test maximum(e_on) > 1e6 * maximum(e_off)
    @test isapprox(maximum(e_on), 3.685318469980373e-11; rtol=1e-12)
end

@testset "end-to-end: synthetic sine-Gaussian burst recovered by q_transform" begin
    # Self-contained (no gwpy needed) sanity check that the full q_values -> plane ->
    # tile -> peak pipeline recovers a known injected burst's time/frequency, mirroring gwpy's own
    # q_transform docstring example (Gaussian-enveloped sinusoid on top of noise).
    Random.seed!(1234)
    fs = 2048.0
    dur = 4.0
    N = Int(fs * dur)
    dt = 1 / fs
    t = (0:(N - 1)) .* dt
    t_glitch = 2.0
    f_glitch = 200.0
    sigma = 0.02
    burst = 4.0 .* exp.(-((t .- t_glitch) .^ 2) ./ (2 * sigma^2)) .* sin.(2pi * f_glitch .* (t .- t_glitch))
    noise = randn(N)
    ts = TimeSeries(noise .+ burst, 0.0, dt)

    scan, peak = q_transform(ts; frange=(20.0, 500.0))
    @test isapprox(peak.time, t_glitch; atol=0.01)
    @test isapprox(peak.frequency, f_glitch; atol=5.0)
    @test peak.snr > 10  # loud injection, should stand out clearly above the noise floor
    @test scan.q == peak.q
end

@testset "qscan_peak search window" begin
    rows = [
        GWData.QTransform.QRow(50.0, [1.0, 2.0, 100.0, 3.0], 0.0, 1.0),
        GWData.QTransform.QRow(60.0, [1.0, 50.0, 2.0, 3.0], 0.0, 1.0),
    ]
    scan = GWData.QTransform.QScan(8.0, (20.0, 500.0), rows)
    peak = qscan_peak(scan)
    @test peak.energy == 100.0
    @test peak.frequency == 50.0
    @test peak.time == 2.0  # index 3 (1-based) -> t = 0 + (3-1)*1.0
    # restricting the search window excludes the global peak, so the row-2 peak wins instead
    peak2 = qscan_peak(scan; search=(0.5, 1.5))
    @test peak2.energy == 50.0
    @test peak2.frequency == 60.0
end

# ================================================================================================
# gwpy cross-validation (intermediate quantities): the exact numbers below were generated by
# `examples/validation/qtile_reference_step.jl` (run against the same
# `LALSuite.jl/examples/psd` CondaPkg gwpy 4.0.1 environment the other GWData.jl comparison
# scripts use) for the SAME deterministic signal and tile parameters
# (q=8.0, frequency=100.0 Hz, duration=4.0s, sampling=4096.0 Hz, mismatch=0.2) as the
# "qtile_energy: raw energy shape and scale sanity" testset above. This is the check that actually
# resolved this task's biggest risk (the `+1` index-offset quirk in `qtile_data_indices` and the
# exact bi-square window normalization in `qtile_window` -- see src/GWData/QTransform.jl's module
# docstring): it compares GWData.jl's Julia output bit-for-bit (to floating-point tolerance)
# against gwpy's own `QTile.get_window()`/`get_data_indices()`/`.transform(..., norm=False)`,
# not just gwpy's final rendered Q-scan image.
# ================================================================================================

@testset "gwpy cross-validation: single QTile intermediate quantities (q=8, f=100 Hz)" begin
    # Reference values from gwpy 4.0.1's `QTile(8.0, 100.0, 4.0, 4096.0, mismatch=0.2)` (via
    # `examples/validation/qtile_reference_step.jl`, run against the same signal.csv
    # this file's "qtile_energy: raw energy shape and scale sanity" testset generates).
    t = QTileParams(8.0, 100.0, 4.0, 4096.0, 0.2)

    # Exact-integer quantities: bit-for-bit match.
    @test _ntiles(t) == 1024
    @test _windowsize(t) == 331
    @test qtile_padding(t) == (346, 347)
    @test isapprox(qprime(t), 2.412090756622109; rtol=1e-12)

    # Data indices: gwpy's `get_data_indices()` (0-based) was [236, 237, 238, 239, 240, ...] for
    # the first 5 entries -- this package's `qtile_data_indices` is 1-based, so it should be
    # exactly one more at every position (the `+1` conversion, on top of gwpy's own already-ported
    # `+1` quirk -- see qtile_data_indices's docstring).
    gwpy_indices0based_first5 = [236, 237, 238, 239, 240]
    @test qtile_data_indices(t)[1:5] == gwpy_indices0based_first5 .+ 1

    # Window values: gwpy's `get_window()` first 5 entries.
    gwpy_window_first5 = [1.5227448228310888e-6, 7.3457411608057715e-6, 1.7452639013314443e-5,
                          3.1764426801053256e-5, 5.020257619596958e-5]
    @test isapprox(qtile_window(t)[1:5], gwpy_window_first5; rtol=1e-9)

    # Raw tile energy at this off-resonance frequency: both implementations agree the result is
    # at the ~1e-38 floating-point leakage floor (19+ orders of magnitude below any real signal),
    # though the exact residual isn't expected to match bit-for-bit (different FFT backends --
    # FFTW vs numpy's pocketfft -- round differently at the level of pure floating-point noise).
    # gwpy's own `.transform(fseries, norm=False)` energy max was 3.596288596541836e-38 at
    # (0-based) index 556; this package's is 3.582559098628859e-38 at (1-based) index 557 -- same
    # index (556 + 1 = 557) and same order of magnitude.
    fs = 4096.0
    dur = 4.0
    N = Int(fs * dur)
    dt = 1 / fs
    tt = (0:(N - 1)) .* dt
    data = sin.(2pi * 35 .* tt) .+ 0.5 .* sin.(2pi * 150 .* tt .+ 0.3) .+ 0.1 .* cos.(2pi * 7 .* tt)
    fseries = normalized_fft(TimeSeries(data, 0.0, dt))
    energy, _ = qtile_energy(fseries.data, t)
    @test argmax(energy) == 557
    @test maximum(energy) < 1e-30
    @test isapprox(maximum(energy), 3.596288596541836e-38; atol=1e-38)  # same floating-point-noise order of magnitude
end

@testset "gwpy cross-validation: single QTile intermediate quantities (q=8, f=150 Hz, on-resonance)" begin
    # Reference values from gwpy 4.0.1's `QTile(8.0, 150.0, 4.0, 4096.0, mismatch=0.2)` (same
    # setup as the f=100 Hz testset above). This is the on-resonance companion case (the
    # signal has a real 0.5-amplitude component at 150 Hz) -- see the "qtile_energy: raw energy
    # shape and scale sanity" testset above for the bit-for-bit energy-value match
    # (3.685318469980373e-11, identical to every digit gwpy printed); this testset covers the
    # remaining tile-parameter/window/index quantities for that same case.
    t = QTileParams(8.0, 150.0, 4.0, 4096.0, 0.2)
    @test _ntiles(t) == 1024
    @test _windowsize(t) == 497
    @test qtile_padding(t) == (263, 264)
    @test isapprox(qprime(t), 2.412090756622109; rtol=1e-12)

    gwpy_indices0based_first5 = [353, 354, 355, 356, 357]
    @test qtile_data_indices(t)[1:5] == gwpy_indices0based_first5 .+ 1

    gwpy_window_first5 = [4.469911220473276e-7, 2.43549286859392e-6, 5.997772542177978e-6,
                          1.1114617847930875e-5, 1.776689443132793e-5]
    @test isapprox(qtile_window(t)[1:5], gwpy_window_first5; rtol=1e-9)
end

# ================================================================================================
# gwpy cross-validation (full end-to-end, real GW150914 data): the intermediate-quantity checks
# above confirm the tiling/window/indexing port bit-for-bit on a synthetic signal; this checks the
# *whole* q_values -> multi-Q search -> peak read-out pipeline on real interferometer data.
# Reference values from `examples/validation/compare_qtransform.jl`, run on real H1 GW150914
# strain (`H-H1_LOSC_4_V1-1126256640-4096.hdf5`, 4s segment straddling the trigger,
# frange=(20,500) Hz, qrange=(4,64), mismatch=0.2, no whitening) against gwpy 4.0.2's `q_scan`:
#   time:      GWData=1.1262594603999023e9  gwpy=1.1262594603999023e9  abs_diff=0.0 s (exact)
#   frequency: GWData=140.25                gwpy=140.25                (exact)
#   q:         GWData=5.65685424949238      gwpy=5.65685424949238      (exact)
#   energy:    GWData=1.5721025376480368e6  gwpy=1.5721025e6           rel_diff=2.39e-8
#   snr:       GWData=1773.1906483218531    gwpy=1773.1906738281       rel_diff=1.44e-8
# The energy/snr differences are NOT algorithmic disagreement: gwpy's `QGram.peak` stores those
# two fields as numpy `float32` (confirmed directly: `peak['energy'].dtype == float32`), while
# GWData.jl computes in Float64 throughout, so the observed ~1e-8 relative differences are exactly
# at float32's ~1.19e-7 machine epsilon -- i.e. this is a bit-for-bit match modulo gwpy's own
# float32 storage, at least as tight as this codebase's established intermediate-quantity bar, and
# far tighter than its own precedent for full end-to-end real-data comparisons (spectrogram/
# resample: ~0.3-0.6% median). time/frequency/q match exactly because gwpy computes those from
# (Float64) index/parameter arithmetic rather than the float32 energy array.
@testset "gwpy cross-validation: full end-to-end GW150914 q_transform (real H1 data)" begin
    strain_path = find_test_data(GW150914_H1_HDF5)

    if strain_path === nothing
        @test_skip "GW150914 H1 HDF5 file not found (see test/testdata.jl)"
    else
        strain, _dq = load_gwosc_strain(strain_path)
        trigger_time = 1126259462.4
        segment = crop(strain, trigger_time - 2.0, trigger_time + 2.0)
        _scan, peak = q_transform(segment; frange=(20.0, 500.0))

        @test peak.time == 1.1262594603999023e9
        @test peak.frequency == 140.25
        @test peak.q == 5.65685424949238
        # float32-epsilon-level agreement with gwpy's q_scan peak (see comment above the testset).
        @test isapprox(peak.energy, 1.5721025e6; rtol=1e-6)
        @test isapprox(peak.snr, 1773.1906738281; rtol=1e-6)
    end
end
