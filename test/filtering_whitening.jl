using Test
using GravitationalWaveData
using DSP: filt, Butterworth, Lowpass, digitalfilter
using FFTW: rfft
using Random

# ================================================================================================
# Unit tests: bandpass/highpass/lowpass/notch filters (a known sine wave should pass or be
# attenuated as expected) and whiten (a synthetic colored-noise realization should come out with an
# approximately flat spectrum).
# ================================================================================================

"""
    _rms(x) -> Float64

Root-mean-square, used throughout to measure how much a tone survives filtering, robust to
transient edge effects that a raw peak-amplitude check would be sensitive to.
"""
_rms(x) = sqrt(sum(abs2, x) / length(x))

"""
    _periodogram_asd(x, dt) -> Vector{Float64}

A minimal (non-Welch-averaged) one-sided periodogram ASD estimate, `sqrt(2*dt/N) * |rfft(x)|`
(standard one-sided PSD density normalization, matching the `whiten` docstring's convention),
used only here to (a) build a self-consistent ASD to feed `whiten` for the flatness test,
independently of `psd`/`asd`, and (b) measure how flat the whitened output's own spectrum is
afterwards.
"""
function _periodogram_asd(x::AbstractVector{<:Real}, dt::Real)
    N = length(x)
    xd = x .- (sum(x) / N)
    X = rfft(xd)
    return sqrt.(2 * dt / N) .* abs.(X)
end

@testset "Filtering & whitening" begin
    @testset "bandpass passes in-band, attenuates out-of-band tones" begin
        fs = 512.0
        dt = 1 / fs
        t = collect(0:dt:(10 - dt))
        f_in, f_low_out, f_high_out = 60.0, 5.0, 200.0
        sig = sin.(2pi * f_in .* t) .+ sin.(2pi * f_low_out .* t) .+ sin.(2pi * f_high_out .* t)
        ts = TimeSeries(sig, 0.0, dt)

        bp = bandpass(ts, 40.0, 80.0; order=4)
        @test length(bp) == length(ts)
        # Well inside the band's steady-state region (skip filter settle-in near the edges): the
        # in-band 60 Hz tone alone has RMS 1/sqrt(2); the two out-of-band tones should be almost
        # entirely gone from that same window, so total RMS should sit close to the single-tone
        # value, not the three-tone value (~1.22).
        core = bp.data[Int(fs * 2):(end - Int(fs * 2))]
        @test isapprox(_rms(core), 1 / sqrt(2); rtol=0.1)

        # Cross-check against a pure single 60 Hz tone processed identically: bandpass output
        # should track it closely in the same interior window.
        pure_in_band = sin.(2pi * f_in .* t)
        pure_ts = TimeSeries(pure_in_band, 0.0, dt)
        bp_pure = bandpass(pure_ts, 40.0, 80.0; order=4)
        @test isapprox(bp.data[(end ÷ 2)], bp_pure.data[(end ÷ 2)]; atol=0.05)
    end

    @testset "highpass attenuates low frequencies, passes high" begin
        fs = 512.0
        dt = 1 / fs
        t = collect(0:dt:(10 - dt))
        low = sin.(2pi * 5.0 .* t)
        high = sin.(2pi * 150.0 .* t)
        ts = TimeSeries(low .+ high, 0.0, dt)

        hp = highpass(ts, 50.0; order=4)
        core = hp.data[(Int(fs) * 2):(end - Int(fs) * 2)]
        @test isapprox(_rms(core), 1 / sqrt(2); rtol=0.15)

        low_only = highpass(TimeSeries(low, 0.0, dt), 50.0; order=4)
        core_low = low_only.data[(Int(fs) * 2):(end - Int(fs) * 2)]
        @test _rms(core_low) < 0.05   # low tone should be almost entirely gone
    end

    @testset "lowpass passes low frequencies, attenuates high" begin
        fs = 512.0
        dt = 1 / fs
        t = collect(0:dt:(10 - dt))
        low = sin.(2pi * 5.0 .* t)
        high = sin.(2pi * 150.0 .* t)
        ts = TimeSeries(low .+ high, 0.0, dt)

        lp = lowpass(ts, 30.0; order=4)
        core = lp.data[(Int(fs) * 2):(end - Int(fs) * 2)]
        @test isapprox(_rms(core), 1 / sqrt(2); rtol=0.15)

        high_only = lowpass(TimeSeries(high, 0.0, dt), 30.0; order=4)
        core_high = high_only.data[(Int(fs) * 2):(end - Int(fs) * 2)]
        @test _rms(core_high) < 0.05   # high tone should be almost entirely gone
    end

    @testset "notch removes a narrow band, leaves the rest" begin
        fs = 512.0
        dt = 1 / fs
        t = collect(0:dt:(10 - dt))
        hum = sin.(2pi * 60.0 .* t)
        signal = sin.(2pi * 100.0 .* t)
        ts = TimeSeries(hum .+ signal, 0.0, dt)

        notched = notch(ts, 60.0; bandwidth=2.0, order=4)
        core = notched.data[(Int(fs) * 2):(end - Int(fs) * 2)]
        # 60 Hz should be gone, 100 Hz (well outside the notch) should survive close to full RMS
        @test isapprox(_rms(core), 1 / sqrt(2); rtol=0.1)

        hum_only = notch(TimeSeries(hum, 0.0, dt), 60.0; bandwidth=2.0, order=4)
        core_hum = hum_only.data[(Int(fs) * 2):(end - Int(fs) * 2)]
        @test _rms(core_hum) < 0.05
    end

    @testset "whiten flattens a synthetic colored-noise realization" begin
        Random.seed!(42)
        fs = 1024.0
        dt = 1 / fs
        N = 2^15   # 32 s -- long enough for the FIR filter's settle-in to be a small fraction
        white = randn(N)
        # Color it with a steep lowpass Butterworth (causal `filt`, not `filtfilt`, so this is a
        # straightforward IIR-colored-noise generator, not itself zero-phase) -- gives a strongly
        # frequency-dependent (falling) spectrum, a meaningfully non-flat test case.
        color_filt = digitalfilter(Lowpass(30.0), Butterworth(4); fs=fs)
        colored = filt(color_filt, white)
        ts = TimeSeries(colored, 0.0, dt)

        asd_vals = _periodogram_asd(colored, dt)
        asd = FrequencySeries(asd_vals, 0.0, 1 / duration(ts))

        whitened = whiten(ts, asd; fduration=2.0)
        @test length(whitened) == length(ts)

        # Compare spectral flatness (stddev/mean of the ASD magnitude across the band actually
        # covered by the coloring filter's passband/rolloff, 5-200 Hz) before vs after whitening --
        # the colored signal's ASD should vary a lot there (steep rolloff around 30 Hz), the
        # whitened one much less so.
        colored_asd = _periodogram_asd(colored, dt)
        whitened_asd = _periodogram_asd(whitened.data, dt)
        freqs = collect(0:(length(colored_asd) - 1)) .* (1 / duration(ts))
        band = findall(f -> 5.0 <= f <= 200.0, freqs)

        cv(x) = (sqrt(sum(abs2, x .- (sum(x) / length(x))) / length(x))) / (sum(x) / length(x))
        cv_colored = cv(colored_asd[band])
        cv_whitened = cv(whitened_asd[band])
        @test cv_whitened < cv_colored / 2

        # Whitened output should also have ~unit variance in its non-edge-corrupted interior
        # (gwpy's own documented normalization target for stationary Gaussian input).
        interior = whitened.data[(Int(fs) * 2):(end - Int(fs) * 2)]
        @test isapprox(sum(abs2, interior) / length(interior), 1.0; rtol=0.3)
    end

    @testset "overlap-save convolve chunked path matches the single-shot path" begin
        Random.seed!(7)
        fir = randn(64)   # even length, as fir_from_transfer always produces
        # Force n > 8*length(fir) so the chunked overlap-save branch executes, and compare its
        # output to the single-shot `_convolve_same` computed over the *whole* array with the same
        # edge-window pre-tapering -- overlap-save is an exact algorithm (not an approximation), so
        # these must agree everywhere, not just approximately.
        n = 8 * length(fir) * 4 + 37   # comfortably > 8*length(fir), and not a round multiple
        data = randn(n)

        chunked = GravitationalWaveData.Filtering._overlap_save_convolve(data, fir, :hann)

        pad = cld(length(fir), 2)
        w = GravitationalWaveData.Filtering._get_window(:hann, length(fir))
        tapered = copy(data)
        tapered[1:pad] .*= w[1:pad]
        tapered[(end - pad + 1):end] .*= w[(end - pad + 1):end]
        direct = GravitationalWaveData.Filtering._convolve_same(tapered, fir)

        @test chunked ≈ direct atol=1e-8
    end
end
