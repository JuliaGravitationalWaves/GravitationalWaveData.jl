using Test
using GWData
using Random
using Statistics: median

@testset "csd/coherence on synthetic signals" begin
    Random.seed!(42)
    dt = 1 / 4096
    sigma = 3.7e-2
    nseg = 4096            # 1s segments
    nseglist = 400          # enough segments (50% overlap) for the mean estimator to converge
    n = nseg * (nseglist + 1) ÷ 2 + nseg  # gives ~nseglist 50%-overlapping segments

    @testset "coherence == 1.0 for a signal against itself" begin
        x = TimeSeries(sigma .* randn(n), 0.0, dt; name="x")
        coh = GWData.coherence(x, x; seglen=nseg, overlap=0.5)
        @test coh isa FrequencySeries{Float64}
        @test coh.f0 == 0.0
        @test coh.df ≈ 1 / (nseg * dt)
        # Pxy(x,x) == Pxx exactly (conj(Fx)*Fx == abs2(Fx) up to floating-point roundoff in the
        # imaginary part), so |Pxy|^2/(Pxx*Pyy) == 1 to machine precision, not just statistically.
        @test all(v -> isapprox(v, 1.0; atol=1e-8), coh.data)
    end

    @testset "coherence is near 0 for independent white noise" begin
        x = TimeSeries(sigma .* randn(n), 0.0, dt; name="x")
        y = TimeSeries(sigma .* randn(n), 0.0, dt; name="y")
        coh = GWData.coherence(x, y; seglen=nseg, overlap=0.5)
        @test all(v -> 0.0 <= v <= 1.0, coh.data)
        # For ~O(nseglist) independent-ish (correlated by 50% overlap, so somewhat fewer effective
        # degrees of freedom) segments, the null-hypothesis coherence estimator has mean ~1/n_eff,
        # small compared to 1. Loose bounds -- this is a sanity check, not a precision claim.
        @test median(coh.data) < 0.1
        @test all(v -> v < 0.5, coh.data)  # no wild outliers indicating a scale/indexing bug
    end

    @testset "coherence tracks a known signal-to-noise ratio" begin
        # y = x + independent noise: analytically, for additive independent noise with
        # signal-power ratio r = var(x)/(var(x)+var(noise)) at every frequency (flat spectra
        # here), the expected coherence is ~r^2 / ... actually for this simple additive model
        # coherence = |Pxy|^2/(Pxx*Pyy) with Pxy=Pxx (since cov(x,y)=var(x)), Pyy=Pxx+Pnoise, so
        # coherence -> Pxx/(Pxx+Pnoise) = r. Use a large signal-to-noise ratio so coherence should
        # be close to 1, and confirm it's clearly above the independent-noise-floor case above.
        x = sigma .* randn(n)
        noise_small = 0.1 * sigma .* randn(n)   # var ratio: signal is 100x noise power
        xt = TimeSeries(x, 0.0, dt; name="x")
        yt = TimeSeries(x .+ noise_small, 0.0, dt; name="y")
        coh = GWData.coherence(xt, yt; seglen=nseg, overlap=0.5)
        expected = var(x) / (var(x) + var(noise_small))
        @test isapprox(median(coh.data), expected; rtol=0.1)
        @test median(coh.data) > 0.9
    end

    @testset "csd Hermitian symmetry: csd(y,x) == conj(csd(x,y))" begin
        x = TimeSeries(sigma .* randn(n), 0.0, dt)
        y = TimeSeries(sigma .* randn(n), 0.0, dt)
        pxy = csd(x, y; seglen=nseg, overlap=0.5)
        pyx = csd(y, x; seglen=nseg, overlap=0.5)
        @test pxy.data ≈ conj.(pyx.data)
    end

    @testset "csd(x,x) matches psd(x; method=:mean) (real part; imaginary part ~0)" begin
        x = TimeSeries(sigma .* randn(n), 0.0, dt; name="x")
        pxx_csd = csd(x, x; seglen=nseg, overlap=0.5)
        pxx_psd = psd(x; seglen=nseg, overlap=0.5, method=:mean)
        @test real.(pxx_csd.data) ≈ pxx_psd.data
        @test maximum(abs.(imag.(pxx_csd.data))) < 1e-8 * maximum(abs.(real.(pxx_csd.data)))
    end

    @testset "csd/coherence return types and axes" begin
        x = TimeSeries(sigma .* randn(n), 0.0, dt)
        y = TimeSeries(sigma .* randn(n), 0.0, dt)
        c = csd(x, y; seglen=nseg, overlap=0.5)
        k = GWData.coherence(x, y; seglen=nseg, overlap=0.5)
        @test eltype(c) <: Complex
        @test eltype(k) <: Real
        @test length(c) == length(k)
        @test c.f0 == k.f0 == 0.0
        @test c.df == k.df
    end

    @testset "error handling: dt/length mismatch" begin
        x = TimeSeries(randn(1000), 0.0, dt)
        mismatched_dt = TimeSeries(randn(1000), 0.0, 2 * dt)
        mismatched_len = TimeSeries(randn(500), 0.0, dt)

        @test_throws ArgumentError csd(x, mismatched_dt; seglen=100)
        @test_throws ArgumentError csd(x, mismatched_len; seglen=100)
        @test_throws ArgumentError GWData.coherence(x, mismatched_dt; seglen=100)
        @test_throws ArgumentError GWData.coherence(x, mismatched_len; seglen=100)

        try
            csd(x, mismatched_dt; seglen=100)
            @test false
        catch e
            @test e isa ArgumentError
            @test occursin("sampling rate", e.msg)
        end

        try
            csd(x, mismatched_len; seglen=100)
            @test false
        catch e
            @test e isa ArgumentError
            @test occursin("lengths differ", e.msg)
        end
    end

    @testset "error handling: segment too short" begin
        short = TimeSeries(randn(100), 0.0, dt)
        short2 = TimeSeries(randn(100), 0.0, dt)
        @test_throws ArgumentError csd(short, short2; seglen=1000)
        @test_throws ArgumentError GWData.coherence(short, short2; seglen=1000)
    end
end

# ================================================================================================
# Cross-validated regression spot-checks on real data: these are GWData.jl's own `csd`/`coherence`
# output at specific frequency bins on the real GW150914 H1/L1 off-source segment, recorded from an
# actual run of `examples/validation/compare_coherence.jl` that confirmed close agreement with
# gwpy's `TimeSeries.csd()`/`.coherence()` on the same segment/parameters:
#
#   CSD magnitude relative difference (20-1024 Hz band, 4017 points): median 9.26e-5, max 7.0e-3
#   Coherence absolute difference (same band):                        median 2.2e-7, max 2.3e-5
#
# (better than this package's own PSD.jl/gwpy baseline of ~0.14%/1.14% median/max -- see
# `compare_psd.jl` -- because coherence's numerator/denominator normalization mostly cancels).
# The values below are GWData.jl's own output at 8 spot-check frequencies from that run (already
# confirmed gwpy-consistent to the precision above), baked in as an exact regression so a future
# change to `csd`/`coherence` that silently alters the numerics fails this Julia-only test, without
# needing gwpy available at test time.
# ================================================================================================

@testset "csd/coherence cross-validated regression (real GW150914 H1/L1)" begin
    h1_path = find_test_data(GW150914_H1_HDF5)
    l1_path = find_test_data(GW150914_L1_HDF5)

    if h1_path === nothing || l1_path === nothing
        @test_skip "GW150914 H1/L1 HDF5 files not found (see test/testdata.jl)"
    else
        h1_strain, _ = load_gwosc_strain(h1_path)
        l1_strain, _ = load_gwosc_strain(l1_path)

        trigger_time = 1126259462.4
        duration = 4.0
        post_trigger = 2.0
        csd_segment = 1024.0
        analysis_start = trigger_time - (duration - post_trigger)
        off_source_end = analysis_start
        off_source_start = off_source_end - csd_segment

        h1_off = crop(h1_strain, off_source_start, off_source_end)
        l1_off = crop(l1_strain, off_source_start, off_source_end)
        nseg = round(Int, duration / h1_strain.dt)

        csd_result = csd(h1_off, l1_off; seglen=nseg, overlap=0.5)
        coh_result = GWData.coherence(h1_off, l1_off; seglen=nseg, overlap=0.5)
        freqs = collect(frequencies(csd_result))

        @test length(csd_result) == length(coh_result)
        @test all(v -> 0.0 <= v <= 1.0 + 1e-9, coh_result.data)
        @test csd_result.df ≈ 0.25   # 1/(4.0s * dt), dt = 1/4096s

        # Spot checks: (freq_hz, csd_re, csd_im, coherence), taken from an actual run of
        # compare_coherence.jl (see header comment above) -- GWData.jl's own output, already
        # confirmed gwpy-consistent at these exact points.
        spot_checks = [
            (20.0, -1.4022590884947186e-45, 1.1463018036657848e-45, 0.0014324602002561884),
            (50.0, -1.1702380592147696e-47, 7.8813172305181e-49, 0.002781040172850796),
            (100.0, -4.1448508130090374e-48, -3.9056604327325714e-48, 0.004018980333789443),
            (150.0, -5.73652293953793e-49, -3.2126534478615254e-49, 8.868235053076232e-5),
            (200.0, 9.93300097750577e-49, 2.9453909280468838e-49, 0.0002480128771854058),
            (300.0, 2.700742032078179e-47, -2.0916227778486978e-47, 0.0010920437097106644),
            (500.0, -7.82851055957936e-46, -7.988810403256803e-46, 0.00149644120736084),
            (1000.0, 4.784614425967114e-48, 1.259465730724406e-47, 0.0008257166716756277),
        ]
        for (f, re_expected, im_expected, coh_expected) in spot_checks
            idx = argmin(abs.(freqs .- f))
            @test freqs[idx] ≈ f
            @test real(csd_result.data[idx]) ≈ re_expected rtol=1e-6
            @test imag(csd_result.data[idx]) ≈ im_expected rtol=1e-6
            @test coh_result.data[idx] ≈ coh_expected rtol=1e-6
        end
    end
end
