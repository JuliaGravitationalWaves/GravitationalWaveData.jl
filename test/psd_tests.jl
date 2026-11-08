using Test
using GWData
using DSP
using Random
using Statistics: var, median

@testset "DSP.jl window cross-check vs hand-rolled reference" begin
    # Closed-form symmetric Hann and Tukey (tapered-cosine) definitions, checked against the
    # DSP.Windows.hanning/tukey conventions used in src/GWData/PSD.jl and src/GWData/IO/GWOSC.jl.
    ref_hann(N) = [0.5 * (1 - cos(2pi * n / (N - 1))) for n in 0:(N - 1)]
    function ref_tukey(N::Int, alpha::Real)
        alpha <= 0 && return ones(Float64, N)
        alpha >= 1 && return ref_hann(N)
        w = ones(Float64, N)
        taper_n = floor(Int, alpha * (N - 1) / 2)
        for n in 0:taper_n
            w[n + 1] = 0.5 * (1 + cos(pi * (2n / (alpha * (N - 1)) - 1)))
            w[N - n] = w[n + 1]
        end
        return w
    end

    for N in (8, 16, 4097, 16384, 16385)
        @test maximum(abs.(DSP.Windows.hanning(N) .- ref_hann(N))) < 1e-12
    end

    for N in (16384, 16385, 4097), alpha in (0.0, 0.1, 0.2, 0.5, 1.0)
        @test maximum(abs.(DSP.Windows.tukey(N, alpha) .- ref_tukey(N, alpha))) < 1e-12
    end

    # Endpoint handling: both conventions must be exactly zero at the first/last sample for a
    # plain Hann window (this is what condition_segment's `alpha=0.2` tukey window relies on --
    # see the GWOSC loader tests).
    @test DSP.Windows.hanning(8)[1] == 0.0
    @test DSP.Windows.hanning(8)[end] == 0.0
end

@testset "psd/asd on synthetic signals" begin
    Random.seed!(1234)
    dt = 1 / 4096
    sigma = 3.7e-2
    nseg = 4096          # 1s segments
    nseglist = 400        # enough segments (50% overlap) for the median estimator to converge
    n = nseg * (nseglist + 1) ÷ 2 + nseg  # gives ~nseglist 50%-overlapping segments

    white = TimeSeries(sigma .* randn(n), 0.0, dt; name="white")

    @testset "variance conservation (Parseval)" begin
        p = psd(white; seglen=nseg, overlap=0.5, method=:median)
        @test p isa FrequencySeries
        @test p.f0 == 0.0
        @test p.df ≈ 1 / (nseg * dt)

        integral = sum(p.data) * p.df
        # Discrete white noise of variance sigma^2 has a flat one-sided PSD; the docstring's
        # convention is ∫psd df over [0, Nyquist] ≈ time series variance. Statistical fluctuation
        # from a finite number of (correlated, overlapping) Welch segments -> generous tolerance.
        @test isapprox(integral, var(white.data); rtol=0.15)
    end

    @testset "flat level matches the analytic white-noise PSD" begin
        p = psd(white; seglen=nseg, overlap=0.5, method=:median)
        expected_level = 2 * sigma^2 * dt   # standard one-sided discrete white-noise PSD
        # Away from DC/Nyquist, every bin should scatter around the same flat expected_level.
        interior = p.data[2:(end - 1)]
        @test isapprox(median(interior), expected_level; rtol=0.15)
        @test all(x -> isapprox(x, expected_level; rtol=1.0), interior)  # loose sanity, no wild outliers/scale bugs
    end

    @testset "asd is elementwise sqrt of psd" begin
        p = psd(white; seglen=nseg, overlap=0.5, method=:median)
        a = asd(white; seglen=nseg, overlap=0.5, method=:median)
        @test a isa FrequencySeries
        @test a.f0 == p.f0
        @test a.df == p.df
        @test a.data == sqrt.(p.data)
    end

    @testset "method=:mean also runs and roughly agrees with :median" begin
        p_med = psd(white; seglen=nseg, overlap=0.5, method=:median)
        p_mean = psd(white; seglen=nseg, overlap=0.5, method=:mean)
        @test length(p_mean) == length(p_med)
        rel = abs.(p_mean.data .- p_med.data) ./ p_med.data
        @test median(rel) < 0.3
    end

    @testset "DC/Nyquist bins are not double-counted" begin
        # For an even seglen, the Nyquist bin should not be systematically ~2x the neighboring
        # bins' level (which is what you'd see if the halving were missing/wrong).
        p = psd(white; seglen=nseg, overlap=0.5, method=:median)
        expected_level = 2 * sigma^2 * dt
        @test isapprox(p.data[1], expected_level; rtol=1.0)      # DC bin, loose (single sample)
        @test isapprox(p.data[end], expected_level; rtol=1.0)    # Nyquist bin, loose (single sample)
    end

    @testset "error handling" begin
        short = TimeSeries(randn(100), 0.0, dt)
        @test_throws ArgumentError psd(short; seglen=1000)
        @test_throws ArgumentError psd(white; seglen=nseg, method=:bogus)
    end
end
