using Test
using GWData
using Aqua
using Statistics

include("testdata.jl")

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(GWData)
end

@testset "GWData.jl" begin
    @testset "TimeSeries construction" begin
        ts = TimeSeries([1.0, 2.0, 3.0, 4.0], 100.0, 0.5)
        @test length(ts) == 4
        @test ts.t0 == 100.0
        @test ts.dt == 0.5
        @test ts.name === nothing
        @test eltype(ts) == Float64

        ts2 = TimeSeries([1.0, 2.0, 3.0, 4.0], 100.0; sample_rate=2.0)
        @test ts2.dt == 0.5
        @test ts == ts2

        named = TimeSeries([1, 2, 3], 0.0, 1.0; name="H1:STRAIN")
        @test named.name == "H1:STRAIN"
        @test eltype(named) == Int

        @test_throws ArgumentError TimeSeries([1.0, 2.0], 0.0, 0.0)
        @test_throws ArgumentError TimeSeries([1.0, 2.0], 0.0, -1.0)
        @test_throws ArgumentError TimeSeries([1.0, 2.0], 0.0; sample_rate=0.0)
        @test_throws ArgumentError TimeSeries([1.0, 2.0], Inf, 1.0)
    end

    @testset "TimeSeries duration/sample_rate/times" begin
        ts = TimeSeries(collect(1.0:10.0), 100.0, 0.1)
        @test duration(ts) == 1.0
        @test sample_rate(ts) == 10.0
        @test collect(times(ts)) == collect(100.0:0.1:100.9)
        @test times(ts)[1] == 100.0
        @test times(ts)[end] == 100.9
    end

    @testset "TimeSeries indexing and crop" begin
        ts = TimeSeries(collect(1.0:10.0), 100.0, 0.1; name="x")
        @test ts[1] == 1.0
        @test ts[end] == 10.0

        sub = ts[3:5]
        @test sub isa TimeSeries
        @test sub.data == [3.0, 4.0, 5.0]
        @test sub.t0 == 100.2
        @test sub.dt == 0.1
        @test sub.name == "x"

        c = crop(ts, 100.2, 100.5)   # half-open [100.2, 100.5): samples at 100.2, 100.3, 100.4
        @test c.data == [3.0, 4.0, 5.0]
        @test c.t0 == 100.2

        # open-ended crop
        @test crop(ts, nothing, 100.3).data == [1.0, 2.0, 3.0]
        @test crop(ts, 100.7, nothing).data == [8.0, 9.0, 10.0]
        @test crop(ts) == ts

        @test_throws ArgumentError crop(ts, 100.5, 100.5)   # empty range
        @test_throws ArgumentError crop(ts, 99.0, 100.5)    # before start
        @test_throws ArgumentError crop(ts, 100.5, 200.0)   # past end
    end

    @testset "TimeSeries arithmetic" begin
        a = TimeSeries([1.0, 2.0, 3.0], 0.0, 1.0)
        b = TimeSeries([10.0, 20.0, 30.0], 0.0, 1.0)

        @test (a + b).data == [11.0, 22.0, 33.0]
        @test (b - a).data == [9.0, 18.0, 27.0]
        @test (a * b).data == [10.0, 40.0, 90.0]
        @test (b / a).data == [10.0, 10.0, 10.0]
        @test (a + b).t0 == 0.0
        @test (a + b).dt == 1.0

        @test (2 * a).data == [2.0, 4.0, 6.0]
        @test (a * 2).data == [2.0, 4.0, 6.0]
        @test (a + 1).data == [2.0, 3.0, 4.0]
        @test (10 - a).data == [9.0, 8.0, 7.0]

        mismatched_dt = TimeSeries([1.0, 2.0, 3.0], 0.0, 2.0)
        @test_throws ArgumentError a + mismatched_dt

        mismatched_t0 = TimeSeries([1.0, 2.0, 3.0], 5.0, 1.0)
        @test_throws ArgumentError a + mismatched_t0

        mismatched_len = TimeSeries([1.0, 2.0], 0.0, 1.0)
        @test_throws ArgumentError a + mismatched_len

        try
            a + mismatched_dt
            @test false
        catch e
            @test e isa ArgumentError
            @test occursin("sampling rate", e.msg)
        end
    end

    @testset "TimeSeries show" begin
        ts = TimeSeries([1.0, 2.0], 0.0, 1.0; name="foo")
        s = sprint(show, ts)
        @test occursin("TimeSeries", s)
        @test occursin("foo", s)
        @test occursin("2 samples", s)
    end

    @testset "FrequencySeries construction" begin
        fs = FrequencySeries(ComplexF64[1, 2, 3, 4], 0.0, 0.25)
        @test length(fs) == 4
        @test fs.f0 == 0.0
        @test fs.df == 0.25
        @test eltype(fs) == ComplexF64

        fs2 = FrequencySeries(ComplexF64[1, 2, 3, 4], 0.25)
        @test fs2.f0 == 0.0
        @test fs2 == fs

        @test_throws ArgumentError FrequencySeries([1.0, 2.0], 0.0, 0.0)
        @test_throws ArgumentError FrequencySeries([1.0, 2.0], -1.0, 1.0)
    end

    @testset "FrequencySeries frequencies/indexing/crop" begin
        fs = FrequencySeries(collect(0.0:9.0), 0.0, 2.0; name="psd")
        @test collect(frequencies(fs)) == collect(0.0:2.0:18.0)

        sub = fs[2:4]
        @test sub.data == [1.0, 2.0, 3.0]
        @test sub.f0 == 2.0
        @test sub.df == 2.0

        c = crop(fs, 4.0, 10.0)
        @test c.data == [2.0, 3.0, 4.0]

        @test crop(fs, nothing, 4.0).data == [0.0, 1.0]
        @test crop(fs, 16.0, nothing).data == [8.0, 9.0]

        @test_throws ArgumentError crop(fs, 5.0, 5.0)
        @test_throws ArgumentError crop(fs, -2.0, 4.0)
        @test_throws ArgumentError crop(fs, 4.0, 1000.0)
    end

    @testset "FrequencySeries arithmetic" begin
        a = FrequencySeries([1.0, 2.0, 3.0], 0.0, 1.0)
        b = FrequencySeries([10.0, 20.0, 30.0], 0.0, 1.0)
        @test (a + b).data == [11.0, 22.0, 33.0]
        @test (2 * a).data == [2.0, 4.0, 6.0]

        mismatched_df = FrequencySeries([1.0, 2.0, 3.0], 0.0, 2.0)
        @test_throws ArgumentError a + mismatched_df

        mismatched_f0 = FrequencySeries([1.0, 2.0, 3.0], 1.0, 1.0)
        @test_throws ArgumentError a + mismatched_f0
    end

    @testset "FrequencySeries show" begin
        fs = FrequencySeries(ComplexF64[1, 2], 0.0, 1.0; name="bar")
        s = sprint(show, fs)
        @test occursin("FrequencySeries", s)
        @test occursin("bar", s)
    end

    @testset "GWOSC loader" begin
        strain_path = find_test_data(GW150914_H1_HDF5)

        if strain_path === nothing
            @test_skip "GW150914 H1 HDF5 file not found (see test/testdata.jl)"
        else
            strain, dq_mask = load_gwosc_strain(strain_path)
            @test strain isa TimeSeries{Float64}
            @test length(strain) == 16_777_216
            @test strain.dt ≈ 0.000244140625
            @test strain.name == "strain"

            @test dq_mask isa TimeSeries
            @test dq_mask.t0 == strain.t0
            @test length(dq_mask) == 4096   # 4096s at 1 Hz

            seg = condition_segment(strain; trigger_time=1126259462.4)
            @test seg isa TimeSeries{Float64}
            @test length(seg) == round(Int, 4.0 / strain.dt)
            @test seg.t0 ≈ 1126259462.4 - 2.0

            # A Tukey (alpha=0.2) window is exactly zero at its first/last sample by
            # construction, regardless of the underlying (unwindowed) strain value there.
            @test seg.data[1] == 0.0
            @test seg.data[end] == 0.0

            @test_throws ArgumentError condition_segment(strain; trigger_time=1.0)  # out of range
        end
    end

    @testset "resample" begin
        @test_throws ArgumentError resample(TimeSeries([1.0, 2.0], 0.0, 1.0), 0.0)
        @test_throws ArgumentError resample(TimeSeries([1.0, 2.0], 0.0, 1.0), -10.0)

        @testset "identity rate is a no-op" begin
            ts = TimeSeries([1.0, 2.0, 3.0, 4.0], 5.0, 0.1; name="x")
            r = resample(ts, sample_rate(ts))
            @test r == ts
            @test r.data !== ts.data   # returns a copy, not an alias
        end

        @testset "downsample-then-upsample round trip" begin
            # A single sinusoid well below both the downsampled and original Nyquist frequencies
            # should round-trip through downsample -> upsample with only small edge-transient
            # error (trimmed away below).
            fs_old = 1000.0
            fs_mid = 100.0
            f_sig = 5.0
            N = 4000
            t = (0:(N - 1)) ./ fs_old
            original = sin.(2pi * f_sig .* t)
            ts = TimeSeries(original, 0.0, 1 / fs_old)

            down = resample(ts, fs_mid)
            @test sample_rate(down) == fs_mid
            @test length(down) == round(Int, N * fs_mid / fs_old)

            back_up = resample(down, fs_old)
            @test sample_rate(back_up) == fs_old
            @test length(back_up) == length(ts)

            # Trim ~10% of samples off each edge to exclude the FIR filter's ramp-up/-down
            # transient (a real, expected artifact of finite-length filtering, not a bug).
            trim = round(Int, 0.1 * N)
            interior = (trim + 1):(N - trim)
            resid = back_up.data[interior] .- original[interior]
            @test maximum(abs.(resid)) < 0.02
        end

        @testset "Nyquist-respecting: above-new-Nyquist content is attenuated, not aliased" begin
            # A pure sinusoid above the new Nyquist frequency (but below the original Nyquist,
            # so it's a legitimate signal at the original rate) should come out heavily
            # attenuated after downsampling -- NOT aliased down into an unattenuated sinusoid at
            # the folded frequency |f - fs_new * round(f / fs_new)|.
            fs_old = 100.0
            fs_new = 20.0            # new Nyquist = 10 Hz
            f_sig = 13.0             # > new Nyquist (10 Hz), < old Nyquist (50 Hz)
            N = 2000
            t = (0:(N - 1)) ./ fs_old
            ts = TimeSeries(sin.(2pi * f_sig .* t), 0.0, 1 / fs_old)

            down = resample(ts, fs_new)
            trim = round(Int, 0.1 * length(down))
            interior_amp = maximum(abs.(down.data[(trim + 1):(end - trim)]))

            # Naive (unfiltered) decimation would alias 13 Hz to |13 - 20| = 7 Hz with amplitude
            # ~1.0 (unattenuated); the anti-aliasing filter should suppress it well below that.
            @test interior_amp < 0.05
        end

        @testset "non-integer-rate ratio (arbitrary resampler path)" begin
            ts = TimeSeries(sin.(2pi * 5.0 .* ((0:999) ./ 1000.0)), 0.0, 1 / 1000.0)
            down = resample(ts, 137.5)
            @test sample_rate(down) == 137.5
            @test down.dt == 1 / 137.5
        end
    end

    @testset "Spectrogram" begin
        fs = 4096.0
        dur = 4.0
        N = round(Int, dur * fs)
        t = (0:(N - 1)) ./ fs
        t0 = 1000.0

        @testset "multi-tone: frequency axis and peak bins" begin
            f1, f2 = 50.0, 200.0
            x = sin.(2π * f1 .* t) .+ sin.(2π * f2 .* t)
            ts = TimeSeries(x, t0, 1 / fs; name="multitone")

            n = 1024
            spec = spectrogram(ts; n=n)
            @test spec isa Spectrogram{Float64}
            @test size(spec, 1) == n ÷ 2 + 1          # one-sided real FFT
            @test spec.f0 == 0.0
            @test spec.df ≈ fs / n

            fvals = collect(frequencies(spec))
            @test fvals[1] == 0.0
            @test fvals[end] ≈ fs / 2

            # Every time column should show its two largest peaks at (approximately) f1 and f2.
            for j in 1:size(spec, 2)
                col = spec[:, j]
                order = sortperm(col.data; rev=true)
                peak_freqs = sort(fvals[order[1:2]])
                @test isapprox(peak_freqs[1], f1; atol=2 * spec.df)
                @test isapprox(peak_freqs[2], f2; atol=2 * spec.df)
            end
        end

        @testset "chirp: peak frequency increases over time" begin
            f_start, f_end = 50.0, 800.0
            k = (f_end - f_start) / dur
            phase = 2π .* (f_start .* t .+ k .* t .^ 2 ./ 2)
            x = sin.(phase)
            ts = TimeSeries(x, t0, 1 / fs; name="chirp")

            n = 512
            spec = spectrogram(ts; n=n)
            fvals = collect(frequencies(spec))

            peak_freq_per_col = [fvals[argmax(spec[:, j].data)] for j in 1:size(spec, 2)]
            # Instantaneous frequency should rise roughly monotonically from near f_start to
            # near f_end. Individual columns can jitter by a bin or two (STFT time/frequency
            # resolution trade-off), so check the trend on quartile-averaged peaks rather than
            # requiring every consecutive pair to be non-decreasing.
            @test peak_freq_per_col[1] < peak_freq_per_col[end]
            @test isapprox(peak_freq_per_col[1], f_start; atol=50.0)
            @test isapprox(peak_freq_per_col[end], f_end; atol=50.0)
            ncols = length(peak_freq_per_col)
            edges = round.(Int, range(0, ncols; length=5))
            quartile_means = [Statistics.mean(peak_freq_per_col[(edges[i] + 1):edges[i + 1]])
                               for i in 1:4]
            @test issorted(quartile_means)
        end

        @testset "GPS time axis" begin
            x = sin.(2π * 100.0 .* t)
            ts = TimeSeries(x, t0, 1 / fs)
            n, noverlap = 1024, 512
            spec = spectrogram(ts; n=n, noverlap=noverlap)

            tvals = collect(times(spec))
            hop = (n - noverlap) / fs
            @test step(times(spec)) ≈ hop
            # DSP.spectrogram's first segment is centered at n/(2*fs) seconds after the start of
            # the input signal.
            @test tvals[1] ≈ t0 + n / (2 * fs)
            @test tvals[end] ≈ t0 + n / (2 * fs) + (size(spec, 2) - 1) * hop
        end

        @testset "indexing" begin
            x = sin.(2π * 100.0 .* t)
            ts = TimeSeries(x, t0, 1 / fs; name="tone")
            spec = spectrogram(ts; n=1024)

            col = spec[:, 1]
            @test col isa FrequencySeries
            @test col.f0 == spec.f0
            @test col.df == spec.df
            @test length(col) == size(spec, 1)

            row = spec[1, :]
            @test row isa TimeSeries
            @test row.t0 == spec.t0
            @test row.dt == spec.dt
            @test length(row) == size(spec, 2)

            @test spec[3, 2] == spec.data[3, 2]
        end

        @testset "show" begin
            x = sin.(2π * 100.0 .* t)
            ts = TimeSeries(x, t0, 1 / fs; name="tone")
            spec = spectrogram(ts; n=1024)
            s = sprint(show, spec)
            @test occursin("Spectrogram", s)
            @test occursin("freq bins", s)
        end
    end
end

include("psd_tests.jl")
include("coherence_tests.jl")
include("filtering_whitening.jl")
include("gps_segments.jl")
include("gwf.jl")
include("test_detector_response.jl")
include("qtransform.jl")
include("plotting_recipes_tests.jl")
include("reference_curves_tests.jl")
include("fetch_tests.jl")
include("posteriors_tests.jl")   # uses fetch_tests.jl's _gwosc_reachable
