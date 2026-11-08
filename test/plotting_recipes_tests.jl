using Test
using GravitationalWaveData
using RecipesBase

# ================================================================================================
# Tests for the `GravitationalWaveDataRecipesBaseExt` package extension (ext/GravitationalWaveDataRecipesBaseExt.jl).
#
# No plotting backend (Plots.jl/Makie/etc.) is needed to test a RecipesBase recipe: calling
# `RecipesBase.apply_recipe(Dict{Symbol,Any}(), obj)` directly runs the `@recipe` block and
# returns a `Vector{RecipesBase.RecipeData}`, whose `.plotattributes` (the attribute defaults the
# recipe set with `-->`/`:=`) and `.args` (the actual x/y/[z] data returned) can be asserted on
# directly -- this is the standard way to unit-test a recipe without a rendering backend as a test
# dependency. Loading `RecipesBase` here (a `[extras]`/test-only dependency, alongside the
# `[weakdeps]` entry in Project.toml) is what triggers Julia's package-extension mechanism to
# compile and load `GravitationalWaveDataRecipesBaseExt` in this test session.
# ================================================================================================

@testset "Plotting recipes (GravitationalWaveDataRecipesBaseExt)" begin
    @testset "TimeSeries recipe" begin
        ts = TimeSeries([1.0, 2.0, 3.0, 4.0], 100.0, 0.5; name="H1:STRAIN")
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), ts)
        @test length(rd) == 1
        attrs = rd[1].plotattributes
        @test attrs[:seriestype] == :path
        @test attrs[:xlabel] == "Time (s) from GPS 100"
        @test attrs[:xformatter](101.5) == "1.5"
        @test attrs[:xformatter](102.0) == "2"
        @test attrs[:ylabel] == "Amplitude"
        @test attrs[:label] == "H1:STRAIN"

        x, y = rd[1].args
        @test collect(x) == collect(times(ts))
        @test y == ts.data
    end

    @testset "Time axis: epoch is round at the scale of the data" begin
        epoch_label(t0, dt, n) = RecipesBase.apply_recipe(Dict{Symbol,Any}(),
                                                           TimeSeries(zeros(n), t0, dt))[1].plotattributes[:xlabel]
        @test epoch_label(1126259454.4, 1.0, 17) == "Time (s) from GPS 1126259450"
        @test epoch_label(1126259462.25, 0.01, 21) == "Time (s) from GPS 1126259462.2"
        @test epoch_label(1126256640.0, 1.0, 4097) == "Time (s) from GPS 1126256000"
        @test epoch_label(1126259462.0, 1.0, 1) == "Time (s) from GPS 1126259462"

        # The label is never in scientific notation, and the x data stays absolute GPS time.
        ts = TimeSeries(zeros(21), 1126259462.25, 0.01)
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), ts)
        @test rd[1].args[1][1] == 1126259462.25
        @test rd[1].plotattributes[:xformatter](1126259462.3) == "0.1"
    end

    @testset "epoch keyword shifts x data" begin
        ts = TimeSeries(zeros(4), 1126259462.0, 0.5)
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(:epoch => 1126259462.4), ts)
        attrs = rd[1].plotattributes
        @test attrs[:xlabel] == "Time (s) from GPS 1126259462.4"
        @test !haskey(attrs, :xformatter)
        @test !haskey(attrs, :epoch)     # consumed here, so Plots.jl never sees an unknown attribute
        @test rd[1].args[1] ≈ [-0.4, 0.1, 0.6, 1.1] atol=1e-6   # Float64 spacing at GPS scale

        spec = Spectrogram(ones(2, 4) .+ [1 2 3 4; 5 6 7 8], 1126259462.0, 0.5, 0.0, 1.0)
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(:epoch => 1126259462), spec)
        @test rd[1].plotattributes[:xlabel] == "Time (s) from GPS 1126259462"
        @test !haskey(rd[1].plotattributes, :epoch)
        @test rd[1].args[1] ≈ [0.0, 0.5, 1.0, 1.5]
    end

    @testset "TimeSeries recipe: unnamed series gets empty label" begin
        ts = TimeSeries([1.0, 2.0], 0.0, 1.0)
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), ts)
        @test rd[1].plotattributes[:label] == ""
    end

    @testset "TimeSeries recipe: complex data plotted as magnitude" begin
        ts = TimeSeries(ComplexF64[3 + 4im, 0 + 5im], 0.0, 1.0)
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), ts)
        _, y = rd[1].args
        @test y == [5.0, 5.0]
        @test eltype(y) <: Real
    end

    @testset "TimeSeries recipe: user override wins over recipe default" begin
        ts = TimeSeries([1.0, 2.0], 0.0, 1.0; name="x")
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(:label => "custom", :ylabel => "Strain"), ts)
        @test rd[1].plotattributes[:label] == "custom"
        @test rd[1].plotattributes[:ylabel] == "Strain"
    end

    @testset "FrequencySeries recipe" begin
        fs = FrequencySeries([1.0, 2.0, 3.0, 4.0], 0.0, 0.25; name="H1 PSD")
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), fs)
        @test length(rd) == 1
        attrs = rd[1].plotattributes
        @test attrs[:seriestype] == :path
        @test attrs[:xlabel] == "Frequency (Hz)"
        @test attrs[:ylabel] == "Amplitude"
        @test attrs[:xscale] == :log10
        @test attrs[:yscale] == :log10
        @test attrs[:label] == "H1 PSD"

        x, y = rd[1].args
        @test collect(x) == collect(frequencies(fs))
        @test y == fs.data
    end

    @testset "FrequencySeries recipe: complex data plotted as magnitude" begin
        fs = FrequencySeries(ComplexF64[3 + 4im, 0 - 5im], 0.0, 1.0)
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), fs)
        _, y = rd[1].args
        @test y == [5.0, 5.0]
        @test eltype(y) <: Real
    end

    @testset "FrequencySeries recipe: user override wins over recipe default" begin
        fs = FrequencySeries([1.0, 2.0], 0.0, 1.0)
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(:xscale => :identity), fs)
        @test rd[1].plotattributes[:xscale] == :identity
        @test rd[1].plotattributes[:yscale] == :log10   # untouched default still applies
    end

    @testset "Spectrogram recipe" begin
        fs_rate = 100.0
        N = 800
        t = (0:(N - 1)) ./ fs_rate
        x = sin.(2π * 10.0 .* t)
        ts = TimeSeries(x, 1000.0, 1 / fs_rate; name="tone")
        # Qualified as `GravitationalWaveData.spectrogram` (not the bare, exported name): by the time this file
        # runs, `test/psd_tests.jl` has already done a plain `using DSP` earlier in the same flat
        # `include()` chain into `Main`, and `DSP` also exports a `spectrogram` -- the bare name
        # is ambiguous at that point regardless of include order within this file itself.
        spec = GravitationalWaveData.spectrogram(ts; n=64, name="tone spectrogram")

        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), spec)
        @test length(rd) == 1
        attrs = rd[1].plotattributes
        @test attrs[:seriestype] == :heatmap
        @test attrs[:xlabel] == "Time (s) from GPS 1000"
        @test attrs[:ylabel] == "Frequency (Hz)"
        @test !haskey(attrs, :colorbar_scale)
        @test attrs[:colorbar_title] == "tone spectrogram (dB)"

        xarg, yarg, zarg = rd[1].args
        @test collect(xarg) == collect(times(spec))
        @test collect(yarg) == collect(frequencies(spec))
        positive = spec.data .> 0
        @test zarg[positive] ≈ 10 .* log10.(spec.data[positive])
        @test all(isnan, zarg[.!positive])
        @test size(zarg) == (size(spec, 1), size(spec, 2))   # (nfreq, ntime): z[i,j] at (x[j], y[i])
    end

    @testset "Spectrogram recipe: color limits ignore outliers and zero power" begin
        data = [1e-40 * (1 + (i + j) / 200) for i in 1:100, j in 1:100]   # about -400 to -397 dB
        data[1, 1] = 1e-20    # lone bright pixel
        data[2, 1] = 1e-80    # lone dark pixel
        data[3, 1] = 0.0      # no power: log10 would be -Inf
        spec = Spectrogram(data, 0.0, 1.0, 1.0, 1.0)
        attrs = RecipesBase.apply_recipe(Dict{Symbol,Any}(), spec)[1].plotattributes
        low, high = attrs[:clims]
        @test -400.0 <= low < high <= -397.0

        # A flat image has no percentile range, so limits are left to the backend, and a caller's
        # limits win.
        flat = Spectrogram(ones(4, 4), 0.0, 1.0, 1.0, 1.0)
        @test !haskey(RecipesBase.apply_recipe(Dict{Symbol,Any}(), flat)[1].plotattributes, :clims)
        chosen = RecipesBase.apply_recipe(Dict{Symbol,Any}(:clims => (-5, 5)), spec)[1].plotattributes
        @test chosen[:clims] == (-5, 5)
    end

    @testset "Spectrogram recipe: unnamed spectrogram gets generic colorbar title" begin
        ts = TimeSeries(sin.(2π * 10.0 .* ((0:511) ./ 100.0)), 0.0, 1 / 100.0)
        spec = GravitationalWaveData.spectrogram(ts; n=64)
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), spec)
        @test rd[1].plotattributes[:colorbar_title] == "Power (dB)"
    end

    @testset "Spectrogram recipe: user override wins over recipe default" begin
        ts = TimeSeries(sin.(2π * 10.0 .* ((0:511) ./ 100.0)), 0.0, 1 / 100.0)
        spec = GravitationalWaveData.spectrogram(ts; n=64)
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(:colorbar_scale => :identity), spec)
        attrs = rd[1].plotattributes
        @test attrs[:colorbar_scale] == :identity
        # Asking for a colorbar scale means raw power, on that scale.
        @test rd[1].args[3] == spec.data
        @test attrs[:colorbar_title] == "Power"
        @test !haskey(attrs, :clims)
    end
end
