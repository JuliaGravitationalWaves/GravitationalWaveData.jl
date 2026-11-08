using Test
using GWData
using RecipesBase

# ================================================================================================
# Tests for the `GWDataRecipesBaseExt` package extension (ext/GWDataRecipesBaseExt.jl).
#
# No plotting backend (Plots.jl/Makie/etc.) is needed to test a RecipesBase recipe: calling
# `RecipesBase.apply_recipe(Dict{Symbol,Any}(), obj)` directly runs the `@recipe` block and
# returns a `Vector{RecipesBase.RecipeData}`, whose `.plotattributes` (the attribute defaults the
# recipe set with `-->`/`:=`) and `.args` (the actual x/y/[z] data returned) can be asserted on
# directly -- this is the standard way to unit-test a recipe without a rendering backend as a test
# dependency. Loading `RecipesBase` here (a `[extras]`/test-only dependency, alongside the
# `[weakdeps]` entry in Project.toml) is what triggers Julia's package-extension mechanism to
# compile and load `GWDataRecipesBaseExt` in this test session.
# ================================================================================================

@testset "Plotting recipes (GWDataRecipesBaseExt)" begin
    @testset "TimeSeries recipe" begin
        ts = TimeSeries([1.0, 2.0, 3.0, 4.0], 100.0, 0.5; name="H1:STRAIN")
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), ts)
        @test length(rd) == 1
        attrs = rd[1].plotattributes
        @test attrs[:seriestype] == :path
        @test attrs[:xlabel] == "GPS Time (s)"
        @test attrs[:ylabel] == "Amplitude"
        @test attrs[:label] == "H1:STRAIN"

        x, y = rd[1].args
        @test collect(x) == collect(times(ts))
        @test y == ts.data
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
        # Qualified as `GWData.spectrogram` (not the bare, exported name): by the time this file
        # runs, `test/psd_tests.jl` has already done a plain `using DSP` earlier in the same flat
        # `include()` chain into `Main`, and `DSP` also exports a `spectrogram` -- the bare name
        # is ambiguous at that point regardless of include order within this file itself.
        spec = GWData.spectrogram(ts; n=64, name="tone spectrogram")

        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), spec)
        @test length(rd) == 1
        attrs = rd[1].plotattributes
        @test attrs[:seriestype] == :heatmap
        @test attrs[:xlabel] == "GPS Time (s)"
        @test attrs[:ylabel] == "Frequency (Hz)"
        @test attrs[:colorbar_scale] == :log10
        @test attrs[:colorbar_title] == "tone spectrogram"

        xarg, yarg, zarg = rd[1].args
        @test collect(xarg) == collect(times(spec))
        @test collect(yarg) == collect(frequencies(spec))
        @test zarg == spec.data
        @test size(zarg) == (size(spec, 1), size(spec, 2))   # (nfreq, ntime): z[i,j] at (x[j], y[i])
    end

    @testset "Spectrogram recipe: unnamed spectrogram gets generic colorbar title" begin
        ts = TimeSeries(sin.(2π * 10.0 .* ((0:511) ./ 100.0)), 0.0, 1 / 100.0)
        spec = GWData.spectrogram(ts; n=64)
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), spec)
        @test rd[1].plotattributes[:colorbar_title] == "Power"
    end

    @testset "Spectrogram recipe: user override wins over recipe default" begin
        ts = TimeSeries(sin.(2π * 10.0 .* ((0:511) ./ 100.0)), 0.0, 1 / 100.0)
        spec = GWData.spectrogram(ts; n=64)
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(:colorbar_scale => :identity), spec)
        @test rd[1].plotattributes[:colorbar_scale] == :identity
    end
end
