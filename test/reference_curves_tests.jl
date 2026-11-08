using Test
using GWData

# ================================================================================================
# Reference ASD/PSD design-sensitivity curves (aLIGO zero-detuning-high-power, AdvVirgo).
#
# The reference values below marked "cross-validated" were obtained by running
# examples/validation/compare_reference_curves.jl against
# lalsimulation.SimNoisePSDaLIGOZeroDetHighPower/lalsimulation.SimNoisePSDAdvVirgo. Both curves
# matched lalsimulation to within ~1e-15 relative difference across 5-2048 Hz (essentially Float64
# machine epsilon, as expected since this is a verbatim port of the same closed-form formulas) --
# not a discrepancy, just floating-point evaluation-order noise.
# ================================================================================================

@testset "ReferenceCurves: aligo_design_psd/asd" begin
    @testset "vector form: basic properties" begin
        f = collect(10.0:10.0:2048.0)
        p = aligo_design_psd(f)
        a = aligo_design_asd(f)

        @test length(p) == length(f)
        @test all(isfinite, p)
        @test all(>(0), p)
        @test a ≈ sqrt.(p)

        # Design curve should have a broad minimum somewhere in the ~100-300 Hz "bucket" band and
        # rise steeply away from it at both low and high frequency (shot noise dominates at high
        # f, suspension/coating thermal + radiation pressure dominate at low f).
        imin = argmin(p)
        fmin = f[imin]
        @test 80.0 <= fmin <= 400.0
        @test p[1] > p[imin]      # 10 Hz is well above the minimum
        @test p[end] > p[imin]    # 2048 Hz is well above the minimum
    end

    @testset "FrequencySeries form" begin
        fs = aligo_design_psd(0.25, 100)
        @test fs isa FrequencySeries
        @test fs.f0 == 0.0
        @test fs.df == 0.25
        @test length(fs) == 100
        @test fs.data[1] == 0.0   # DC bin set to exactly 0.0, not evaluated (see module docstring)
        @test all(>=(0), fs.data)
        @test all(isfinite, fs.data)

        # matches the plain-vector form away from the DC bin
        f = collect(frequencies(fs))
        @test fs.data[2:end] ≈ aligo_design_psd(f)[2:end]

        fs2 = aligo_design_psd(100.0, 0.25, 50)   # explicit f0
        @test fs2.f0 == 100.0
        @test fs2.data ≈ aligo_design_psd(collect(frequencies(fs2)))

        asd_fs = aligo_design_asd(0.25, 100)
        @test asd_fs.data ≈ sqrt.(fs.data)
        @test asd_fs.f0 == fs.f0 && asd_fs.df == fs.df
    end

    @testset "cross-validated spot checks" begin
        # A handful of frequencies from examples/validation's compare_reference_curves.jl run,
        # reproduced here as a standalone regression check (rtol gives generous headroom above the
        # ~1e-15 max relative difference observed against lalsimulation in that run).
        @test aligo_design_psd([10.0])[1] ≈ 9.466734353467914e-45 rtol=1e-10
        @test aligo_design_psd([50.0])[1] ≈ 2.475897902328087e-47 rtol=1e-10
        @test aligo_design_psd([250.0])[1] ≈ 1.3012458026243426e-47 rtol=1e-10   # near the minimum
        @test aligo_design_psd([1000.0])[1] ≈ 2.731538274974741e-47 rtol=1e-10
        @test aligo_design_psd([2000.0])[1] ≈ 7.847384300710755e-47 rtol=1e-10
    end
end

@testset "ReferenceCurves: avirgo_design_psd/asd" begin
    @testset "vector form: basic properties" begin
        f = collect(10.0:10.0:2048.0)
        p = avirgo_design_psd(f)
        a = avirgo_design_asd(f)

        @test length(p) == length(f)
        @test all(isfinite, p)
        @test all(>(0), p)
        @test a ≈ sqrt.(p)

        imin = argmin(p)
        fmin = f[imin]
        @test 100.0 <= fmin <= 500.0
        @test p[1] > p[imin]
        @test p[end] > p[imin]
    end

    @testset "FrequencySeries form" begin
        fs = avirgo_design_psd(0.25, 100)
        @test fs isa FrequencySeries
        @test fs.f0 == 0.0
        @test fs.df == 0.25
        @test fs.data[1] == 0.0
        @test all(>=(0), fs.data)
        @test all(isfinite, fs.data)

        f = collect(frequencies(fs))
        @test fs.data[2:end] ≈ avirgo_design_psd(f)[2:end]

        asd_fs = avirgo_design_asd(0.25, 100)
        @test asd_fs.data ≈ sqrt.(fs.data)
    end

    @testset "cross-validated spot checks" begin
        @test avirgo_design_psd([10.0])[1] ≈ 1.2641275030900394e-42 rtol=1e-10
        @test avirgo_design_psd([50.0])[1] ≈ 4.1489025238697056e-47 rtol=1e-10
        @test avirgo_design_psd([300.0])[1] ≈ 1.1478495865573972e-47 rtol=1e-10   # near the minimum
        @test avirgo_design_psd([1000.0])[1] ≈ 3.495040221230733e-46 rtol=1e-10
        @test avirgo_design_psd([2000.0])[1] ≈ 2.3166528016959212e-45 rtol=1e-10
    end
end
