# Tests for published posterior-sample download and reading (GravitationalWaveData.Posteriors).

using HDF5
using Statistics: median
using GravitationalWaveData.GWOSC: _url_filename
using GravitationalWaveData.Posteriors: _is_hdf5_url

const GW150914_PE_RELEASE = "IGWN-GWTC2p1-v2-GW150914_095045_PEDataRelease_mixed_cosmo.h5"

@testset "Posteriors: URL handling (no network)" begin
    @test _url_filename("https://gwosc.org/archive/data/O1/1126170624/H-H1_LOSC_4_V1-1126256640-4096.hdf5") ==
          "H-H1_LOSC_4_V1-1126256640-4096.hdf5"
    @test _url_filename("https://zenodo.org/api/records/6513631/files/$(GW150914_PE_RELEASE)/content") ==
          GW150914_PE_RELEASE
    @test _url_filename("https://example.org/a/b.h5?download=1") == "b.h5"

    @test _is_hdf5_url("https://zenodo.org/api/records/1/files/x.h5/content")
    @test _is_hdf5_url("https://dcc.ligo.org/public/0157/P1800370/005/GW170817_GWTC-1.hdf5")
    @test !_is_hdf5_url("https://dcc.ligo.org/public/0169/P2000223/007/GW190412.tar")
    @test !_is_hdf5_url("https://git.ligo.org/some/repo")
end

@testset "Posteriors: PESummary reader (synthetic file)" begin
    path = joinpath(mktempdir(), "synthetic_pesummary.h5")
    rows = [(chirp_mass=30.0 + i, mass_ratio=0.5 + i / 100) for i in 1:5]
    h5open(path, "w") do file
        write(file, "C01:Mixed/posterior_samples", rows)
        write(file, "C01:IMRPhenomXPHM/posterior_samples", rows[1:2])
        write(file, "version/pesummary", "test")   # non-analysis group, must be ignored
    end

    @test list_posterior_analyses(path) == ["C01:IMRPhenomXPHM", "C01:Mixed"]

    samples = load_posterior_samples(path)   # defaults to C01:Mixed
    @test sort(collect(keys(samples))) == ["chirp_mass", "mass_ratio"]
    @test samples["chirp_mass"] == [31.0, 32.0, 33.0, 34.0, 35.0]

    @test length(load_posterior_samples(path, "C01:IMRPhenomXPHM")["mass_ratio"]) == 2
    @test_throws ArgumentError load_posterior_samples(path, "C01:NotThere")
end

@testset "Posteriors: real GW150914 GWTC-2.1 release" begin
    path = joinpath(gravitationalwavedata_cache_dir("posteriors"), GW150914_PE_RELEASE)
    if !isfile(path)
        @test_skip "GW150914 PE release not cached (run fetch_posterior_samples(\"GW150914\"))"
    else
        analyses = list_posterior_analyses(path)
        @test analyses == ["C01:IMRPhenomXPHM", "C01:Mixed", "C01:SEOBNRv4PHM"]
        samples = load_posterior_samples(path, "C01:IMRPhenomXPHM")
        @test length(samples["chirp_mass"]) == 147_634
        # GWTC-2.1 GW150914 source-frame chirp mass is ~28.6 Msun.
        @test 27.0 < median(samples["chirp_mass_source"]) < 30.5
    end
end

@testset "Posteriors: live GWOSC lookup" begin
    if !_gwosc_reachable()
        @test_skip "GWOSC (https://gwosc.org) unreachable -- skipping live posterior lookup"
    else
        # Resolves to the GWTC-2.1 release; served from the cache if already downloaded.
        path = fetch_posterior_samples("GW150914")
        @test basename(path) == GW150914_PE_RELEASE
        @test "C01:Mixed" in list_posterior_analyses(path)
        @test_throws ArgumentError fetch_posterior_samples("GW000000")
    end
end
