# Tests for `fetch_open_data` (GWData.GWOSC).
#
# Split into two groups:
#   - Network-independent: the local-caching logic (`_gwosc_download`) and cache-directory
#     resolution, tested against synthetic in-memory/on-disk data only.
#   - Network-dependent: real live calls to https://gwosc.org (no stubbing/mocking), gated behind a
#     short reachability probe so the suite degrades to `@test_skip` (not a hard failure) in a
#     network-isolated environment rather than hanging or erroring.

using Downloads
using GWData.GWOSC: _gwosc_download

function _gwosc_reachable()
    try
        io = IOBuffer()
        Downloads.download("https://gwosc.org/archive/1/2/json/", io; timeout=15)
        return true
    catch
        return false
    end
end

@testset "GWOSC fetch: cache directory resolution" begin
    withenv("GWDATA_CACHE_DIR" => "/tmp/gwd", "XDG_CACHE_HOME" => "/tmp/xdg") do
        @test gwdata_cache_dir() == "/tmp/gwd"
        @test gwdata_cache_dir("strain") == joinpath("/tmp/gwd", "strain")
    end
    withenv("GWDATA_CACHE_DIR" => nothing, "XDG_CACHE_HOME" => "/tmp/xdg") do
        @test gwdata_cache_dir("strain") == joinpath("/tmp/xdg", "gwdata", "strain")
    end
    withenv("GWDATA_CACHE_DIR" => nothing, "XDG_CACHE_HOME" => nothing) do
        @test gwdata_cache_dir() == joinpath(homedir(), ".cache", "gwdata")
    end
end

@testset "GWOSC fetch: local caching (no network)" begin
    @testset "existing cached file is returned without touching the network" begin
        cache_dir = mktempdir()
        fake_url = "https://gwosc.invalid.example/does/not/exist/fake-1234-1.hdf5"
        dest = joinpath(cache_dir, "fake-1234-1.hdf5")
        write(dest, "cached-bytes")

        got = _gwosc_download(fake_url, cache_dir)
        @test got == dest
        @test read(dest, String) == "cached-bytes"   # untouched -- proves no download was attempted
    end

    @testset "force=true re-downloads, and a failed download doesn't corrupt the cache" begin
        cache_dir = mktempdir()
        fake_url = "https://gwosc.invalid.example/does/not/exist/fake-5678-1.hdf5"
        dest = joinpath(cache_dir, "fake-5678-1.hdf5")
        write(dest, "original-bytes")

        @test_throws ErrorException _gwosc_download(fake_url, cache_dir; force=true)
        @test read(dest, String) == "original-bytes"   # untouched on failure
        @test !isfile(dest * ".part")                  # no leftover partial-download file
    end
end

@testset "GWOSC fetch: argument validation (no network)" begin
    @test_throws ArgumentError fetch_open_data("H1", 100.0, 100.0)   # gps_start == gps_end
    @test_throws ArgumentError fetch_open_data("H1", 200.0, 100.0)   # gps_start > gps_end
end

@testset "GWOSC fetch: live GWOSC network calls" begin
    if !_gwosc_reachable()
        @test_skip "GWOSC (https://gwosc.org) unreachable from this environment -- skipping live fetch tests"
    else
        cache_dir = mktempdir()

        @testset "GPS-range fetch reproduces the known-good GW150914 file" begin
            path = fetch_open_data("H1", 1126256640, 1126256640 + 4096; cache_dir=cache_dir)
            @test isfile(path)

            # Strong regression check: byte-identical to a previously obtained copy of the same
            # file, when one is available locally (see test/testdata.jl).
            known_path = find_test_data(GW150914_H1_HDF5)
            if known_path === nothing
                @test_skip "reference GW150914 H1 file not found; skipping byte-identical check"
            else
                @test read(path) == read(known_path)
            end

            # Independently of that, the fetched file must itself parse correctly and
            # have the expected, known shape for this file.
            strain, dq_mask = load_gwosc_strain(path)
            @test length(strain) == 16_777_216
            @test strain.dt ≈ 0.000244140625
            @test all(isfinite, strain.data)
            @test length(dq_mask) == 4096
        end

        @testset "repeat GPS-range fetch is served from the local cache" begin
            path1 = fetch_open_data("H1", 1126256640, 1126256640 + 4096; cache_dir=cache_dir)
            mtime1 = mtime(path1)
            sleep(0.05)
            path2 = fetch_open_data("H1", 1126256640, 1126256640 + 4096; cache_dir=cache_dir)
            @test path2 == path1
            @test mtime(path2) == mtime1   # file wasn't re-downloaded/re-written
        end

        @testset "event-name fetch: GW190412 32s live smoke test" begin
            path = fetch_open_data("GW190412", "H1"; sample_rate=4096, duration=32,
                                    cache_dir=cache_dir)
            @test isfile(path)

            strain, dq_mask = load_gwosc_strain(path)
            @test length(strain) == 32 * 4096
            @test strain.dt ≈ 1 / 4096
            @test all(isfinite, strain.data)
        end

        @testset "GPS range spanning more than one science-run segment errors clearly" begin
            @test_throws ArgumentError fetch_open_data("H1", 1126256640, 1126256640 + 5000;
                                                         cache_dir=cache_dir)
        end

        @testset "unreachable host raises a clear error, not a hang" begin
            @test_throws ErrorException fetch_open_data("H1", 1126256640, 1126256640 + 4096;
                                                          host="https://gwosc.invalid.example",
                                                          cache_dir=cache_dir, force=true)
        end

        @testset "unknown event name raises a clear error" begin
            @test_throws ArgumentError fetch_open_data("GW_not_a_real_event_name_xyz", "H1";
                                                         cache_dir=cache_dir)
        end
    end
end
