# Tests for the GWF (IGWD frame) reader in `src/GWData/IO/GWF.jl`.
#
# Primary validation: a real GWOSC frame file, `H-H1_LOSC_4_V1-1126256640-4096.gwf` (see
# `data/README.md` for how to obtain it), covering the same GPS segment (1126256640, 4096s, H1,
# 4 kHz) as the HDF5 file `load_gwosc_strain` reads. This lets the GWF-parsed strain be checked
# directly against the HDF5 loader on identical data, not just internally.

@testset "GWF frame-file reader" begin
    gwf_path = something(find_test_data(GW150914_H1_GWF), "")
    hdf5_path = something(find_test_data(GW150914_H1_HDF5), "")

    @testset "internal helpers" begin
        # _read_string on a small hand-built buffer: STRING = INT_2U length (incl. trailing '\0')
        # + that many bytes of ASCII ending in '\0'. Independent of any real/synthetic frame file.
        buf = UInt8[0x04, 0x00, UInt8('a'), UInt8('b'), UInt8('c'), 0x00, 0xFF]
        s, p = GWData.GWF._read_string(buf, 1)
        @test s == "abc"
        @test p == 7   # consumed 2 (length prefix) + 4 (3 chars + null) = 6 bytes, next pos = 7
        @test buf[p] == 0xFF   # byte after the string untouched

        # FrVect type-code -> Julia type table (LIGO-T970130 Appendix C), spot-checked against
        # the codes actually seen in the real file below (type=2 for STRAIN, type=4 for the masks).
        @test GWData.GWF._FRVECT_ELTYPE[0x0002] == Float64
        @test GWData.GWF._FRVECT_ELTYPE[0x0004] == Int32

        # Raw (uncompressed) round trip through _decompress_frvect.
        raw = reinterpret(UInt8, Float64[1.0, 2.0, 3.0]) |> collect
        @test GWData.GWF._decompress_frvect(raw, 256, Float64, 3) == [1.0, 2.0, 3.0]

        # An unimplemented compression scheme errors clearly rather than misreading data.
        @test_throws ArgumentError GWData.GWF._decompress_frvect(UInt8[1, 2, 3, 4], 261, Int16, 2)
    end

    if !isfile(gwf_path)
        @test_skip "GW150914 H1 .gwf file not found (see data/README.md and test/testdata.jl)"
    else
        @testset "list_gwf_channels" begin
            chans = list_gwf_channels(gwf_path)
            @test "H1:LOSC-STRAIN" in chans
            @test "H1:LOSC-DQMASK" in chans
            @test "H1:LOSC-INJMASK" in chans
        end

        @testset "load_gwf_strain: structure" begin
            ts = load_gwf_strain(gwf_path)
            @test ts isa TimeSeries{Float64}
            @test ts.name == "H1:LOSC-STRAIN"
            @test length(ts) == 16_777_216
            @test ts.dt ≈ 0.000244140625
            @test ts.t0 == 1126256640.0

            # Explicit channel= argument selects the same data.
            ts2 = load_gwf_strain(gwf_path; channel="H1:LOSC-STRAIN")
            @test ts2 == ts

            @test_throws ArgumentError load_gwf_strain(gwf_path; channel="NOT-A-REAL-CHANNEL")
        end

        if !isfile(hdf5_path)
            @test_skip "GW150914 H1 HDF5 file not found; skipping GWF-vs-HDF5 cross-check"
        else
            @testset "GWF strain matches HDF5 strain bit-for-bit" begin
                gwf_ts = load_gwf_strain(gwf_path)
                hdf5_ts, _dq = load_gwosc_strain(hdf5_path)

                @test length(gwf_ts) == length(hdf5_ts)
                @test gwf_ts.t0 == hdf5_ts.t0
                @test gwf_ts.dt == hdf5_ts.dt
                # gzip round-trips REAL_8 samples exactly (lossless compression of exact bit
                # patterns) -- this is a genuine bit-for-bit check, not a tolerance-based one.
                @test gwf_ts.data == hdf5_ts.data
            end
        end
    end
end
