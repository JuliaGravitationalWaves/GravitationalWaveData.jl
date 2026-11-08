using Test
using Dates
using GWData

@testset "GPSTime" begin
    @testset "epoch anchor" begin
        @test utc_to_gps(GPS_EPOCH) == 0.0
        @test gps_to_utc(0.0) == GPS_EPOCH
        @test leap_seconds(GPS_EPOCH) == 0
        @test leap_seconds(0.0) == 0
    end

    @testset "round trip across leap-second boundaries" begin
        # A handful of GPS times spanning several historical leap-second insertions. Deliberately
        # avoid values that land exactly on a leap-second insertion instant itself (e.g. a GPS
        # time that is a whole multiple of 86400 seconds past the epoch can coincide with a
        # 00:00:00 UTC leap-second effective date) -- that exact 1-second window genuinely has no
        # unique `DateTime` representation (it *is* the inserted leap second, `23:59:60`, which
        # `Dates.DateTime` cannot represent), so round-tripping through it is not expected to be
        # exact; see the caveat in `gps_to_utc`'s docstring. The fractional-second offsets below
        # keep every value away from any leap-second boundary (always exactly on a whole second).
        for gps in (0.0, 100_000.3, 46_830_000.3, 500_000_000.7, 1_000_000_000.2,
                    1_126_259_462.4, 1_700_000_000.6)
            utc = gps_to_utc(gps)
            @test utc_to_gps(utc) ≈ gps atol=1e-3   # millisecond DateTime resolution
        end
    end

    @testset "leap_seconds monotonically non-decreasing, matches known offsets" begin
        # GPS-UTC was 0 at the 1980 epoch and is 18s today (no leap second since 2016-12-31,
        # effective 2017-01-01) -- see GPSTime.jl's LEAP_SECOND_TABLE header for the source.
        @test leap_seconds(GPS_EPOCH) == 0
        @test leap_seconds(DateTime(2016, 12, 31)) == 17
        @test leap_seconds(DateTime(2017, 1, 1)) == 18
        @test leap_seconds(DateTime(2026, 8, 28)) == 18   # "today", per the task's stated date

        gps_offsets = [leap_seconds(gps) for gps in 0.0:1e8:2e9]
        @test issorted(gps_offsets)
    end

    @testset "GW150914 trigger time spot check" begin
        # GW150914's GPS trigger time is publicly documented (GWOSC event portal,
        # https://www.gw-openscience.org/eventapi/html/GWTC-1-confident/GW150914/, and LIGO's own
        # discovery-paper press materials) as GPS 1126259462.4 == 2015-09-14 09:50:45.4 UTC.
        # Independently cross-checked directly against `gwpy.time.tconvert(1126259462.4)`
        # (the sibling LALSuite.jl repo's `examples/psd/.CondaPkg` env, read-only), which returns
        # exactly `2015-09-14 09:50:45.400000+00:00` -- matches bit-for-bit.
        # At that date (after the 2015-07-01 leap second, before the next one on 2017-01-01),
        # GPS-UTC = 17s, so this is also a check on that specific historical offset.
        gw150914_gps = 1126259462.4
        gw150914_utc = DateTime(2015, 9, 14, 9, 50, 45, 400)

        @test leap_seconds(gw150914_gps) == 17
        @test gps_to_utc(gw150914_gps) == gw150914_utc
        @test utc_to_gps(gw150914_utc) ≈ gw150914_gps atol=1e-3
    end

    @testset "GW170817 trigger time spot check" begin
        # GW170817's GPS trigger time, also publicly documented on the GWOSC event portal
        # (https://www.gw-openscience.org/eventapi/html/GWTC-1-confident/GW170817/): GPS
        # 1187008882.4 == 2017-08-17 12:41:04.4 UTC. Independently cross-checked against
        # `gwpy.time.tconvert(1187008882.4)`, which returns exactly
        # `2017-08-17 12:41:04.400000+00:00` -- matches bit-for-bit. This postdates the
        # 2017-01-01 leap second (GPS-UTC=18s), exercising a different table row than the
        # GW150914 check above.
        gw170817_gps = 1187008882.4
        gw170817_utc = DateTime(2017, 8, 17, 12, 41, 4, 400)

        @test leap_seconds(gw170817_gps) == 18
        @test gps_to_utc(gw170817_gps) == gw170817_utc
        @test utc_to_gps(gw170817_utc) ≈ gw170817_gps atol=1e-3
    end
end

@testset "Segment" begin
    @testset "construction" begin
        s = Segment(0.0, 10.0)
        @test s.start == 0.0
        @test s.stop == 10.0
        @test duration(s) == 10.0
        @test_throws ArgumentError Segment(10.0, 10.0)
        @test_throws ArgumentError Segment(10.0, 5.0)
    end

    @testset "in / overlaps / equality" begin
        s = Segment(10.0, 20.0)
        @test 10.0 in s
        @test 15.0 in s
        @test !(20.0 in s)   # half-open: stop excluded
        @test !(9.9 in s)

        @test overlaps(Segment(0.0, 10.0), Segment(5.0, 15.0))
        @test !overlaps(Segment(0.0, 10.0), Segment(10.0, 20.0))   # touching, not overlapping
        @test !overlaps(Segment(0.0, 10.0), Segment(20.0, 30.0))

        @test Segment(1.0, 2.0) == Segment(1.0, 2.0)
        @test Segment(1.0, 2.0) != Segment(1.0, 3.0)
    end
end

@testset "SegmentList" begin
    @testset "construction, sorting, coalescing" begin
        empty_sl = SegmentList()
        @test isempty(empty_sl)
        @test length(empty_sl) == 0
        @test duration(empty_sl) == 0.0

        # Out-of-order, overlapping, and touching segments should all sort + merge.
        sl = SegmentList([Segment(10.0, 20.0), Segment(0.0, 5.0), Segment(5.0, 8.0),
                           Segment(18.0, 25.0)])
        @test length(sl) == 2
        @test sl[1] == Segment(0.0, 8.0)     # [0,5) and [5,8) touch -> merged
        @test sl[2] == Segment(10.0, 25.0)   # [10,20) and [18,25) overlap -> merged
        @test duration(sl) == 8.0 + 15.0

        # coalesce_segments is an alias for the SegmentList constructor.
        @test coalesce_segments([Segment(0.0, 5.0), Segment(5.0, 8.0)]) ==
              SegmentList([Segment(0.0, 5.0), Segment(5.0, 8.0)])
        @test coalesce_segments([Segment(0.0, 5.0), Segment(5.0, 8.0)]) ==
              SegmentList([Segment(0.0, 8.0)])
    end

    @testset "in / iteration" begin
        sl = SegmentList([Segment(0.0, 5.0), Segment(10.0, 15.0)])
        @test 2.0 in sl
        @test !(7.0 in sl)
        @test 12.0 in sl
        @test collect(sl) == [Segment(0.0, 5.0), Segment(10.0, 15.0)]
    end

    @testset "union" begin
        a = SegmentList([Segment(0.0, 5.0), Segment(10.0, 15.0)])
        b = SegmentList([Segment(4.0, 12.0)])
        u = union(a, b)
        @test u == SegmentList([Segment(0.0, 15.0)])   # bridges the gap between a's two segments

        c = SegmentList([Segment(20.0, 25.0)])
        @test union(a, c) == SegmentList([Segment(0.0, 5.0), Segment(10.0, 15.0),
                                           Segment(20.0, 25.0)])

        # union with an empty list is a no-op (up to coalescing of the non-empty side).
        @test union(a, SegmentList()) == a
    end

    @testset "intersect" begin
        a = SegmentList([Segment(0.0, 10.0), Segment(20.0, 30.0)])
        b = SegmentList([Segment(5.0, 25.0)])
        @test intersect(a, b) == SegmentList([Segment(5.0, 10.0), Segment(20.0, 25.0)])

        # disjoint lists intersect to empty
        @test isempty(intersect(SegmentList([Segment(0.0, 5.0)]),
                                 SegmentList([Segment(5.0, 10.0)])))   # touching, not overlapping

        # a list intersected with itself is itself
        @test intersect(a, a) == a

        # intersect with empty is empty
        @test isempty(intersect(a, SegmentList()))
    end

    @testset "gaps" begin
        sl = SegmentList([Segment(0.0, 5.0), Segment(10.0, 15.0), Segment(20.0, 22.0)])

        # no bound: only internal gaps
        @test gaps(sl) == SegmentList([Segment(5.0, 10.0), Segment(15.0, 20.0)])

        # a single segment has no internal gaps
        @test isempty(gaps(SegmentList([Segment(0.0, 5.0)])))

        # with an explicit bound: leading/trailing gaps included too
        g = gaps(sl, Segment(-5.0, 30.0))
        @test g == SegmentList([Segment(-5.0, 0.0), Segment(5.0, 10.0), Segment(15.0, 20.0),
                                 Segment(22.0, 30.0)])

        # bound narrower than the list's extent: segments get clipped, not just filtered
        g2 = gaps(sl, Segment(3.0, 12.0))
        @test g2 == SegmentList([Segment(5.0, 10.0)])

        # bound with no coverage at all -> the whole bound is one big gap
        @test gaps(SegmentList(), Segment(0.0, 10.0)) == SegmentList([Segment(0.0, 10.0)])

        # bound fully covered -> no gaps
        @test isempty(gaps(SegmentList([Segment(0.0, 10.0)]), Segment(2.0, 8.0)))
    end
end
