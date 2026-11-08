using Test
using GWData
using GWData.Detectors: H1, L1, V1

# ================================================================================================
# Sidereal time, detector geometry, antenna response.
#
# The reference values below marked "cross-validated" were obtained by running
# examples/validation/compare_gmst.jl and examples/validation/compare_antenna_response.jl
# against lal.GreenwichMeanSiderealTime/lal.CachedDetectors/lal.ComputeDetAMResponse. GMST matched
# lal bit-for-bit (0.0 absolute difference) at every GPS time tried, including a leap-second
# boundary and fractional-second inputs. Antenna response matched to within ~1.7e-8 absolute /
# ~1.2e-6 relative -- attributable to LALSuite's own `LALDetector.response` being stored as REAL4
# (single precision) internally (see src/GWData/Detectors.jl's module docstring), not a discrepancy
# in this port.
# ================================================================================================

@testset "SiderealTime: gmst" begin
    # cross-validated (bit-for-bit match against lal.GreenwichMeanSiderealTime).
    @test GWData.gmst(1126259462.4) ≈ 36137.05523464343 rtol=0 atol=0  # GW150914
    @test GWData.gmst(1187008882.4) ≈ 40566.97324959053 rtol=0 atol=0  # GW170817
    @test GWData.gmst(630763213.0) ≈ 4.894961212823059 rtol=0 atol=0   # J2000.0 epoch

    # Leap-second boundary (2015-Jul-01, GPS 1119744016, TAI-UTC 35->36): GPS second 1119744015
    # is UTC 23:59:59, GPS second 1119744016 is the inserted leap second UTC 23:59:60 (still
    # 2015-06-30, per XLALGPSToUTC's `utc->tm_sec += 1` NOT renormalizing the rest of the civil
    # time struct), and GPS second 1119744017 is UTC 2015-07-01 00:00:00. GMST must advance by
    # ~1 sidereal second from 1119744015 to 1119744016 (a normal step), but 1119744016
    # (23:59:60) and 1119744017 (00:00:00 next day) land on the exact same Julian Day under
    # XLALConvertCivilTimeToJD's calendar arithmetic (incrementing the day-of-month field by one
    # vs. adding a full day's worth of seconds are equivalent under that formula) and so must be
    # bit-identical -- confirmed against lal in examples/validation/compare_gmst.jl (0.0
    # absolute difference at all three GPS times), not an artifact of this port.
    g0 = GWData.gmst(1119744015.0)
    g1 = GWData.gmst(1119744016.0)
    g2 = GWData.gmst(1119744017.0)
    sidereal_step = g0 - GWData.gmst(1119744014.0)  # a normal, non-boundary 1-second GMST step
    @test g1 - g0 ≈ sidereal_step atol=1e-9   # 23:59:59 -> 23:59:60: an ordinary step
    @test g2 == g1                             # 23:59:60 == 00:00:00 (next day), bit-for-bit

    # GMST increases monotonically with GPS time over a full day (no reduction mod 2*pi, matching
    # LALSuite's own convention).
    @test GWData.gmst(1126259462.4 + 86400.0) > GWData.gmst(1126259462.4)

    # Sidereal day (~23h56m4.0905s) is slightly shorter than a solar day: GMST should advance by
    # slightly MORE than 2*pi over exactly one solar (86400s) day.
    delta = GWData.gmst(1126259462.4 + 86400.0) - GWData.gmst(1126259462.4)
    @test delta > 2 * pi
    @test delta ≈ 2 * pi * (86400.0 / 86164.0905) rtol=1e-3
end

@testset "SiderealTime: leap second table sanity" begin
    @test_throws ArgumentError GWData.SiderealTime._leap_seconds_tai_utc(-100000000)  # before 1980
    @test GWData.SiderealTime._leap_seconds_tai_utc(1119744015) == 35   # just before the 2015 leap second
    @test GWData.SiderealTime._leap_seconds_tai_utc(1119744016) == 36   # at/after it
    @test GWData.SiderealTime._leap_seconds_tai_utc(1440000000) == 37   # holds the last table value
end

@testset "Detectors: geometry sanity" begin
    for det in (H1, L1, V1)
        loc = geometry(det).location
        r_earth = sqrt(sum(x -> x^2, loc))
        # Earth's mean radius is ~6.371e6 m; all real detector sites are within a few hundred
        # meters of sea level, so the vertex location vector's norm must be close to that.
        @test 6.3e6 < r_earth < 6.4e6

        # response = 0.5*(X ox X - Y ox Y) for unit vectors X, Y => exactly traceless
        # (up to floating-point roundoff in the constants' sum-of-squares).
        resp = geometry(det).response
        tr = resp[1][1] + resp[2][2] + resp[3][3]
        @test isapprox(tr, 0.0; atol=1e-10)

        # response must be symmetric (it's an outer-product combination by construction).
        for i in 1:3, j in 1:3
            @test resp[i][j] ≈ resp[j][i]
        end
    end

    # Different detectors must have different locations/orientations.
    @test geometry(H1).location != geometry(L1).location
    @test geometry(H1).response != geometry(L1).response
end

# --- Antenna-response sanity checks: source directly overhead vs. along an arm (a natural
# --- "horizon-like" in-plane direction for this geometry), derived from first principles rather
# --- than an external reference, and exact for any detector whose two arms are exactly
# --- perpendicular unit vectors (true for H1/L1/V1 to within their tiny real-world arm-alignment
# --- imperfections).
#
# The wave's propagation direction, in XLALComputeDetAMResponse's own (gha, dec) parametrization,
# is Z = (-cos(gha)*cos(dec), sin(gha)*cos(dec), -sin(dec)) (see DetResponse.c's
# XLALComputeDetAMResponseParts/ExtraModes -- XLALComputeDetAMResponse itself doesn't return Z,
# but uses the same X/Y wave-frame construction). Inverting this gives the (ra, dec) that makes
# the source arrive from a chosen Earth-fixed unit direction `v` at a given `gps`.
function _radec_for_propagation_direction(v, gps)
    g = GWData.gmst(gps)
    dec = -asin(clamp(v[3], -1.0, 1.0))
    gha = atan(v[2], -v[1])
    return g - gha, dec
end

_cross(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
_normalize(a) = (n = sqrt(sum(x -> x^2, a)); a ./ n)

# H1's arm-direction unit vectors (LALDetectors.h LAL_LHO_4K_ARM_{X,Y}_DIRECTION_*, same
# constants src/GWData/Detectors.jl builds `H1`'s response tensor from).
const _H1_ARM_X = (-0.22389266154, 0.79983062746, 0.55690487831)
const _H1_ARM_Y = (-0.91397818574, 0.02609403989, -0.40492342125)

@testset "Detectors: isbits geometry and generic element types" begin
    g = geometry(H1)
    @test isbitstype(typeof(g))
    @test g isa DetectorGeometry{Float64}
    @test geometry(g) === g
    @test H1.prefix == "H1" && L1.prefix == "L1" && V1.prefix == "V1"
    @test sprint(show, H1) == "Detector(\"H1\")"

    ra, dec, psi, gps = 1.0, 0.5, 0.3, 1126259462.4
    gm = gmst(gps)
    # The Detector, its geometry and the gps/gmst entry points all agree exactly.
    @test antenna_response_gmst(H1, ra, dec, psi, gm) === antenna_response_gmst(g, ra, dec, psi, gm)
    @test compute_antenna_response(H1, ra, dec, psi, gps) === antenna_response_gmst(g, ra, dec, psi, gm)
    @test time_delay_from_geocenter(L1, ra, dec, gps) === time_delay_from_geocenter_gmst(geometry(L1), ra, dec, gm)

    # Float32 geometry + Float32 inputs stay Float32 (no silent promotion on GPU paths). `gmst` is
    # not reduced mod 2pi (~3.6e4 rad here), so it must be reduced before rounding to Float32.
    g32 = DetectorGeometry{Float32}(H1)
    @test isbitstype(typeof(g32))
    gm32 = Float32(mod2pi(gm))
    fp32, fc32 = antenna_response_gmst(g32, Float32(ra), Float32(dec), Float32(psi), gm32)
    dt32 = time_delay_from_geocenter_gmst(g32, Float32(ra), Float32(dec), gm32)
    @test fp32 isa Float32 && fc32 isa Float32 && dt32 isa Float32
    fp, fc = antenna_response_gmst(g, ra, dec, psi, gm)
    @test isapprox(fp32, fp; atol=1e-5) && isapprox(fc32, fc; atol=1e-5)
    @test isapprox(dt32, time_delay_from_geocenter_gmst(g, ra, dec, gm); rtol=1e-5)
end

@testset "compute_antenna_response: overhead vs. along-an-arm sanity" begin
    gps = 1126259462.4
    up = _normalize(_cross(_H1_ARM_X, _H1_ARM_Y))

    # Directly overhead (propagation along the arm-plane normal): for a differential-arm
    # detector with (near-)perpendicular unit arm vectors, F+^2 + Fx^2 = 1 EXACTLY here,
    # independent of psi -- maximal, polarization-insensitive sensitivity at normal incidence.
    ra_up, dec_up = _radec_for_propagation_direction(up, gps)
    for psi in (0.0, 0.3, 1.0, pi / 4)
        fp, fc = compute_antenna_response(H1, ra_up, dec_up, psi, gps)
        @test fp^2 + fc^2 ≈ 1.0 atol=1e-8
    end

    # Nadir (straight through the Earth, -up): same magnitude by symmetry of the response tensor.
    ra_dn, dec_dn = _radec_for_propagation_direction(.-(up), gps)
    fp, fc = compute_antenna_response(H1, ra_dn, dec_dn, 0.0, gps)
    @test fp^2 + fc^2 ≈ 1.0 atol=1e-8

    # Along either arm direction (propagation direction lying in the plane of the arms -- the
    # "horizon" of the detector's own local frame, as opposed to overhead): a strictly and
    # substantially smaller, but still exact, response of F+^2 + Fx^2 = 1/4, independent of psi.
    ra_x, dec_x = _radec_for_propagation_direction(_normalize(_H1_ARM_X), gps)
    for psi in (0.0, 0.3, 1.0, pi / 4)
        fp, fc = compute_antenna_response(H1, ra_x, dec_x, psi, gps)
        @test fp^2 + fc^2 ≈ 0.25 atol=1e-8
    end

    ra_y, dec_y = _radec_for_propagation_direction(_normalize(_H1_ARM_Y), gps)
    fp, fc = compute_antenna_response(H1, ra_y, dec_y, 0.0, gps)
    @test fp^2 + fc^2 ≈ 0.25 atol=1e-8

    # Overhead response is strictly, substantially larger than the in-plane ("horizon") response.
    fp_up, fc_up = compute_antenna_response(H1, ra_up, dec_up, 0.0, gps)
    fp_x, fc_x = compute_antenna_response(H1, ra_x, dec_x, 0.0, gps)
    @test fp_up^2 + fc_up^2 > fp_x^2 + fc_x^2
end

@testset "compute_antenna_response: cross-validated spot checks" begin
    # A handful of (detector, ra, dec, psi, gps) points from examples/validation's
    # compare_antenna_response.jl run, reproduced here as a standalone regression check (rtol
    # gives headroom above the ~1.2e-6 max relative difference observed against lal in that run).
    @test all(isapprox.(
        compute_antenna_response(H1, 1.0, 0.5, 0.3, 1126259462.4),
        (-0.40139482390491776, 0.7361070064503537); rtol=1e-4))
    @test all(isapprox.(
        compute_antenna_response(L1, 3.0, -1.0, 1.2, 1187008882.4),
        (0.5934610478282492, 0.15788814143897917); rtol=1e-4))
    @test all(isapprox.(
        compute_antenna_response(V1, 5.9, 1.4, 2.9, 1440000000.0),
        (-0.2153138779273707, 0.6571851719114025); rtol=1e-4))
end

# ================================================================================================
# time_delay_from_geocenter, ported from LALSuite's XLALTimeDelayFromEarthCenter (via its
# XLALArrivalTimeDiff helper, specialized to detector2 = earth_center = (0,0,0)) --
# lal/lib/tools/TimeDelay.c, fetched 2026-09-01 from
# https://lscsoft.docs.ligo.org/lalsuite/dev/lal/_time_delay_8c_source.html (the canonical
# git.ligo.org raw-source URL 404s; this versioned doxygen mirror serves the identical `dev`
# (HEAD-of-`master`) source). Sign convention (matches LALSuite exactly): positive when the
# wavefront arrives at the detector AFTER the geocenter.
#
# cross-validated (examples/validation/compare_time_delay.jl) against
# lal.TimeDelayFromEarthCenter across H1/L1/V1 x 8 (ra, dec, gps) points (generic spread + near-
# pole + leap-second-boundary + J2000.0 edge cases, same point convention as the antenna-response
# comparison above): **bit-for-bit exact, 0.0 absolute/relative difference at every one of the 24
# points** -- even tighter than compute_antenna_response's ~1e-6 (which is limited by LALSuite's
# own REAL4-precision response tensor; XLALTimeDelayFromEarthCenter has no such precision loss).
# ================================================================================================

@testset "time_delay_from_geocenter: cross-validated spot checks" begin
    @test time_delay_from_geocenter(H1, 0.0, 0.0, 1126259462.4) ≈ -0.013676276399281918 rtol=0 atol=0
    @test time_delay_from_geocenter(H1, 1.0, 0.5, 1126259462.4) ≈ -0.017787561403509365 rtol=0 atol=0
    @test time_delay_from_geocenter(H1, 3.0, -1.0, 1187008882.4) ≈ 0.018516313630107306 rtol=0 atol=0
    @test time_delay_from_geocenter(H1, 0.0, pi / 2 - 0.001, 1126259462.4) ≈ -0.015358785211000089 rtol=0 atol=0
    @test time_delay_from_geocenter(H1, 2.5, -(pi / 2 - 0.001), 1126259462.4) ≈ 0.015352867647758311 rtol=0 atol=0
    @test time_delay_from_geocenter(H1, 4.2, 0.05, 1119744016.0) ≈ -0.002980577300737952 rtol=0 atol=0
    @test time_delay_from_geocenter(H1, 5.9, 1.4, 1440000000.0) ≈ -0.014080172308287965 rtol=0 atol=0
    @test time_delay_from_geocenter(H1, 0.7, -0.3, 630763213.0) ≈ 0.011748038931490389 rtol=0 atol=0

    @test time_delay_from_geocenter(L1, 0.0, 0.0, 1126259462.4) ≈ -0.011791897822641196 rtol=0 atol=0
    @test time_delay_from_geocenter(L1, 1.0, 0.5, 1126259462.4) ≈ -0.02111577008450708 rtol=0 atol=0
    @test time_delay_from_geocenter(L1, 3.0, -1.0, 1187008882.4) ≈ 0.011831565958731374 rtol=0 atol=0
    @test time_delay_from_geocenter(L1, 4.2, 0.05, 1119744016.0) ≈ -0.011638529866987834 rtol=0 atol=0

    @test time_delay_from_geocenter(V1, 0.0, 0.0, 1126259462.4) ≈ 0.013522698419534863 rtol=0 atol=0
    @test time_delay_from_geocenter(V1, 1.0, 0.5, 1126259462.4) ≈ -0.006067968101849518 rtol=0 atol=0
    @test time_delay_from_geocenter(V1, 5.9, 1.4, 1440000000.0) ≈ -0.013267666101930383 rtol=0 atol=0
    @test time_delay_from_geocenter(V1, 0.7, -0.3, 630763213.0) ≈ 0.009147822106336849 rtol=0 atol=0
end

@testset "time_delay_from_geocenter: structural sanity" begin
    const_c = 299792458.0
    gps = 1126259462.4

    # Direct port of the (ra, dec) <-> source-direction-unit-vector relationship inverting
    # XLALArrivalTimeDiff's own ehat_src = (cos(dec)*cos(gha), -cos(dec)*sin(gha), sin(dec)),
    # gha = gmst - ra -- independent little helper (not reusing this function's own internals),
    # so this checks the port's *geometry*, not just its self-consistency.
    function radec_for_source_direction(u, t)
        g = GWData.gmst(t)
        dec = asin(clamp(u[3], -1.0, 1.0))
        gha = atan(-u[2], u[1])
        return g - gha, dec
    end

    for det in (H1, L1, V1)
        loc = geometry(det).location
        r = sqrt(sum(x -> x^2, loc))
        up = loc ./ r   # unit vector from geocenter toward this detector

        # Source directly overhead the detector (ehat_src == up, i.e. detector lies exactly
        # between the geocenter and the source): the wave reaches the detector first, so the
        # delay is exactly -r/c -- the most negative value possible for this detector, since
        # |dot(ehat_src, location)| <= r for any unit ehat_src (Cauchy-Schwarz).
        ra_up, dec_up = radec_for_source_direction(up, gps)
        @test time_delay_from_geocenter(det, ra_up, dec_up, gps) ≈ -r / const_c atol=1e-9

        # Source at nadir (ehat_src == -up, geocenter lies between source and detector): the wave
        # reaches the geocenter first, then the detector -- delay is exactly +r/c, the most
        # positive value possible.
        ra_dn, dec_dn = radec_for_source_direction(.-up, gps)
        @test time_delay_from_geocenter(det, ra_dn, dec_dn, gps) ≈ r / const_c atol=1e-9

        # |delay| is bounded by r/c (light-travel time across the full geocenter-to-detector
        # distance) for any sky position -- Cauchy-Schwarz on the dot product again, checked here
        # across a spread of arbitrary (ra, dec) rather than just the two extremal directions above.
        for (ra, dec) in ((0.3, 0.7), (2.1, -0.4), (5.5, 0.0), (1.0, 1.4))
            @test abs(time_delay_from_geocenter(det, ra, dec, gps)) <= r / const_c + 1e-12
        end

        # Antipodal sky positions give exactly opposite-sign, equal-magnitude delays (ehat_src
        # flips sign, and the dot product is linear in ehat_src).
        ra1, dec1 = 1.3, 0.6
        ra2, dec2 = ra1 + pi, -dec1   # antipodal direction on the sky
        d1 = time_delay_from_geocenter(det, ra1, dec1, gps)
        d2 = time_delay_from_geocenter(det, ra2, dec2, gps)
        @test d1 ≈ -d2 atol=1e-12
    end

    # Different detectors give different delays for the same generic sky position/time (distinct
    # locations, no accidental degeneracy at this point).
    @test time_delay_from_geocenter(H1, 1.0, 0.5, gps) != time_delay_from_geocenter(L1, 1.0, 0.5, gps)
    @test time_delay_from_geocenter(H1, 1.0, 0.5, gps) != time_delay_from_geocenter(V1, 1.0, 0.5, gps)
end
