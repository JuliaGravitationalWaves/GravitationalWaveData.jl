module Detectors

using ..SiderealTime: gmst

export AbstractDetector, Detector, DetectorGeometry, geometry, H1, L1, V1
export compute_antenna_response, antenna_response_gmst, time_delay_from_geocenter, time_delay_from_geocenter_gmst

# ================================================================================================
# Detector geometry (location vector + response tensor), ported VERBATIM from LALSuite's
# hardcoded detector table, and antenna response (F+, Fx), ported from XLALComputeDetAMResponse.
#
# Source (fetched 2026-08-28 from the canonical LALSuite repository, commit at HEAD of the
# `master` branch at https://git.ligo.org/lscsoft/lalsuite):
#   - lal/lib/tools/LALDetectors.h -- lalCachedDetectors' underlying #define constants:
#       LAL_{LHO_4K,LLO_4K,VIRGO}_VERTEX_LOCATION_{X,Y,Z}_SI  (Earth-centered Cartesian vertex
#         position, meters)
#       LAL_{LHO_4K,LLO_4K,VIRGO}_ARM_{X,Y}_DIRECTION_{X,Y,Z} (Earth-centered Cartesian unit
#         vectors along each arm)
#     These are LALSuite's own hardcoded constants (not re-derived from lat/lon/azimuth here) --
#     copied faithfully rather than recomputed, so this module can't drift from LALSuite's numbers
#     even if e.g. WGS84 ellipsoid parameters were ported slightly differently elsewhere.
#   - lal/lib/tools/CreateDetector.c -- confirms the response-tensor convention used to build
#     `lalCachedDetectors[...].response` for an `LALDETECTORTYPE_IFODIFF` detector (H1/L1/V1 are
#     all this type): d^ab = (1/2)(x^a x^b - y^a y^b) where x, y are the unit arm-direction
#     vectors above. LALSuite stores this tensor as REAL4 (single precision) internally; this
#     port keeps it as Float64 since nothing here needs bit-identical float32 rounding, and the
#     cross-validation tolerance (see examples/validation/) comfortably absorbs the
#     resulting ~1e-7 relative difference.
#   - lal/lib/tools/DetResponse.c -- XLALComputeDetAMResponse (the antenna-pattern formula
#     itself, Anderson et al. PRD 63 042003 (2001) eqs. B4/B5/B7).
#
# K1 (KAGRA) is not included: LALDetectors.h only gives KAGRA's geodetic
# (lat/lon/elevation/arm-azimuth/arm-altitude) parameters, not the hardcoded Cartesian vertex/arm
# constants H1/L1/V1 get -- reproducing K1 bit-for-bit would additionally require porting LALSuite's
# WGS84-ellipsoid geodetic-to-Cartesian conversion (`getCartesianComponents` in CreateDetector.c)
# and cross-validating that conversion itself, which has not been done.
#
# Two detector types share one interface (`AbstractDetector`, unwrapped by `geometry`):
#   - `DetectorGeometry{T}` holds only the numbers (location, response tensor). It is isbits, so it
#     can be built inside a GPU kernel and passed through Enzyme-differentiated code.
#   - `Detector` adds the host-side `prefix` ("H1") used for file names, config keys and
#     calibration priors. Its `String` field makes it non-isbits: GPU/Enzyme code takes
#     `geometry(det)` instead.
# Every function below is written once against `DetectorGeometry` and type-generic in its scalar
# arguments, with all tuple indices literal (runtime tuple indexing does not compile on GPU).
# ================================================================================================

"""
    AbstractDetector

Supertype of [`Detector`](@ref) (named, host-side) and [`DetectorGeometry`](@ref) (isbits, usable
inside GPU kernels and Enzyme-differentiated code). Every detector function accepts either;
[`geometry`](@ref) returns the underlying `DetectorGeometry`.
"""
abstract type AbstractDetector end

"""
    DetectorGeometry{T}(location, response)

A gravitational-wave interferometer's geometry: Earth-centered Cartesian `location` vector (meters)
and `response` tensor (dimensionless, the `d^{ab}` of Anderson et al. PRD 63 042003 (2001), stored
as `response[i][j]`). Plain bits, so it can be passed to or built inside GPU kernels.

`DetectorGeometry{Float32}(H1)` converts a detector's geometry to another element type.
"""
struct DetectorGeometry{T<:Real} <: AbstractDetector
    location::NTuple{3,T}
    response::NTuple{3,NTuple{3,T}}
end

DetectorGeometry{T}(g::DetectorGeometry) where {T<:Real} =
    DetectorGeometry{T}(map(T, g.location), map(r -> map(T, r), g.response))
DetectorGeometry{T}(d::AbstractDetector) where {T<:Real} = DetectorGeometry{T}(geometry(d))

"""
    Detector(prefix, geometry::DetectorGeometry{Float64})

A named interferometer: `prefix` is the two-character GWOSC/channel-name prefix (e.g. `"H1"`),
`geometry` its [`DetectorGeometry`](@ref). See [`H1`](@ref)/[`L1`](@ref)/[`V1`](@ref) for the three
ready-made instances this package provides, ported verbatim from LALSuite's `lalCachedDetectors`
table (see this file's module docstring for the exact source).

Not isbits (the `prefix` is a `String`): pass `geometry(detector)` into GPU kernels or
Enzyme-differentiated code.
"""
struct Detector <: AbstractDetector
    prefix::String
    geometry::DetectorGeometry{Float64}
end

"""
    geometry(detector::AbstractDetector) -> DetectorGeometry

The isbits [`DetectorGeometry`](@ref) of `detector` (itself, if it already is one).
"""
@inline geometry(g::DetectorGeometry) = g
@inline geometry(d::Detector) = d.geometry

Base.show(io::IO, d::Detector) = print(io, "Detector(\"", d.prefix, "\")")

# --- Response tensor construction: d^{ab} = (1/2)(x^a x^b - y^a y^b), CreateDetector.c's
# --- LALDETECTORTYPE_IFODIFF case, applied here to LALSuite's own hardcoded arm-direction
# --- constants (rather than re-deriving the arm vectors from lat/lon/azimuth/altitude).
function _response_tensor(x::NTuple{3,Float64}, y::NTuple{3,Float64})
    d = ntuple(i -> ntuple(j -> 0.5 * (x[i] * x[j] - y[i] * y[j]), 3), 3)
    return d
end

# LHO_4k (H1). LALDetectors.h constants:
#   LAL_LHO_4K_VERTEX_LOCATION_{X,Y,Z}_SI, LAL_LHO_4K_ARM_{X,Y}_DIRECTION_{X,Y,Z}
const _H1_LOCATION = (-2.16141492636e+06, -3.83469517889e+06, 4.60035022664e+06)
const _H1_ARM_X = (-0.22389266154, 0.79983062746, 0.55690487831)
const _H1_ARM_Y = (-0.91397818574, 0.02609403989, -0.40492342125)

# LLO_4k (L1). LALDetectors.h constants:
#   LAL_LLO_4K_VERTEX_LOCATION_{X,Y,Z}_SI, LAL_LLO_4K_ARM_{X,Y}_DIRECTION_{X,Y,Z}
const _L1_LOCATION = (-7.42760447238e+04, -5.49628371971e+06, 3.22425701744e+06)
const _L1_ARM_X = (-0.95457412153, -0.14158077340, -0.26218911324)
const _L1_ARM_Y = (0.29774156894, -0.48791033647, -0.82054461286)

# VIRGO (V1). LALDetectors.h constants:
#   LAL_VIRGO_VERTEX_LOCATION_{X,Y,Z}_SI, LAL_VIRGO_ARM_{X,Y}_DIRECTION_{X,Y,Z}
const _V1_LOCATION = (4.54637409900e+06, 8.42989697626e+05, 4.37857696241e+06)
const _V1_ARM_X = (-0.70045821479, 0.20848948619, 0.68256166277)
const _V1_ARM_Y = (-0.05379255368, -0.96908180549, 0.24080451708)

"""
    H1 :: Detector

LIGO Hanford Observatory's 4km interferometer (LALSuite's "LHO_4k"). Geometry (vertex location,
arm-direction unit vectors used to build `response`) ported verbatim from LALSuite's
`lalCachedDetectors` table (`LAL_LHO_4K_*` constants in `lal/lib/tools/LALDetectors.h`) -- see this
file's module docstring for the exact source and provenance.
"""
const H1 = Detector("H1", DetectorGeometry(_H1_LOCATION, _response_tensor(_H1_ARM_X, _H1_ARM_Y)))

"""
    L1 :: Detector

LIGO Livingston Observatory's 4km interferometer (LALSuite's "LLO_4k"). Geometry ported verbatim
from LALSuite's `lalCachedDetectors` table (`LAL_LLO_4K_*` constants in
`lal/lib/tools/LALDetectors.h`) -- see this file's module docstring for the exact source and
provenance. See [`H1`](@ref) for the general construction.
"""
const L1 = Detector("L1", DetectorGeometry(_L1_LOCATION, _response_tensor(_L1_ARM_X, _L1_ARM_Y)))

"""
    V1 :: Detector

Virgo Observatory's 3km interferometer (LALSuite's "VIRGO"). Geometry ported verbatim from
LALSuite's `lalCachedDetectors` table (`LAL_VIRGO_*` constants in `lal/lib/tools/LALDetectors.h`)
-- see this file's module docstring for the exact source and provenance. See [`H1`](@ref) for the
general construction.
"""
const V1 = Detector("V1", DetectorGeometry(_V1_LOCATION, _response_tensor(_V1_ARM_X, _V1_ARM_Y)))

"""
    compute_antenna_response(detector::AbstractDetector, ra, dec, psi, gps_time) -> (Fplus, Fcross)

Antenna pattern functions `F+`/`Fx` for a gravitational wave from sky position `(ra, dec)`
(equatorial right ascension/declination, radians), polarization angle `psi` (radians), arriving
at `gps_time` (GPS seconds). Direct port of `XLALComputeDetAMResponse`
(`lal/lib/tools/DetResponse.c`), implementing Eqs. (B4), (B5), (B7) of Anderson et al., PRD 63,
042003 (2001).

Computes [`gmst`](@ref GravitationalWaveData.SiderealTime.gmst)`(gps_time)` (relatively expensive: leap-second
table scan + Julian-day conversion) and calls [`antenna_response_gmst`](@ref), which callers that
already hold the GMST -- and GPU kernels -- should use directly.
"""
function compute_antenna_response(detector::AbstractDetector, ra::Real, dec::Real, psi::Real, gps_time::Real)
    return antenna_response_gmst(detector, ra, dec, psi, gmst(gps_time))
end

"""
    antenna_response_gmst(detector::AbstractDetector, ra, dec, psi, gmst) -> (Fplus, Fcross)

[`compute_antenna_response`](@ref) with the Greenwich Mean Sidereal Time (radians) given instead of
the GPS time. Allocation-free and type-generic (the result's type follows the promoted argument
and geometry types), so it is safe inside GPU kernels and Enzyme-differentiated code when called on
a [`DetectorGeometry`](@ref).

[`gmst`](@ref GravitationalWaveData.SiderealTime.gmst) is not reduced modulo 2pi (it is ~4e4 rad for current GPS
times): pass `Float32(mod2pi(gmst(gps_time)))` for single-precision evaluation, or the hour angle
loses ~1e-3 rad.
"""
@inline antenna_response_gmst(detector::Detector, ra, dec, psi, g) =
    antenna_response_gmst(geometry(detector), ra, dec, psi, g)

@inline function antenna_response_gmst(detector::DetectorGeometry, ra, dec, psi, g)
    gha = g - ra  # Greenwich hour angle of the source, radians

    cosgha = cos(gha)
    singha = sin(gha)
    cosdec = cos(dec)
    sindec = sin(dec)
    cospsi = cos(psi)
    sinpsi = sin(psi)

    # Eq. (B4) of [ABCF2001].
    X = (
        -cospsi * singha - sinpsi * cosgha * sindec,
        -cospsi * cosgha + sinpsi * singha * sindec,
        sinpsi * cosdec,
    )
    # Eq. (B5) of [ABCF2001].
    Y = (
        sinpsi * singha - cospsi * cosgha * sindec,
        sinpsi * cosgha + cospsi * singha * sindec,
        cospsi * cosdec,
    )

    # Eq. (B7) of [ABCF2001], contracted with the response tensor. Manually unrolled: a runtime
    # loop index into a tuple lowers to a Julia runtime call GPUCompiler cannot compile.
    r1, r2, r3 = detector.response
    DX1 = r1[1] * X[1] + r1[2] * X[2] + r1[3] * X[3]
    DY1 = r1[1] * Y[1] + r1[2] * Y[2] + r1[3] * Y[3]
    DX2 = r2[1] * X[1] + r2[2] * X[2] + r2[3] * X[3]
    DY2 = r2[1] * Y[1] + r2[2] * Y[2] + r2[3] * Y[3]
    DX3 = r3[1] * X[1] + r3[2] * X[2] + r3[3] * X[3]
    DY3 = r3[1] * Y[1] + r3[2] * Y[2] + r3[3] * Y[3]
    fplus = X[1] * DX1 - Y[1] * DY1 + X[2] * DX2 - Y[2] * DY2 + X[3] * DX3 - Y[3] * DY3
    fcross = X[1] * DY1 + Y[1] * DX1 + X[2] * DY2 + Y[2] * DX2 + X[3] * DY3 + Y[3] * DX3

    return (fplus, fcross)
end

# --- Speed of light, m/s (exact SI value). Ported from LALSuite's `LAL_C_SI`
# --- (`lal/lib/std/LALConstants.h`). Also used by `ReferenceCurves.jl`.
const _C_SI = 299792458.0

"""
    time_delay_from_geocenter(detector::AbstractDetector, ra, dec, gps_time) -> Float64

Light-travel-time offset (seconds) between the geocenter and `detector`, for a gravitational wave
arriving from sky position `(ra, dec)` (equatorial right ascension/declination, radians) at
`gps_time` (GPS seconds). Direct port of `XLALTimeDelayFromEarthCenter`
(`lal/lib/tools/TimeDelay.c`), via its `XLALArrivalTimeDiff` helper specialized to
`detector2 = earth_center = (0,0,0)`.

Sign convention (matches LALSuite exactly): **positive when the wavefront arrives at the detector
after arriving at the geocenter** -- i.e. `t_detector = t_geocenter + time_delay_from_geocenter(...)`.
Computes [`gmst`](@ref GravitationalWaveData.SiderealTime.gmst)`(gps_time)` and calls
[`time_delay_from_geocenter_gmst`](@ref).
"""
function time_delay_from_geocenter(detector::AbstractDetector, ra::Real, dec::Real, gps_time::Real)
    return time_delay_from_geocenter_gmst(detector, ra, dec, gmst(gps_time))
end

"""
    time_delay_from_geocenter_gmst(detector::AbstractDetector, ra, dec, gmst) -> delay

[`time_delay_from_geocenter`](@ref) with the Greenwich Mean Sidereal Time (radians) given instead of
the GPS time. Allocation-free and type-generic, like [`antenna_response_gmst`](@ref).
"""
@inline time_delay_from_geocenter_gmst(detector::Detector, ra, dec, g) =
    time_delay_from_geocenter_gmst(geometry(detector), ra, dec, g)

@inline function time_delay_from_geocenter_gmst(detector::DetectorGeometry{T}, ra, dec, g) where {T}
    gha = g - ra  # Greenwich hour angle of the source, radians

    cosdec = cos(dec)
    # Unit vector from the geocenter to the source, in the Earth-fixed frame (XLALArrivalTimeDiff's
    # `ehat_src`).
    ehat_src = (cosdec * cos(gha), -cosdec * sin(gha), sin(dec))

    loc = detector.location
    # XLALArrivalTimeDiff(detector, earth_center=(0,0,0), ra, dec, gps) computes
    # dot(ehat_src, earth_center - detector) / c == -dot(ehat_src, detector.location) / c.
    dp = ehat_src[1] * loc[1] + ehat_src[2] * loc[2] + ehat_src[3] * loc[3]
    return -dp / T(_C_SI)
end

end # module Detectors
