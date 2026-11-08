module SiderealTime

using Dates

export gmst, gps_to_utc_seconds

# ================================================================================================
# Greenwich Mean Sidereal Time, ported VERBATIM from LALSuite's C implementation.
#
# Why not a generic astronomy package: modern sidereal-time formulae (IAU-2006/ERFA-style
# precession-nutation models) are more astronomically precise than the ~1982-vintage formula
# LALSuite uses (itself borrowed from the NOVAS-C library), but they do NOT agree with LALSuite
# bit-for-bit -- disagreements land at a level (sub-arcsecond to arcsecond-ish in Earth rotation
# angle) that would silently corrupt sky-location-dependent antenna-response terms downstream.
# Every existing GW pipeline (LALSuite, and everything built on it, including gwpy/PyCBC/Bilby)
# uses LALSuite's exact GMST convention, so THIS is the "correct" one for interoperability
# purposes, not a more modern formula. Hence: port the LALSuite source directly, do not use
# e.g. AstroTime.jl/ERFA.jl's sidereal-time routines.
#
# Source (fetched 2026-08-28 from the canonical LALSuite repository, commit at HEAD of the
# `master` branch at https://git.ligo.org/lscsoft/lalsuite):
#   - lal/lib/date/XLALSiderealTime.c   -- XLALGreenwichSiderealTime/XLALGreenwichMeanSiderealTime
#   - lal/lib/date/XLALCivilTime.c      -- XLALGPSToUTC, XLALConvertCivilTimeToJD, leap-second
#                                           bookkeeping (XLALLeapSeconds/delta_tai_utc)
#   - lal/lib/date/XLALLeapSeconds.h    -- the static `leaps[]` table itself
#   - lal/lib/date/Date.h               -- XLAL_EPOCH_J2000_0_JD / XLAL_EPOCH_UNIX_GPS /
#                                           XLAL_EPOCH_GPS_TAI_UTC constants
# ================================================================================================

const XLAL_EPOCH_J2000_0_JD = 2451545.0
const XLAL_EPOCH_UNIX_GPS = 315964800
const XLAL_EPOCH_GPS_TAI_UTC = 19

# Leap-second table, ported verbatim from lal/lib/date/XLALLeapSeconds.h. Columns: Julian Day
# (UTC) at which the leap second takes effect, GPS second at which it takes effect, and the
# resulting TAI-UTC offset in seconds. LALSuite's own comment: "the table below must be updated
# whenever a leap second is added" -- as of the fetch date above, LALSuite's own upstream table
# stops at 2017-Jan-01 (no leap second has been added since, in reality, as of this writing), so
# reproducing this exact table (not extending it) is what LALSuite-consistency requires.
struct LeapEntry
    jd::Float64
    gpssec::Int
    taiutc::Int
end

const LEAPS = LeapEntry[
    LeapEntry(2444239.5,    -43200, 19),  # 1980-Jan-01
    LeapEntry(2444786.5,  46828800, 20),  # 1981-Jul-01
    LeapEntry(2445151.5,  78364801, 21),  # 1982-Jul-01
    LeapEntry(2445516.5, 109900802, 22),  # 1983-Jul-01
    LeapEntry(2446247.5, 173059203, 23),  # 1985-Jul-01
    LeapEntry(2447161.5, 252028804, 24),  # 1988-Jan-01
    LeapEntry(2447892.5, 315187205, 25),  # 1990-Jan-01
    LeapEntry(2448257.5, 346723206, 26),  # 1991-Jan-01
    LeapEntry(2448804.5, 393984007, 27),  # 1992-Jul-01
    LeapEntry(2449169.5, 425520008, 28),  # 1993-Jul-01
    LeapEntry(2449534.5, 457056009, 29),  # 1994-Jul-01
    LeapEntry(2450083.5, 504489610, 30),  # 1996-Jan-01
    LeapEntry(2450630.5, 551750411, 31),  # 1997-Jul-01
    LeapEntry(2451179.5, 599184012, 32),  # 1999-Jan-01
    LeapEntry(2453736.5, 820108813, 33),  # 2006-Jan-01
    LeapEntry(2454832.5, 914803214, 34),  # 2009-Jan-01
    LeapEntry(2456109.5, 1025136015, 35), # 2012-Jul-01
    LeapEntry(2457204.5, 1119744016, 36), # 2015-Jul-01
    LeapEntry(2457754.5, 1167264017, 37), # 2017-Jan-01
]

"""
    _leap_seconds_tai_utc(gpssec::Integer) -> Int

Port of `XLALLeapSeconds` (`lal/lib/date/XLALCivilTime.c`): returns TAI-UTC (in seconds) in
effect at the given integer GPS second, by scanning the module's `LEAPS` table.
"""
function _leap_seconds_tai_utc(gpssec::Integer)
    gpssec < LEAPS[1].gpssec &&
        throw(ArgumentError("don't know leap seconds before GPS time $(LEAPS[1].gpssec)"))
    leap = 1
    for i in 2:length(LEAPS)
        if gpssec < LEAPS[i].gpssec
            break
        end
        leap = i
    end
    return LEAPS[leap].taiutc
end

"""
    _delta_tai_utc(gpssec::Integer) -> Int

Port of `delta_tai_utc` (`lal/lib/date/XLALCivilTime.c`): the change in TAI-UTC exactly AT
`gpssec` relative to `gpssec - 1` (nonzero only in the single GPS second a leap second is
inserted/removed; 0 otherwise). Used to detect whether UTC needs a 60th ("leap") second.
"""
function _delta_tai_utc(gpssec::Integer)
    for leap in 2:length(LEAPS)
        if gpssec == LEAPS[leap].gpssec
            return LEAPS[leap].taiutc - LEAPS[leap-1].taiutc
        end
    end
    return 0
end

"""
    UTCTime

Broken-down UTC civil time, mirroring the fields of `XLALGPSToUTC`'s C `struct tm` output that
[`gmst`](@ref)'s Julian-day conversion actually needs (year, month, day, hour, minute, second).
`sec` can be 60 for a UTC leap second, matching LALSuite's convention.
"""
struct UTCTime
    year::Int   # full year, e.g. 2015 (not "years since 1900" as in C's struct tm)
    month::Int  # 1-12
    day::Int    # 1-31
    hour::Int
    min::Int
    sec::Int    # 0-60 (60 only during an inserted UTC leap second)
end

"""
    gps_to_utc_seconds(gpssec::Integer) -> UTCTime

Port of `XLALGPSToUTC` (`lal/lib/date/XLALCivilTime.c`): converts an integer GPS second count to
a broken-down UTC civil time, accounting for leap seconds via the module's `LEAPS` table.

LALSuite computes this via `unixsec = gpssec - leapsec + XLAL_EPOCH_GPS_TAI_UTC +
XLAL_EPOCH_UNIX_GPS` fed to `gmtime_r`; Julia's `Dates` provides the equivalent Unix-epoch ->
civil-calendar conversion (proleptic Gregorian, matching glibc's `gmtime_r` in the range GW data
covers) via `Dates.unix2datetime`.
"""
function gps_to_utc_seconds(gpssec::Integer)
    leapsec = _leap_seconds_tai_utc(gpssec)
    unixsec = gpssec - leapsec + XLAL_EPOCH_GPS_TAI_UTC + XLAL_EPOCH_UNIX_GPS
    dt = Dates.unix2datetime(unixsec)
    sec = Dates.second(dt)
    if _delta_tai_utc(gpssec) > 0
        sec += 1  # inserted UTC leap second: 23:59:59 -> 23:59:60
    end
    return UTCTime(Dates.year(dt), Dates.month(dt), Dates.day(dt), Dates.hour(dt), Dates.minute(dt), sec)
end

"""
    _convert_civil_time_to_jd(utc::UTCTime) -> Float64

Port of `XLALConvertCivilTimeToJD` (`lal/lib/date/XLALCivilTime.c`) -- the "short" Julian Day
formula valid for dates since 1900 March (LALSuite's own comment: not expected to be used before
that). Julian days start at noon, hence the `- 0.5`.
"""
function _convert_civil_time_to_jd(utc::UTCTime)
    year = utc.year
    month = utc.month
    day = utc.day
    sec_of_day = utc.sec + 60 * (utc.min + 60 * utc.hour)
    jd = 367 * year - div(7 * (year + div(month + 9, 12)), 4) + div(275 * month, 9) + day + 1721014
    jd += sec_of_day / 86400.0 - 0.5
    return Float64(jd)
end

"""
    gmst(gps_seconds::Real) -> Float64

Greenwich Mean Sidereal Time, in radians, at the given GPS time. Direct port of
`XLALGreenwichMeanSiderealTime`/`XLALGreenwichSiderealTime` (with `equation_of_equinoxes = 0`)
from `lal/lib/date/XLALSiderealTime.c` -- see this file's module docstring for why a generic
astronomy package is deliberately NOT used here.

The result is NOT reduced modulo 2*pi (matching LALSuite's own convention -- callers that need a
value in `[0, 2*pi)` should reduce it themselves, e.g. via `mod2pi`).

`gps_seconds` may carry a fractional part (sub-second precision is handled the same way
LALSuite's `LIGOTimeGPS` integer-seconds/nanoseconds split is: the integer part feeds the
UTC/leap-second/Julian-day machinery which only has whole-second resolution, and the fractional
part is added back in as `t_lo`).
"""
function gmst(gps_seconds::Real)
    gps_int = floor(Int, gps_seconds)
    gps_frac = gps_seconds - gps_int  # fractional second, in [0, 1)

    utc = gps_to_utc_seconds(gps_int)
    julian_day = _convert_civil_time_to_jd(utc)

    t_hi = (julian_day - XLAL_EPOCH_J2000_0_JD) / 36525.0
    t_lo = gps_frac / (36525.0 * 86400.0)
    t = t_hi + t_lo

    sidereal_time = (-6.2e-6 * t + 0.093104) * t * t + 67310.54841
    sidereal_time += 8640184.812866 * t_lo
    sidereal_time += 3155760000.0 * t_lo
    sidereal_time += 8640184.812866 * t_hi
    sidereal_time += 3155760000.0 * t_hi

    return sidereal_time * pi / 43200.0
end

end # module SiderealTime
