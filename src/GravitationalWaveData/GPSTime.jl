module GPSTime

export GPS_EPOCH, gps_to_utc, utc_to_gps, leap_seconds

using Dates

# ================================================================================================
# GPS <-> UTC conversion.
#
# GPS time is a continuous count of SI seconds since its epoch (1980-01-06T00:00:00 UTC) with NO
# leap seconds inserted after that epoch -- it never stops, skips, or repeats a second. UTC, by
# contrast, has had a leap second inserted periodically (always at 23:59:60 on 30 June or 31
# December) to keep it within 0.9s of UT1 (mean solar time). This means GPS - UTC grows by exactly
# 1 second at each leap-second insertion, and any correct GPS<->UTC conversion needs the actual
# historical (and, for future dates, currently-known) table of when those insertions happened --
# NOT a single hardcoded constant, which silently goes stale the next time IERS announces one.
#
# TAI (International Atomic Time) is the usual bridge: TAI - UTC is exactly the leap-second count
# maintained by IERS Bulletin C, and GPS - TAI = -19s is a FIXED constant for all time (GPS's
# epoch was defined so that GPS = TAI - 19s at that instant, and GPS never adjusts afterward).
# So GPS - UTC = (TAI - UTC) - 19.
#
# `LEAP_SECOND_TABLE` below lists, for each UTC instant at which a new leap second took effect,
# the cumulative (TAI - UTC) offset from that instant onward. Source: IERS Bulletin C historical
# leap-second announcements (https://hpiers.obspm.fr/iers/bul/bulc/bulletinc.dat), cross-checked
# against the USNO/NIST published TAI-UTC table
# (https://data.iana.org/time-zones/data/leap-seconds.list carries the same dates in NTP epoch
# form). Vintage: table current as of Bulletin C 71 (July 2025) confirming NO leap second is to be
# inserted at the end of December 2025 -- the most recent actual insertion remains 2016-12-31
# (effective 2017-01-01, TAI-UTC=37s, i.e. GPS-UTC=18s), unchanged since. If IERS announces a new
# leap second (Bulletin C is issued twice yearly, ~5-6 months ahead of each possible insertion
# date), extend this table with the new (effective_utc_date, new_TAI_minus_UTC) entry -- do not
# just bump a constant elsewhere, since `gps_to_utc`/`utc_to_gps` both derive everything from this
# one table.
# ================================================================================================

"""
    GPS_EPOCH :: DateTime

GPS time's epoch: 1980-01-06T00:00:00 UTC (the start of GPS week 0). By definition, GPS time equals
UTC at this instant (`GPS - UTC = 0`).
"""
const GPS_EPOCH = DateTime(1980, 1, 6, 0, 0, 0)

# GPS - TAI is a fixed constant (never adjusted after GPS's epoch was defined).
const _GPS_MINUS_TAI = -19

# (utc_effective_date, TAI-UTC offset effective from that UTC instant onward). The first row is
# the anchor at GPS_EPOCH itself (TAI-UTC=19 there, consistent with GPS-UTC=0: 19 + (-19) = 0).
# Every subsequent row is a real historical leap-second insertion (see header comment for source).
const LEAP_SECOND_TABLE = [
    (DateTime(1980, 1, 6), 19),
    (DateTime(1981, 7, 1), 20),
    (DateTime(1982, 7, 1), 21),
    (DateTime(1983, 7, 1), 22),
    (DateTime(1985, 7, 1), 23),
    (DateTime(1988, 1, 1), 24),
    (DateTime(1990, 1, 1), 25),
    (DateTime(1991, 1, 1), 26),
    (DateTime(1992, 7, 1), 27),
    (DateTime(1993, 7, 1), 28),
    (DateTime(1994, 7, 1), 29),
    (DateTime(1996, 1, 1), 30),
    (DateTime(1997, 7, 1), 31),
    (DateTime(1999, 1, 1), 32),
    (DateTime(2006, 1, 1), 33),
    (DateTime(2009, 1, 1), 34),
    (DateTime(2012, 7, 1), 35),
    (DateTime(2015, 7, 1), 36),
    (DateTime(2017, 1, 1), 37),
]

# GPS-seconds-since-epoch at which each table row's offset takes effect, precomputed once. A leap
# second insertion is a single physical instant; in GPS-elapsed-seconds terms that instant is at
# `naive_calendar_diff(effective_date, GPS_EPOCH) + new_offset` (the "+new_offset" accounts for the
# fact that, by the time UTC's calendar reads `effective_date`, GPS has already ticked through the
# extra inserted second that UTC's calendar arithmetic alone doesn't see).
const _GPS_BREAKPOINTS = [
    (Dates.value(d - GPS_EPOCH) / 1000 + (off + _GPS_MINUS_TAI), off + _GPS_MINUS_TAI)
    for (d, off) in LEAP_SECOND_TABLE
]

"""
    leap_seconds(gps_or_utc) -> Int

`GPS - UTC` in whole seconds, as of the given instant (either a GPS time in seconds, `Real`, or a
UTC `DateTime`). This is `0` at `GPS_EPOCH` and increases by `1` at each historical leap-second
insertion (see `LEAP_SECOND_TABLE`).
"""
function leap_seconds(gps::Real)
    off = first(_GPS_BREAKPOINTS)[2]
    for (bp_gps, o) in _GPS_BREAKPOINTS
        bp_gps <= gps || break
        off = o
    end
    return off
end

function leap_seconds(t::DateTime)
    off = LEAP_SECOND_TABLE[1][2] + _GPS_MINUS_TAI
    for (d, tai_minus_utc) in LEAP_SECOND_TABLE
        d <= t || break
        off = tai_minus_utc + _GPS_MINUS_TAI
    end
    return off
end

"""
    utc_to_gps(t::DateTime) -> Float64

Convert a UTC date/time to GPS seconds (elapsed SI seconds since `GPS_EPOCH`), leap-second aware
via `LEAP_SECOND_TABLE`. Sub-millisecond precision is not representable (`DateTime`'s own
resolution limit); this is not a concern for GW analysis timescales.

Caveat shared by essentially every such converter: `DateTime` cannot represent the leap second
itself (`23:59:60`) since Julia's `Dates` module has no notion of leap seconds -- only ordinary
instants strictly before or after an insertion are meaningful inputs here.
"""
function utc_to_gps(t::DateTime)
    naive_seconds = Dates.value(t - GPS_EPOCH) / 1000
    return naive_seconds + leap_seconds(t)
end

"""
    gps_to_utc(gps::Real) -> DateTime

Convert GPS seconds (elapsed SI seconds since `GPS_EPOCH`) to a UTC date/time, leap-second aware
via `LEAP_SECOND_TABLE`. Rounds to the nearest millisecond (`DateTime`'s resolution limit).
"""
function gps_to_utc(gps::Real)
    utc_seconds = gps - leap_seconds(gps)
    return GPS_EPOCH + Millisecond(round(Int, utc_seconds * 1000))
end

end
