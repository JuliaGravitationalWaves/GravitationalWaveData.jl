# # GPS Time & Segments
#
# Gravitational-wave pipelines timestamp everything in **GPS seconds** rather than UTC, and track
# "when was the detector taking good data" with half-open **segments**, not individual timestamps.
# This tutorial covers both corners of `GravitationalWaveData.jl`: [`GPSTime`](@ref GravitationalWaveData.gps_to_utc) for
# converting between the two time systems, and [`Segments`](@ref GravitationalWaveData.Segment) for building and
# combining data-availability intervals.
#
# `GravitationalWaveData.jl` doesn't fetch data-quality segments over the network, so the example builds a
# small, hand-crafted "instrument locked" scenario instead: it runs instantly, offline, and shows
# the segment algebra (`union`, `intersect`, `gaps`) used to find, for example, double-coincident
# observing time.

using GravitationalWaveData
using Dates

# ## Part 1: GPS time <-> UTC
#
# GPS time is a continuous SI-second count from the GPS epoch (1980-01-06T00:00:00 UTC) with **no**
# leap seconds inserted after that epoch, unlike UTC. `GravitationalWaveData.jl` converts between the two via a
# real historical leap-second table (see `docs/CONVENTIONS.md` and `src/GravitationalWaveData/GPSTime.jl` for the
# full rationale) rather than a single hardcoded offset, since that offset has changed 18 times
# since 1980 and will change again whenever IERS announces a new leap second.
#
# The GPS epoch itself is the simplest round-trip: `GPS - UTC = 0` there by definition.

GPS_EPOCH

#-

utc_to_gps(GPS_EPOCH)

#-

gps_to_utc(0.0)

# `leap_seconds` reports the current `GPS - UTC` offset (in whole seconds) as of a given GPS time
# or UTC `DateTime`. It is `0` at the epoch and has stepped up by 1 at each of the 18 leap seconds
# inserted since:

leap_seconds(GPS_EPOCH)

# Two real, published GW event trigger times make good round-trip examples, since they straddle
# one of those leap-second insertions (2017-01-01, when `GPS - UTC` stepped from 17s to 18s) and so
# exercise two different rows of the underlying table:

gw150914_gps = 1_126_259_462.4   # GW150914, before the 2017-01-01 leap second
leap_seconds(gw150914_gps)

#-

gw150914_utc = gps_to_utc(gw150914_gps)

#-

utc_to_gps(gw150914_utc) ≈ gw150914_gps

#-

gw170817_gps = 1_187_008_882.4   # GW170817, after the 2017-01-01 leap second
leap_seconds(gw170817_gps)

#-

gps_to_utc(gw170817_gps)

# Both round-trip exactly (up to `DateTime`'s millisecond resolution limit — GPS time itself can
# carry sub-millisecond precision, `DateTime` cannot). The one input `gps_to_utc`/`utc_to_gps` can't
# round-trip meaningfully is an instant that falls exactly on a leap-second insertion itself (the
# inserted leap second `23:59:60` has no `DateTime` representation at all, since Julia's `Dates`
# module has no notion of leap seconds) — not a concern for real analysis GPS times, which are
# essentially never chosen to land on that exact boundary.

# ## Part 2: Segment and SegmentList
#
# `Segment(start, stop)` is a half-open `[start, stop)` GPS-time interval — the same half-open
# convention `TimeSeries`/`FrequencySeries` cropping uses elsewhere in this package, so a boundary
# sample is never double-counted depending on which API touched it. A `SegmentList` is a sorted,
# automatically-coalesced collection of `Segment`s: overlapping *or touching* segments merge into
# one on construction, so every `SegmentList` you hold is already in canonical form.
#
# Suppose H1 and L1 each log a handful of "interferometer locked" intervals (GPS-second offsets
# from an arbitrary reference `t0`, chosen here so the window brackets GW150914's own trigger time)
# over a two-hour stretch:

t0 = gw150914_gps - 3600.0   # window starts 1 hour before the GW150914 trigger
trigger_offset = gw150914_gps - t0
trigger_offset

#-

h1_locked = SegmentList([
    Segment(0.0, 1200.0),      # locked for the first 20 minutes
    Segment(1500.0, 3600.0),   # relocks after a 5-minute break...
    Segment(3600.0, 5400.0),   # ...and holds straight through (touches the previous segment)
    Segment(6000.0, 7200.0),
])

# Note that the middle two segments above were given as separate, *touching* intervals
# (`[1500,3600)` then `[3600,5400)`) but `SegmentList` merged them into one on construction, exactly
# as `docs/CONVENTIONS.md` describes — touching segments coalesce even though `overlaps` (below)
# would say they don't overlap:

h1_locked

#-

l1_locked = SegmentList([
    Segment(600.0, 2400.0),
    Segment(3000.0, 4800.0),
    Segment(5000.0, 7200.0),
])

l1_locked

# `overlaps` answers a sharper question than "do these segments' lists cover common time" — it
# tests two individual `Segment`s for a genuinely shared point, which merely *touching* segments do
# not have:

overlaps(Segment(0.0, 10.0), Segment(5.0, 15.0))   # true: share (5, 10)

#-

overlaps(Segment(0.0, 10.0), Segment(10.0, 20.0))   # false: touch at 10, but half-open excludes it

# ### Combining lists: union, intersect, duration
#
# `union` gives every second covered by either instrument — useful for "when was *any* data being
# taken" — and `intersect` gives only the double-coincident time both interferometers need for a
# confident detection:

either_locked = union(h1_locked, l1_locked)

#-

both_locked = intersect(h1_locked, l1_locked)

# `duration` totals the covered time of a `SegmentList` (or the length of a single `Segment`), so we
# can quantify the double-coincident duty cycle over this two-hour window directly:

analysis_window = Segment(0.0, 7200.0)
duration(both_locked) / duration(analysis_window)

# Was GW150914's own trigger time (`trigger_offset` above) inside a double-coincident stretch?
# `Segment`/`SegmentList` support `in` for exactly this query:

trigger_offset in both_locked

# ### Gaps
#
# `gaps` returns the complement of a `SegmentList` — the vetoed/no-data time. Without a `bound`,
# only the *internal* gaps (strictly between the first and last segment) are meaningful:

gaps(h1_locked)

# Passing an explicit `bound::Segment` computes the complement within that specific window instead,
# including any leading/trailing gap against the window's own edges — e.g. against the full
# two-hour `analysis_window` defined above, H1's `[0.0, 7200.0)` list happens to run edge-to-edge
# already, so widen the bound slightly to see leading/trailing gaps appear:

gaps(h1_locked, Segment(-600.0, 7800.0))

# The same query on the double-coincident list shows exactly which stretches of the analysis window
# are *not* usable for a joint H1+L1 search — the complement of `both_locked` within
# `analysis_window`:

gaps(both_locked, analysis_window)

# ### Tying the two halves together
#
# Every offset above is a plain GPS-second count relative to `t0`; converting a segment boundary
# back to an absolute, human-readable UTC timestamp is just `gps_to_utc(t0 + offset)`, using exactly
# the leap-second-aware conversion from Part 1. For instance, the start of the first double-locked
# stretch:

first_both = both_locked[1]
gps_to_utc(t0 + first_both.start)
