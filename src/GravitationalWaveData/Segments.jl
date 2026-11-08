module Segments

export Segment, SegmentList, coalesce_segments, gaps, duration, overlaps

import ..Types: duration

# ================================================================================================
# Segment / SegmentList: half-open [start, stop) GPS-time intervals, for data-quality vetoes.
#
# Design follows the standard LIGO `segment`/`segmentlist` convention:
# a `SegmentList` is always kept SORTED and COALESCED (overlapping or touching segments merged into
# one) as a class invariant, maintained on construction rather than checked ad hoc by every
# operation. This is the same half-open-interval convention `TimeSeries`/`FrequencySeries` use for
# `crop` in `Types.jl`, for consistency across the package.
#
# Typical use: decode the per-second `dq_mask` `TimeSeries` returned by `load_gwosc_strain`
# bit-by-bit into per-flag `SegmentList`s of "good" intervals, then combine them via
# `intersect`/`union` across flags/detectors before applying them as a veto.
# ================================================================================================

"""
    Segment(start, stop)

A half-open GPS-time interval `[start, stop)`, `start < stop` required. Same half-open convention
`TimeSeries`/`FrequencySeries` use for `crop`.
"""
struct Segment
    start::Float64
    stop::Float64

    function Segment(start::Real, stop::Real)
        start < stop ||
            throw(ArgumentError("Segment requires start < stop, got start=$start, stop=$stop"))
        return new(Float64(start), Float64(stop))
    end
end

"""
    duration(s::Segment) -> Float64

Length of the interval, `s.stop - s.start`.
"""
duration(s::Segment) = s.stop - s.start

Base.in(t::Real, s::Segment) = s.start <= t < s.stop
Base.isless(a::Segment, b::Segment) = isless((a.start, a.stop), (b.start, b.stop))
Base.:(==)(a::Segment, b::Segment) = a.start == b.start && a.stop == b.stop

"""
    overlaps(a::Segment, b::Segment) -> Bool

Whether the two half-open intervals share any point. Merely touching (`a.stop == b.start`) does
NOT count as overlapping here (no shared point in a half-open interval) -- but touching segments
are still merged by `coalesce_segments`/`SegmentList`, which is the more permissive, canonicalizing
notion of "adjacent enough to combine" that segment-list libraries conventionally use.
"""
overlaps(a::Segment, b::Segment) = a.start < b.stop && b.start < a.stop

function Base.show(io::IO, s::Segment)
    print(io, "Segment($(s.start), $(s.stop))")
end

# ================================================================================================
# SegmentList
# ================================================================================================

"""
    SegmentList(segments)
    SegmentList()

A sorted, coalesced list of `Segment`s: constructing from any iterable of `Segment`s sorts by
start time and merges overlapping-or-touching segments into one, as a maintained invariant (so
every `SegmentList` a caller holds is always already in canonical form -- no separate "did you
remember to coalesce?" step).
"""
struct SegmentList
    segments::Vector{Segment}

    global _unsafe_segmentlist(segs::Vector{Segment}) = new(segs)
end

function SegmentList(segs)
    return _unsafe_segmentlist(_coalesce(collect(Segment, segs)))
end
SegmentList() = _unsafe_segmentlist(Segment[])

"""
    coalesce_segments(segments) -> SegmentList

Sort `segments` (any iterable of `Segment`) and merge overlapping-or-touching ones. Equivalent to
`SegmentList(segments)` -- provided as an explicitly-named entry point since "coalescing" is a
named operation callers reach for directly (e.g. after collecting raw candidate segments from a
DQ-bitmask scan).
"""
coalesce_segments(segs) = SegmentList(segs)

function _coalesce(segs::Vector{Segment})
    isempty(segs) && return Segment[]
    sorted = sort(segs)
    out = Segment[sorted[1]]
    for s in @view sorted[2:end]
        last = out[end]
        if s.start <= last.stop   # overlapping or touching -> merge
            out[end] = Segment(last.start, max(last.stop, s.stop))
        else
            push!(out, s)
        end
    end
    return out
end

Base.length(sl::SegmentList) = length(sl.segments)
Base.isempty(sl::SegmentList) = isempty(sl.segments)
Base.iterate(sl::SegmentList, state...) = iterate(sl.segments, state...)
Base.getindex(sl::SegmentList, i) = sl.segments[i]
Base.firstindex(sl::SegmentList) = firstindex(sl.segments)
Base.lastindex(sl::SegmentList) = lastindex(sl.segments)
Base.:(==)(a::SegmentList, b::SegmentList) = a.segments == b.segments
Base.in(t::Real, sl::SegmentList) = any(s -> t in s, sl.segments)

"""
    duration(sl::SegmentList) -> Float64

Total covered time, i.e. the sum of each (disjoint, by construction) segment's duration.
"""
duration(sl::SegmentList) = isempty(sl) ? 0.0 : sum(duration, sl.segments)

function Base.show(io::IO, sl::SegmentList)
    print(io, "SegmentList($(length(sl)) segment", length(sl) == 1 ? "" : "s",
          ", total duration=$(duration(sl))s)")
end

"""
    union(a::SegmentList, b::SegmentList) -> SegmentList

All time covered by either list, coalesced.
"""
function Base.union(a::SegmentList, b::SegmentList)
    return SegmentList(vcat(a.segments, b.segments))
end

"""
    intersect(a::SegmentList, b::SegmentList) -> SegmentList

Time covered by both lists. Implemented as a linear two-pointer sweep over each (sorted, disjoint)
list's segments, since both inputs are already in canonical form.
"""
function Base.intersect(a::SegmentList, b::SegmentList)
    out = Segment[]
    i, j = 1, 1
    asegs, bsegs = a.segments, b.segments
    while i <= length(asegs) && j <= length(bsegs)
        sa, sb = asegs[i], bsegs[j]
        lo = max(sa.start, sb.start)
        hi = min(sa.stop, sb.stop)
        lo < hi && push!(out, Segment(lo, hi))
        if sa.stop <= sb.stop
            i += 1
        else
            j += 1
        end
    end
    return _unsafe_segmentlist(out)   # already sorted & disjoint by construction
end

"""
    gaps(sl::SegmentList) -> SegmentList
    gaps(sl::SegmentList, bound::Segment) -> SegmentList

The complement of `sl`: time NOT covered by any of its segments.

Without `bound`, only the "internal" gaps strictly between `sl`'s first and last segment are
returned (there is no well-defined gap before the first or after the last segment when the list's
own extent is the only reference frame).

With an explicit `bound`, the complement is computed within `[bound.start, bound.stop)` -- e.g.
`gaps(good_data_segments, Segment(analysis_start, analysis_end))` to find the vetoed time within a
specific analysis window, including any leading/trailing gap against that window's edges. `sl`'s
segments are clipped to `bound` first (segments partially or fully outside `bound` are handled
correctly).
"""
function gaps(sl::SegmentList)
    segs = sl.segments
    out = Segment[]
    for k in 1:(length(segs) - 1)
        segs[k].stop < segs[k + 1].start && push!(out, Segment(segs[k].stop, segs[k + 1].start))
    end
    return _unsafe_segmentlist(out)
end

function gaps(sl::SegmentList, bound::Segment)
    lo, hi = bound.start, bound.stop
    clipped = Segment[]
    for s in sl.segments
        cs, ce = max(s.start, lo), min(s.stop, hi)
        cs < ce && push!(clipped, Segment(cs, ce))
    end

    out = Segment[]
    cursor = lo
    for s in clipped
        cursor < s.start && push!(out, Segment(cursor, s.start))
        cursor = s.stop
    end
    cursor < hi && push!(out, Segment(cursor, hi))
    return _unsafe_segmentlist(out)
end

end
