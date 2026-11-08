module GWOSC

export load_gwosc_strain, condition_segment, fetch_open_data, gwdata_cache_dir

using HDF5
using DSP: DSP
using Downloads: Downloads
using JSON: JSON
using ..Types: TimeSeries

# ================================================================================================
# GWOSC/LOSC strain file I/O. File layout used: `/strain/Strain` (Float64 dataset with
# `Xspacing`/`Xstart` attributes), `/meta/GPSstart`, `/quality/simple/DQmask` (per-second
# bitmask, `Xspacing=1.0`).
# ================================================================================================

"""
    load_gwosc_strain(path) -> (strain::TimeSeries{Float64}, dq_mask::TimeSeries)

Read a standard GWOSC/LOSC HDF5 strain file, returning two `TimeSeries`:

  - `strain`: the `/strain/Strain` dataset, with `t0`/`dt` taken from that dataset's own
    `Xstart`/`Xspacing` attributes (not `/meta/GPSstart` directly -- both agree in every file
    checked so far, but the dataset's own attributes are the authoritative source for the sample
    grid).
  - `dq_mask`: the per-second data-quality bitmask `/quality/simple/DQmask`, as its own
    `TimeSeries` anchored to the same GPS epoch as `strain` but with its own (typically 1.0s)
    sample spacing read from the file rather than assumed.

The DQ mask is returned as a second `TimeSeries` rather than bundled as metadata on `strain`, since
it is itself a uniformly-sampled GPS-anchored series (just at a coarser rate) and gets
`times()`/`crop()` for free this way -- e.g. `crop(dq_mask, t0, t1)` to pull the DQ bits covering an
analysis segment.
"""
function load_gwosc_strain(path::AbstractString)
    strain_data, dt, gps_start, dq_data, dq_dt = h5open(path, "r") do file
        strain_ds = file["strain/Strain"]
        strain_data = read(strain_ds)
        dt = read_attribute(strain_ds, "Xspacing")
        gps_start = read_attribute(strain_ds, "Xstart")
        dq_ds = file["quality/simple/DQmask"]
        dq_data = read(dq_ds)
        dq_dt = read_attribute(dq_ds, "Xspacing")
        return strain_data, dt, gps_start, dq_data, dq_dt
    end
    strain = TimeSeries(strain_data, gps_start, dt; name="strain")
    dq_mask = TimeSeries(dq_data, gps_start, dq_dt; name="DQmask")
    return strain, dq_mask
end

# ================================================================================================
# Trigger-relative segmenting + Tukey windowing (`DSP.Windows.tukey`, symmetric convention:
# `padding=0, zerophase=false`).
# ================================================================================================

"""
    condition_segment(ts::TimeSeries; trigger_time, duration=4.0, post_trigger=2.0, alpha=0.2)
        -> TimeSeries

Slice a trigger-relative analysis segment out of a full strain `TimeSeries` and apply a Tukey
window (`DSP.Windows.tukey`, `alpha` default matching common LIGO/bilby conditioning practice):
`duration` total seconds, ending `post_trigger` seconds after
`trigger_time` (so the segment starts at `trigger_time - (duration - post_trigger)`).

Only the DATA segment is windowed here -- the frequency-domain waveform TEMPLATE evaluated
downstream is analytic and is not (and should not be) windowed; this asymmetry (tapered data,
un-tapered analytic template) is standard practice in matched-filter GW parameter estimation, not
an oversight.
"""
function condition_segment(ts::TimeSeries; trigger_time::Real, duration::Real=4.0,
                            post_trigger::Real=2.0, alpha::Real=0.2)
    dt = ts.dt
    seg_start_time = trigger_time - (duration - post_trigger)
    i0 = round(Int, (seg_start_time - ts.t0) / dt) + 1
    nseg = round(Int, duration / dt)
    i1 = i0 + nseg - 1
    (i0 >= 1 && i1 <= length(ts)) ||
        throw(ArgumentError("Requested segment [$seg_start_time, $(seg_start_time + duration)) " *
                             "falls outside the strain TimeSeries' coverage " *
                             "[$(ts.t0), $(ts.t0 + length(ts) * dt))"))

    windowed = ts.data[i0:i1] .* DSP.Windows.tukey(nseg, alpha)
    seg_t0 = ts.t0 + (i0 - 1) * dt
    return TimeSeries(windowed, seg_t0, dt; name=ts.name)
end

# ================================================================================================
# `fetch_open_data`: download an official GWOSC strain HDF5 file over HTTPS via stdlib
# `Downloads.jl`, returning a local path for `load_gwosc_strain`. Endpoints used (all under
# https://gwosc.org, GWOSC's public API host, documented at https://gwosc.org/api/v1/docs/):
#
#   - `/archive/{gpsstart}/{gpsend}/json/` -- lists which named science-run datasets (e.g. "O1",
#     "O1_16KHZ") and/or catalog events overlap a GPS interval.
#   - `/archive/links/{run}/{detector}/{gpsstart}/{gpsend}/json/` -- for one named run + detector,
#     lists the actual strain files (`"strain"`: array of `{GPSstart, detector, sampling_rate,
#     duration, format, url, ...}`) overlapping that interval, in every offered format
#     (hdf5/gwf/txt.gz).
#   - `/eventapi/json/allevents/` -- catalog of every confirmed-event dataset (keyed
#     `"{commonName}-v{version}"`), each with a `commonName`, `version`, and a `jsonurl` pointing at
#     that event/version's own detailed JSON.
#   - an event's own `jsonurl` (e.g. `/eventapi/json/O3_Discovery_Papers/GW190412/v1`) -- same
#     `"strain"` array shape as the run-links endpoint above, but for that event's short
#     (typically 32s, up to 16 kHz) data-release segments plus the full science-run file.
# ================================================================================================

const GWOSC_DEFAULT_HOST = "https://gwosc.org"

# --- HTTP fetch helpers ----------------------------------------------------------------------

"""
    _gwosc_fetch_json(url) -> parsed JSON

`Downloads.download` `url` into memory and parse it as JSON. Network/HTTP failures (GWOSC
unreachable, DNS failure, 404, ...) are converted into a clear `ErrorException` rather than left
as a raw `Downloads.RequestError`/hang, per this package's convention of failing loudly.
"""
function _gwosc_fetch_json(url::AbstractString)
    io = IOBuffer()
    try
        Downloads.download(url, io)
    catch e
        throw(ErrorException("Failed to reach GWOSC at $url ($e). This means either GWOSC " *
                              "(https://gwosc.org) is down/unreachable, or this environment has no " *
                              "outbound HTTPS access -- not a bug in the request itself."))
    end
    return JSON.parse(String(take!(io)); dicttype=Dict{String,Any})
end

"""
    _url_filename(url) -> String

The file name a download URL refers to: its last path component, or the one before it for
Zenodo-style `.../files/<name>/content` URLs.
"""
function _url_filename(url::AbstractString)
    parts = split(rstrip(split(url, '?')[1], '/'), '/')
    zenodo = parts[end] == "content" && length(parts) > 1
    return String(zenodo ? parts[end - 1] : parts[end])
end

"""
    _gwosc_download(url, cache_dir; force=false) -> path

Download `url` (a GWOSC strain file or other data-release URL) into `cache_dir`, named after the
URL's own filename (see `_url_filename`; GWOSC and catalog filenames already encode
detector/GPS-start/duration/rate/version, e.g. `H-H1_LOSC_4_V1-1126256640-4096.hdf5`, so this is a
natural content-addressed cache key). If a same-named file already exists in `cache_dir` and
`force=false`, the network is skipped entirely and the cached path is returned. Downloads to a
`.part` sibling file first and `mv`s it into place on success, so a failed/interrupted download
never leaves a corrupt file at the cache path.
"""
function _gwosc_download(url::AbstractString, cache_dir::AbstractString; force::Bool=false)
    mkpath(cache_dir)
    dest = joinpath(cache_dir, _url_filename(url))
    (!force && isfile(dest) && filesize(dest) > 0) && return dest
    tmp = dest * ".part"
    try
        Downloads.download(url, tmp)
    catch e
        isfile(tmp) && rm(tmp; force=true)
        throw(ErrorException("Failed to download $url ($e). Either the host is down/unreachable, " *
                              "or this environment has no outbound HTTPS access."))
    end
    mv(tmp, dest; force=true)
    return dest
end

"""
    gwdata_cache_dir(subdirs...) -> String

Root of GWData.jl's download cache, optionally joined with `subdirs`. Resolved in order:

  1. `ENV["GWDATA_CACHE_DIR"]`, if set;
  2. `\$XDG_CACHE_HOME/gwdata`, if `XDG_CACHE_HOME` is set;
  3. `~/.cache/gwdata`.

The directory is a plain tree of files under their upstream names, so other tools can read it
directly. Subdirectories: `strain/` holds GWOSC strain files downloaded by
[`fetch_open_data`](@ref), and `posteriors/` holds parameter-estimation releases downloaded by
`fetch_posterior_samples`. The path is not created by this function.
"""
function gwdata_cache_dir(subdirs::AbstractString...)
    root = if haskey(ENV, "GWDATA_CACHE_DIR")
        ENV["GWDATA_CACHE_DIR"]
    elseif haskey(ENV, "XDG_CACHE_HOME")
        joinpath(ENV["XDG_CACHE_HOME"], "gwdata")
    else
        joinpath(homedir(), ".cache", "gwdata")
    end
    return joinpath(root, subdirs...)
end

_default_gwosc_cache_dir() = gwdata_cache_dir("strain")

"""
    _gwosc_event_candidates(event, version, host) -> Vector{Pair{String,Any}}

GWOSC catalog entries (`key => metadata`, metadata including `jsonurl`) for `event`, matched
against the event `commonName`, the dataset key, or the `GWYYMMDD` prefix of the common name,
preferring exact `commonName` matches. Restricted to `version` if given, and ordered
highest-version-first. Throws `ArgumentError` if nothing matches.
"""
function _gwosc_event_candidates(event::AbstractString, version::Union{Integer,Nothing},
                                 host::AbstractString)
    allevents = _gwosc_fetch_json("$host/eventapi/json/allevents/")
    events_dict = get(allevents, "events", Dict{String,Any}())

    matched = Pair{String,Any}[]
    for kv in events_dict
        key, meta = kv
        common = get(meta, "commonName", "")
        prefix = split(common, "_")[1]
        (event == key || event == common || event == prefix) && push!(matched, kv)
    end
    isempty(matched) &&
        throw(ArgumentError("No GWOSC event dataset matches $(repr(event)) " *
                             "(queried $host/eventapi/json/allevents/)"))
    exact = filter(kv -> get(kv[2], "commonName", "") == event, matched)
    pool = isempty(exact) ? matched : exact
    if version !== nothing
        pool = filter(kv -> get(kv[2], "version", nothing) == version, pool)
        isempty(pool) &&
            throw(ArgumentError("No GWOSC event dataset matches $(repr(event)) at version " *
                                 "$version"))
    end
    return sort(pool; by=kv -> get(kv[2], "version", 0), rev=true)
end

_matches_version(url::AbstractString, version::Integer) =
    occursin("_V$(version)-", url) || occursin("_R$(version)-", url)

"""
    fetch_open_data(detector, gps_start, gps_end; sample_rate=4096, version=nothing,
                     cache_dir=<user cache>, host="https://gwosc.org", force=false) -> String

Download the official GWOSC science-run strain HDF5 file whose coverage fully contains
`[gps_start, gps_end)` for `detector` (e.g. `"H1"`, `"L1"`, `"V1"`), returning a local file path
directly consumable by [`load_gwosc_strain`](@ref) -- this function only gets the bytes onto disk,
it doesn't parse them (that's already solved).

GWOSC's public science-run archive is split into fixed-length (typically 4096s) segment files, so
this only succeeds when the requested range falls within a single such segment; it raises a clear
`ArgumentError` (not a partial/silently-wrong download) if no single hdf5 file at `sample_rate`
fully covers the request. `version`, if given, additionally requires the chosen file's name to
carry that data-release version tag (`_V<version>-`/`_R<version>-`, GWOSC's own file-naming
convention).

Downloaded files are cached under `cache_dir` (default `gwdata_cache_dir("strain")`, i.e.
`~/.cache/gwdata/strain` unless `GWDATA_CACHE_DIR` or `XDG_CACHE_HOME` is set; see
[`gwdata_cache_dir`](@ref)) keyed by GWOSC's own filename; a repeat call for the same file re-does
the (cheap, small-JSON) locator query but skips re-downloading the (potentially large) strain file
itself, unless `force=true`. See this module's docstring above for the exact GWOSC endpoints
queried.
"""
function fetch_open_data(detector::AbstractString, gps_start::Real, gps_end::Real;
                          sample_rate::Real=4096, version::Union{Integer,Nothing}=nothing,
                          cache_dir::AbstractString=_default_gwosc_cache_dir(),
                          host::AbstractString=GWOSC_DEFAULT_HOST, force::Bool=false)
    gps_start < gps_end ||
        throw(ArgumentError("gps_start ($gps_start) must be < gps_end ($gps_end)"))
    gi0, gi1 = round(Int, gps_start), round(Int, gps_end)

    dataset_json = _gwosc_fetch_json("$host/archive/$gi0/$gi1/json/")
    runs = get(dataset_json, "runs", Dict{String,Any}())
    names = sort(collect(keys(runs)))
    non_history = filter(!=("history"), names)
    candidates = isempty(non_history) ? names : non_history
    isempty(candidates) &&
        throw(ArgumentError("GWOSC has no science-run dataset covering GPS [$gi0, $gi1) at all " *
                             "(queried $host/archive/$gi0/$gi1/json/)"))
    # Ordering heuristic only (every candidate is still verified against its actual reported
    # `sampling_rate` below): try runs whose name doesn't look like a different-sample-rate variant
    # of the requested one first (e.g. skip "O1_16KHZ" first when `sample_rate` is the 4096 Hz
    # default).
    wants_16k = sample_rate > 4096
    sort!(candidates; by=name -> (occursin("16KHZ", uppercase(name)) == wants_16k) ? 0 : 1)

    for run in candidates
        links_json = _gwosc_fetch_json("$host/archive/links/$run/$detector/$gi0/$gi1/json/")
        for entry in get(links_json, "strain", Any[])
            get(entry, "format", "") == "hdf5" || continue
            get(entry, "detector", "") == detector || continue
            round(Int, get(entry, "sampling_rate", -1)) == round(Int, sample_rate) || continue
            entry_start = round(Int, entry["GPSstart"])
            entry_end = entry_start + round(Int, entry["duration"])
            (entry_start <= gi0 && entry_end >= gi1) || continue
            url = entry["url"]
            (version === nothing || _matches_version(url, version)) || continue
            return _gwosc_download(url, cache_dir; force=force)
        end
    end
    throw(ArgumentError("No single GWOSC hdf5 file for detector $(repr(detector)) fully covers " *
                         "GPS [$gi0, $gi1) at $sample_rate Hz (checked dataset(s): " *
                         "$(join(candidates, ", "))). GWOSC science-run files are fixed-length " *
                         "segments -- request a range within one segment, or check the detector/" *
                         "GPS range/sample_rate/version against https://gwosc.org."))
end

"""
    fetch_open_data(event, detector; sample_rate=4096, duration=nothing, version=nothing,
                     cache_dir=<user cache>, host="https://gwosc.org", force=false) -> String

Download the official GWOSC event-data-release strain HDF5 file for `event` (e.g. `"GW190412"`,
matched against GWOSC's event `commonName`/dataset key/`GWYYMMDD` prefix, following the highest
available catalog version unless `version` is given) and `detector`, returning a local file path
directly consumable by [`load_gwosc_strain`](@ref).

Events typically offer several segment lengths (e.g. a short ~32s release alongside the full
science-run file); `duration` (seconds) selects one exactly, defaulting to the *shortest* available
match at `sample_rate` (the common case: a small smoke-test/quick-look segment rather than the full
4096s file). See this module's docstring above for the exact GWOSC endpoints queried.
"""
function fetch_open_data(event::AbstractString, detector::AbstractString;
                          sample_rate::Real=4096, duration::Union{Real,Nothing}=nothing,
                          version::Union{Integer,Nothing}=nothing,
                          cache_dir::AbstractString=_default_gwosc_cache_dir(),
                          host::AbstractString=GWOSC_DEFAULT_HOST, force::Bool=false)
    ordered = _gwosc_event_candidates(event, version, host)
    # `version` here is GWOSC's per-(event,catalog) release version, not a global ranking of data
    # quality -- e.g. GW190412 has versions from "O3_Discovery_Papers", "GWTC-2", and independent
    # reanalysis catalogs like "IAS-O3a", and not every catalog publishes a strain-file release
    # (some are parameter-estimation-only). So rather than committing to the single
    # highest-`version` candidate and failing if THAT catalog has no matching strain file, try
    # candidates highest-version-first and fall through to the next until one actually has a
    # matching strain entry.
    tried_urls = String[]
    for (_, meta) in ordered
        jsonurl = meta["jsonurl"]
        push!(tried_urls, jsonurl)
        event_json = _gwosc_fetch_json(jsonurl)
        event_entry = first(values(get(event_json, "events", Dict{String,Any}())))
        matches = filter(get(event_entry, "strain", Any[])) do e
            get(e, "format", "") == "hdf5" || return false
            get(e, "detector", "") == detector || return false
            round(Int, get(e, "sampling_rate", -1)) == round(Int, sample_rate) || return false
            duration === nothing || round(Int, get(e, "duration", -1)) == round(Int, duration) ||
                return false
            return true
        end
        isempty(matches) && continue
        chosen = duration === nothing ? matches[argmin([e["duration"] for e in matches])] :
                  matches[1]
        return _gwosc_download(chosen["url"], cache_dir; force=force)
    end
    throw(ArgumentError("No GWOSC hdf5 strain file found for event $(repr(event)), detector " *
                         "$(repr(detector)) at $sample_rate Hz" *
                         (duration === nothing ? "" : " and duration=$duration") *
                         " (checked catalog release(s): $(join(tried_urls, ", ")))"))
end

end
