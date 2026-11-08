module Posteriors

export fetch_posterior_samples, list_posterior_analyses, load_posterior_samples

using HDF5
using ..GWOSC: GWOSC_DEFAULT_HOST, gravitationalwavedata_cache_dir, _gwosc_fetch_json, _gwosc_download,
               _gwosc_event_candidates

# ================================================================================================
# Published parameter-estimation (PE) posterior samples. GWOSC's per-event JSON (the `jsonurl` of
# each catalog entry, see `GWOSC.jl`) lists every analysis of the event under `"parameters"`; PE
# analyses have `pipeline_type == "pe"`, a `data_url` pointing at the data release, and one of them
# is flagged `is_preferred`. For GWTC-2.1 and GWTC-3 the data release is a single PESummary HDF5
# file on Zenodo per event (e.g. `IGWN-GWTC2p1-v2-GW150914_095045_PEDataRelease_mixed_cosmo.h5`)
# with one group per analysis (`C01:IMRPhenomXPHM`, `C01:SEOBNRv4PHM`, `C01:Mixed`, ...), each
# holding a compound `posterior_samples` dataset. Older catalogs (GWTC-1, GWTC-2) publish other
# formats (DCC-hosted `.hdf5`/`.tar`); `fetch_posterior_samples` only selects HDF5 releases, and the
# reader below only handles the PESummary layout.
# ================================================================================================

_is_hdf5_url(url::AbstractString) =
    occursin(r"\.(h5|hdf5)(/content)?$", split(url, '?')[1])

"""
    fetch_posterior_samples(event; version=nothing, cache_dir=gravitationalwavedata_cache_dir("posteriors"),
                            host="https://gwosc.org", force=false) -> String

Download the preferred published posterior-sample release for `event` (e.g. `"GW150914"`, matched as
in [`fetch_open_data`](@ref GravitationalWaveData.GWOSC.fetch_open_data)'s event method) and return the local file path.
GWOSC catalog versions are tried highest-first, and the first one whose preferred PE analysis has an
HDF5 data release is used, so recent events resolve to their GWTC-2.1/GWTC-3 PESummary files.
`version` restricts the search to one GWOSC event version.

Files are cached under `cache_dir` by their upstream file name; a cached file is returned without
re-downloading unless `force=true`. Read the result with [`list_posterior_analyses`](@ref) and
[`load_posterior_samples`](@ref).
"""
function fetch_posterior_samples(event::AbstractString;
                                 version::Union{Integer,Nothing}=nothing,
                                 cache_dir::AbstractString=gravitationalwavedata_cache_dir("posteriors"),
                                 host::AbstractString=GWOSC_DEFAULT_HOST, force::Bool=false)
    tried_urls = String[]
    for (_, meta) in _gwosc_event_candidates(event, version, host)
        jsonurl = meta["jsonurl"]
        push!(tried_urls, jsonurl)
        event_json = _gwosc_fetch_json(jsonurl)
        event_entry = first(values(get(event_json, "events", Dict{String,Any}())))
        params = values(get(event_entry, "parameters", Dict{String,Any}()))
        pe = [p for p in params if get(p, "pipeline_type", "") == "pe" &&
              _is_hdf5_url(something(get(p, "data_url", ""), ""))]
        isempty(pe) && continue
        preferred = filter(p -> get(p, "is_preferred", false) === true, pe)
        chosen = isempty(preferred) ? first(pe) : first(preferred)
        return _gwosc_download(chosen["data_url"], cache_dir; force=force)
    end
    throw(ArgumentError("No HDF5 posterior-sample release found for event $(repr(event)) " *
                         "(checked catalog release(s): $(join(tried_urls, ", ")))"))
end

"""
    list_posterior_analyses(path) -> Vector{String}

Names of the analyses in a PESummary posterior file (the top-level groups that contain a
`posterior_samples` dataset), e.g. `["C01:IMRPhenomXPHM", "C01:Mixed", "C01:SEOBNRv4PHM"]`.
"""
function list_posterior_analyses(path::AbstractString)
    return h5open(path, "r") do file
        sort([k for k in keys(file)
              if file[k] isa HDF5.Group && haskey(file[k], "posterior_samples")])
    end
end

"""
    load_posterior_samples(path, analysis=nothing) -> Dict{String,Vector{Float64}}

Read one analysis's posterior samples from a PESummary file as a column dictionary, parameter
name => samples (e.g. `"chirp_mass"`, `"mass_ratio"`, `"luminosity_distance"`, `"ra"`, `"dec"`,
...; names follow the bilby/PESummary conventions of the release). `analysis` defaults to
`"C01:Mixed"` when present (the combined-waveform result), otherwise the only analysis in the
file; see [`list_posterior_analyses`](@ref).
"""
function load_posterior_samples(path::AbstractString,
                                analysis::Union{AbstractString,Nothing}=nothing)
    analyses = list_posterior_analyses(path)
    isempty(analyses) &&
        throw(ArgumentError("$path has no PESummary `<analysis>/posterior_samples` datasets"))
    if analysis === nothing
        analysis = "C01:Mixed" in analyses ? "C01:Mixed" :
                   length(analyses) == 1 ? only(analyses) :
                   throw(ArgumentError("$path has several analyses $(analyses); pass one " *
                                        "explicitly"))
    elseif !(analysis in analyses)
        throw(ArgumentError("Analysis $(repr(analysis)) not in $path (available: $(analyses))"))
    end
    rows = h5open(file -> read(file["$analysis/posterior_samples"]), path, "r")
    isempty(rows) && return Dict{String,Vector{Float64}}()
    return Dict(String(name) => Float64[getfield(r, name) for r in rows]
                for name in fieldnames(eltype(rows)))
end

end
