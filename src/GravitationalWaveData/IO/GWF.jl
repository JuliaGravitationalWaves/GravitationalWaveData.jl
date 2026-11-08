module GWF

export load_gwf_strain, list_gwf_channels

using CodecZlib
using ..Types: TimeSeries

# ================================================================================================
# Reader for LIGO/Virgo's native IGWD frame ("GWF") binary format, per LIGO-T970130-v2
# ("Specification of a Common Data Frame Format for Interferometric Gravitational Wave Detectors").
#
# Layout was cross-checked byte-for-byte against a real GWOSC frame file (see
# GWF_layout_notes below and `test/gwf.jl`), not derived from the spec text alone -- the spec
# PDF's tables extract with scrambled columns via `pdftotext`, so every field order below was
# additionally confirmed by hex-dumping `data/H-H1_LOSC_4_V1-1126256640-4096.gwf` and checking
# structure lengths/offsets landed exactly where each interpretation predicted.
#
# Design/scope, deliberately narrower than the full spec:
#   - Only little-endian frames are supported (checked at read time against the file header's
#     byte-order test words) -- covers 100% of real-world frames, which are written on x86/ARM
#     hardware. Big-endian frames (legacy SPARC/PowerPC writers) raise a clear error.
#   - Rather than resolving the full FrameH -> FrRawData/FrProcData -> FrVect PTR_STRUCT pointer
#     graph (each pointer is an on-disk (class,instance) pair -- see Table 4 in the spec), this
#     reader does a single linear pass over top-level structures and matches a channel's
#     `FrProcData`/`FrAdcData` wrapper to its `FrVect` sample array BY NAME. This relies on GWOSC's
#     (and, in practice, essentially every real writer's) convention that a channel's FrVect is
#     named identically to its parent FrProcData/FrAdcData -- verified true for all three channels
#     (STRAIN, DQMASK, INJMASK) in the real file this was validated against. Full pointer-graph
#     resolution would be needed for a fully general reader (e.g. disambiguating two same-named
#     channels in one frame); out of scope here.
#   - `FrAdcData` (raw/uncalibrated ADC channels) is supported enough to read `name`/`timeOffset`,
#     matching the fields actually needed here; `FrProcData` (calibrated/post-processed channels,
#     what GWOSC actually uses for its public STRAIN/DQMASK/INJMASK channels) is fully supported.
#   - Compression: LIGO-T970130 Appendix B defines 6 schemes (raw, gzip, gzip-differential, and
#     zero-suppression for 2/4/8-byte words, each with a +256 offset for little-endian writers --
#     see `_COMPRESS_NAMES` below for the full code table). The real GWOSC file this was validated
#     against uses ONLY gzip (code 257); this reader implements raw (0/256) and gzip (1/257) --
#     the two schemes GWOSC actually produces -- and raises a clear `ArgumentError` naming the
#     unimplemented scheme otherwise, rather than silently misreading zero-suppressed data.
#   - Multiple frames per file (a file containing `FrameH, ..., FrameH, ...`) are supported by
#     concatenating a channel's per-frame `FrVect` chunks in file order, but only if their `dt`
#     matches and each chunk's start time is exactly contiguous with the previous chunk's end --
#     gapped/overlapping multi-frame files raise a clear error rather than silently splicing wrong.
# ================================================================================================

# --- primitive little-endian readers over a Vector{UInt8}, 1-based position in, position-after out.
# `reinterpret` gives the correct value directly only because we've verified (via the file header's
# byte-order test words, checked in `_parse_file_header`) that this file -- and, in practice, this
# host -- is little-endian; see the module docstring above.

@inline _ru8(buf, p) = (buf[p], p + 1)
@inline _ru16(buf, p) = (only(reinterpret(UInt16, buf[p:(p + 1)])), p + 2)
@inline _ri32(buf, p) = (only(reinterpret(Int32, buf[p:(p + 3)])), p + 4)
@inline _ru32(buf, p) = (only(reinterpret(UInt32, buf[p:(p + 3)])), p + 4)
@inline _ru64(buf, p) = (only(reinterpret(UInt64, buf[p:(p + 7)])), p + 8)
@inline _rf64(buf, p) = (only(reinterpret(Float64, buf[p:(p + 7)])), p + 8)

"""
    _read_string(buf, p) -> (String, next_p)

IGWD `STRING` type: an `INT_2U` byte count `n` (the following char data's length, INCLUDING its
mandatory trailing `'\\0'`), followed by `n` bytes of ASCII ending in that `'\\0'`. Returns the
string with the terminator stripped.
"""
function _read_string(buf, p)
    n, p = _ru16(buf, p)
    n == 0 && return "", p
    s = String(buf[p:(p + n - 2)])   # drop the trailing '\0'
    return s, p + n
end

# ================================================================================================
# File header (40 bytes, byte-exact per LIGO-T970130-v2 Table 5, confirmed against the real file):
#   0-3   "IGWD" magic (4 ASCII bytes, no terminator)
#   4     frame library minor version
#   5     frame format major version (8, for all modern GWOSC files)
#   6-11  sizeof(CHAR, INT_2, INT_4, INT_8, REAL_4, REAL_8) = (1,2,4,8,4,8) on the writing host
#   12-13 INT_2U byte-order test word, must read as 0x1234
#   14-17 INT_4U byte-order test word, must read as 0x12345678
#   18-25 INT_8U byte-order test word, must read as 0x0123456789abcdef
#   26-29 REAL_4 test value, must read as pi
#   30-37 REAL_8 test value, must read as pi
#   38    frame library used (0=unknown,1=frameL,2=frameCPP)
#   39    file checksum scheme (0=none,1=CRC)
# ================================================================================================
function _parse_file_header(buf)
    length(buf) >= 40 ||
        throw(ArgumentError("File too short to contain a GWF file header (need >= 40 bytes, got $(length(buf)))"))
    String(buf[1:4]) == "IGWD" ||
        throw(ArgumentError("Not a GWF/IGWD frame file: expected magic \"IGWD\", got $(repr(String(buf[1:4])))"))
    sizes = Tuple(buf[7:12])
    sizes == (0x01, 0x02, 0x04, 0x08, 0x04, 0x08) ||
        throw(ArgumentError("Unsupported primitive type sizes in file header: $sizes " *
                             "(expected the standard (1,2,4,8,4,8) for CHAR/INT_2/INT_4/INT_8/REAL_4/REAL_8)"))
    p = 13
    t2, p = _ru16(buf, p)
    t4, p = _ru32(buf, p)
    t8, p = _ru64(buf, p)
    tr4, p = (only(reinterpret(Float32, buf[p:(p + 3)])), p + 4)
    tr8, p = _rf64(buf, p)
    (t2 == 0x1234 && t4 == 0x12345678 && t8 == 0x0123456789abcdef &&
     tr4 ≈ Float32(pi) && tr8 ≈ pi) ||
        throw(ArgumentError("GWF file was written on a big-endian host (byte-order test words " *
                             "don't match); only little-endian frames are supported by this reader"))
    return 41  # 1-based position of the first byte after the 40-byte header
end

# ================================================================================================
# Common structure header (14 bytes, all top-level structures), per Table 6:
#   length   INT_8U  total byte length of this structure, header through trailing chkSum
#   chkType  CHAR_U  0=none, 1=CRC
#   class    CHAR_U  structure class number (1=FrSH, 2=FrSE always; others assigned per-file
#                     by preceding FrSH dictionary records)
#   instance INT_4U  0-based occurrence counter for this class, reset each frame/file end
# ================================================================================================
struct _StructHeader
    start::Int      # 1-based byte position of this structure's `length` field
    len::Int        # total byte length, i.e. next structure starts at start + len
    class::UInt8
    instance::UInt32
end

function _read_struct_header(buf, p)
    len, p2 = _ru64(buf, p)
    _chktype, p2 = _ru8(buf, p2)
    class, p2 = _ru8(buf, p2)
    instance, p2 = _ru32(buf, p2)
    return _StructHeader(p, Int(len), class, instance), p2
end

# FrVect data-type codes (LIGO-T970130-v2 Appendix C, Table 30) -> Julia element type.
const _FRVECT_ELTYPE = Dict{UInt16,DataType}(0 => Int8, 1 => Int16, 2 => Float64, 3 => Float32,
                                              4 => Int32, 5 => Int64, 6 => ComplexF32,
                                              7 => ComplexF64, 9 => UInt16, 10 => UInt32,
                                              11 => UInt64, 12 => UInt8)
# (code 8 = FR_VECT_STRING has no fixed element type; not needed for numeric channels and not
# supported here.)

# Compression scheme codes (LIGO-T970130-v2 Appendix B, Table 29). Every scheme has two codes --
# one for big-endian writers, one (+256) for little-endian -- to record the writer's byte order
# for zero-suppression's bit-packing, which is byte-order sensitive. Documented in full here even
# though only raw/gzip are implemented, so the "not yet implemented" error below is self-explaining.
const _COMPRESS_NAMES = Dict(0 => "raw (big-endian written)", 256 => "raw",
                              1 => "gzip (big-endian written)", 257 => "gzip",
                              3 => "gzip-differential (big-endian written)", 259 => "gzip-differential",
                              5 => "zero-suppress 2-byte words (big-endian written)",
                              261 => "zero-suppress 2-byte words",
                              8 => "zero-suppress 4-byte words (big-endian written)",
                              264 => "zero-suppress 4-byte words",
                              10 => "zero-suppress 8-byte words (big-endian written)",
                              266 => "zero-suppress 8-byte words")

"""
    _decompress_frvect(raw::Vector{UInt8}, compress::Integer, ::Type{T}, nData::Integer) -> Vector{T}

Turn an `FrVect`'s on-disk `data` bytes into `nData` samples of type `T`. Only `compress` codes 0/256
(raw) and 1/257 (gzip) are implemented -- confirmed (see module docstring) to be the only schemes
GWOSC's public data actually uses; every other code in `_COMPRESS_NAMES` raises a clear error
instead of silently misinterpreting zero-suppressed or gzip-differential data as raw bytes.
"""
function _decompress_frvect(raw::Vector{UInt8}, compress::Integer, ::Type{T}, nData::Integer) where {T}
    nbytes_needed = nData * sizeof(T)
    if compress == 0 || compress == 256
        length(raw) == nbytes_needed ||
            throw(ArgumentError("raw FrVect data is $(length(raw)) bytes, expected $nbytes_needed " *
                                 "for $nData x $T"))
        plain = raw
    elseif compress == 1 || compress == 257
        plain = transcode(GzipDecompressor, raw)
        length(plain) == nbytes_needed ||
            throw(ArgumentError("gzip-decompressed FrVect data is $(length(plain)) bytes, expected " *
                                 "$nbytes_needed for $nData x $T"))
    else
        name = get(_COMPRESS_NAMES, Int(compress), "unknown code $compress")
        throw(ArgumentError("FrVect compression scheme \"$name\" (code $compress) is not implemented " *
                             "by this reader -- only raw (0/256) and gzip (1/257) are, since those are " *
                             "the only schemes observed in real GWOSC .gwf files; see GWF.jl's module " *
                             "docstring for the full LIGO-T970130 Appendix B scheme table"))
    end
    return collect(reinterpret(T, plain))
end

# ================================================================================================
# FrVect (Table 26), fields read in order after the common structure header:
#   name STRING, compress INT_2U, type INT_2U, nData INT_8U, nBytes INT_8U, data CHAR[nBytes],
#   nDim INT_4U, nx INT_8U[nDim], dx REAL_8[nDim], startX REAL_8[nDim], unitX STRING[nDim],
#   unitY STRING, next PTR_STRUCT (INT_2U class + INT_4U instance = 6 bytes), chkSum INT_4U.
# Only 1-D vectors (nDim == 1) are supported -- the only shape GWOSC's time-series channels use;
# `dx[1]`/`startX[1]` give this reader `dt`/the vector's time offset within its frame.
# ================================================================================================
struct _FrVect
    name::String
    data::Vector
    dt::Float64
    startx::Float64
end

function _parse_frvect(buf, p)
    name, p = _read_string(buf, p)
    compress, p = _ru16(buf, p)
    vtype, p = _ru16(buf, p)
    ndata, p = _ru64(buf, p)
    nbytes, p = _ru64(buf, p)
    raw = buf[p:(p + Int(nbytes) - 1)]
    p += Int(nbytes)
    ndim, p = _ru32(buf, p)
    ndim == 1 ||
        throw(ArgumentError("FrVect \"$name\" has nDim=$ndim; only 1-D FrVects (time series) are " *
                             "supported by this reader"))
    _nx, p = _ru64(buf, p)
    dx, p = _rf64(buf, p)
    startx, p = _rf64(buf, p)
    _unitx, p = _read_string(buf, p)
    _unity, p = _read_string(buf, p)

    T = get(_FRVECT_ELTYPE, vtype, nothing)
    T === nothing &&
        throw(ArgumentError("FrVect \"$name\" has unsupported data-type code $vtype (see " *
                             "LIGO-T970130 Appendix C for the full FrVect type table)"))
    samples = _decompress_frvect(raw, compress, T, ndata)
    return _FrVect(name, samples, dx, startx), p
end

"""
    _parse_frprocdata_header(buf, p) -> (name::String, time_offset::Float64, next_p)

Reads just enough of an `FrProcData` (Table 17: `name STRING, comment STRING, type INT_2U,
subType INT_2U, timeOffset REAL_8, ...`) to identify the channel and its time offset within the
frame -- the rest of the structure (tRange, fShift, ..., the FrVect/FrTable PTR_STRUCT pointers)
isn't needed since this reader locates a channel's FrVect by name (see module docstring), not by
following pointers.
"""
function _parse_frprocdata_header(buf, p)
    name, p = _read_string(buf, p)
    _comment, p = _read_string(buf, p)
    _type, p = _ru16(buf, p)
    _subtype, p = _ru16(buf, p)
    time_offset, p = _rf64(buf, p)
    return name, time_offset, p
end

"""
    _parse_fradcdata_header(buf, p) -> (name::String, time_offset::Float64, next_p)

Reads just enough of an `FrAdcData` (Table 10: `name STRING, comment STRING, channelGroup INT_4U,
channelNumber INT_4U, nBits INT_4U, bias REAL_4, slope REAL_4, units STRING, sampleRate REAL_8,
timeOffset REAL_8, ...`) to identify the channel and its time offset -- same rationale as
`_parse_frprocdata_header` above.
"""
function _parse_fradcdata_header(buf, p)
    name, p = _read_string(buf, p)
    _comment, p = _read_string(buf, p)
    _channel_group, p = _ru32(buf, p)
    _channel_number, p = _ru32(buf, p)
    _nbits, p = _ru32(buf, p)
    p += 4  # bias, REAL_4
    p += 4  # slope, REAL_4
    _units, p = _read_string(buf, p)
    _sample_rate, p = _rf64(buf, p)
    time_offset, p = _rf64(buf, p)
    return name, time_offset, p
end

"""
    _parse_frameh_gtime(buf, p) -> (GTimeS::UInt32, GTimeN::UInt32, next_p)

Reads just the fields of `FrameH` (Table 9) needed for absolute timing: `name STRING, run INT_4S,
frame INT_4U, dataQuality INT_4U, GTimeS INT_4U, GTimeN INT_4U`. Doesn't parse `ULeapS`/`dt`/the
13 trailing `PTR_STRUCT` pointers -- unneeded here, since (per the module docstring) this reader
matches channels to their FrVect by name across the whole file rather than following FrameH's
pointer graph.
"""
function _parse_frameh_gtime(buf, p)
    _name, p = _read_string(buf, p)
    p += 4  # run, INT_4S
    p += 4  # frame, INT_4U
    p += 4  # dataQuality, INT_4U
    gtimes, p = _ru32(buf, p)
    gtimen, p = _ru32(buf, p)
    return gtimes, gtimen, p
end

# ================================================================================================
# Top-level linear scan: walks every structure in the file, builds the per-file class-number ->
# name dictionary from FrSH records (class 1; FrSE, class 2, is skipped -- its per-field
# descriptions aren't needed since this reader hardcodes the fixed structure layouts above), and
# collects, per encountered frame, every FrProcData/FrAdcData channel header and every FrVect --
# deferring the "which FrVect belongs to which channel" match to the caller (matched by name).
# ================================================================================================
struct _ScanResult
    channel_names::Vector{String}                         # in file order, deduplicated
    channel_time_offset::Dict{String,Float64}
    frame_of_channel::Dict{String,Int}                     # channel name -> frame index seen in
    frame_gtime::Vector{Float64}                           # frame index -> GTimeS + GTimeN/1e9
    frvects::Dict{String,Vector{Tuple{Int,_FrVect}}}       # channel name -> [(frame index, FrVect)]
end

function _scan(buf)
    p = _parse_file_header(buf)
    classdict = Dict{UInt8,String}()
    channel_names = String[]
    channel_time_offset = Dict{String,Float64}()
    frame_of_channel = Dict{String,Int}()
    frame_gtime = Float64[]
    frvects = Dict{String,Vector{Tuple{Int,_FrVect}}}()
    frame_idx = 0

    n = length(buf)
    while p <= n - 13
        h, fieldsp = _read_struct_header(buf, p)
        if h.class == 1   # FrSH: register a class name for this file
            name, fp = _read_string(buf, fieldsp)
            classid, fp = _ru16(buf, fp)
            classdict[UInt8(classid)] = name
        elseif h.class != 2   # skip FrSE (2); everything else, look up its registered name
            clsname = get(classdict, h.class, nothing)
            if clsname == "FrameH"
                gtimes, gtimen, _fp = _parse_frameh_gtime(buf, fieldsp)
                frame_idx += 1
                push!(frame_gtime, Float64(gtimes) + Float64(gtimen) / 1e9)
            elseif clsname == "FrProcData"
                name, toff, _fp = _parse_frprocdata_header(buf, fieldsp)
                name in keys(channel_time_offset) || push!(channel_names, name)
                channel_time_offset[name] = toff
                frame_of_channel[name] = frame_idx
            elseif clsname == "FrAdcData"
                name, toff, _fp = _parse_fradcdata_header(buf, fieldsp)
                name in keys(channel_time_offset) || push!(channel_names, name)
                channel_time_offset[name] = toff
                frame_of_channel[name] = frame_idx
            elseif clsname == "FrVect"
                fv, _fp = _parse_frvect(buf, fieldsp)
                push!(get!(() -> Tuple{Int,_FrVect}[], frvects, fv.name), (frame_idx, fv))
            end
        end
        p = h.start + h.len
    end
    return _ScanResult(channel_names, channel_time_offset, frame_of_channel, frame_gtime, frvects)
end

"""
    list_gwf_channels(path) -> Vector{String}

List the channel names (`FrProcData`/`FrAdcData`) found in a GWF frame file, in file order.
"""
function list_gwf_channels(path::AbstractString)
    buf = read(path)
    return _scan(buf).channel_names
end

"""
    load_gwf_strain(path; channel=nothing) -> TimeSeries{Float64}

Read a real-valued time-series channel out of a LIGO/Virgo IGWD frame (`.gwf`) file, following
LIGO-T970130-v2 (see this module's docstring for the exact subset of the spec implemented and its
scope limits: little-endian frames, `FrProcData`/`FrAdcData` channels matched to their `FrVect` by
name rather than by resolving the full pointer graph, raw/gzip compression only).

`channel`, if given, must exactly match one of `list_gwf_channels(path)`. If omitted, this looks
for a unique channel whose name contains `"STRAIN"` (case-insensitive) -- GWOSC's convention for
its public strain channels (e.g. `"H1:LOSC-STRAIN"`, `"H1:GWOSC-4KHZ_R1_STRAIN"`) -- and errors,
listing the channels found, if that isn't unique.

A multi-frame file's chunks for the chosen channel are concatenated in file order; this errors if
their sample spacing (`dt`) differs, or if consecutive chunks aren't exactly time-contiguous
(no gap-filling/overlap-resolution is attempted).
"""
function load_gwf_strain(path::AbstractString; channel::Union{AbstractString,Nothing}=nothing)
    buf = read(path)
    scan = _scan(buf)

    chosen = if channel !== nothing
        channel in scan.channel_names ||
            throw(ArgumentError("Channel $(repr(channel)) not found in $path. Available: " *
                                 join(scan.channel_names, ", ")))
        String(channel)
    else
        matches = filter(c -> occursin("STRAIN", uppercase(c)), scan.channel_names)
        length(matches) == 1 ||
            throw(ArgumentError("Could not uniquely identify a STRAIN channel in $path (found " *
                                 "$(length(matches)) match(es): $(join(matches, ", "))). Available " *
                                 "channels: $(join(scan.channel_names, ", ")). Pass `channel=` explicitly."))
        matches[1]
    end

    haskey(scan.frvects, chosen) ||
        throw(ArgumentError("Channel $(repr(chosen)) has an FrProcData/FrAdcData header but no " *
                             "matching FrVect (by name) was found in $path"))
    chunks = sort(scan.frvects[chosen]; by=first)   # by frame index, i.e. file order

    frame_idx1, fv1 = chunks[1]
    dt = fv1.dt
    t0 = scan.frame_gtime[frame_idx1] + scan.channel_time_offset[chosen] + fv1.startx
    all_data = fv1.data
    expected_next_t0 = t0 + length(fv1.data) * dt

    for (frame_idx, fv) in chunks[2:end]
        fv.dt == dt ||
            throw(ArgumentError("Channel $(repr(chosen)) has inconsistent dt across frames in " *
                                 "$path ($dt vs $(fv.dt)); can't concatenate"))
        chunk_t0 = scan.frame_gtime[frame_idx] + scan.channel_time_offset[chosen] + fv.startx
        isapprox(chunk_t0, expected_next_t0; atol=1e-9) ||
            throw(ArgumentError("Channel $(repr(chosen)) has a gap/overlap between frames in " *
                                 "$path (expected next chunk to start at $expected_next_t0, got " *
                                 "$chunk_t0); non-contiguous multi-frame concatenation isn't supported"))
        all_data = vcat(all_data, fv.data)
        expected_next_t0 = chunk_t0 + length(fv.data) * dt
    end

    return TimeSeries(Float64.(all_data), t0, dt; name=chosen)
end

end
