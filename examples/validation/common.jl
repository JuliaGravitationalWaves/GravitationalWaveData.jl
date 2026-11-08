# Shared helpers for the comparison drivers in this directory (see README.md).
#
# Every comparison runs as two separate Julia processes that exchange results through CSV
# files: one in the GravitationalWaveData.jl project, one in this directory's reference project (PythonCall
# with gwpy and LALSuite's Python bindings from conda-forge). They must stay separate because
# Julia's HDF5_jll (a GravitationalWaveData.jl dependency) and Python's h5py each ship a `libhdf5`, and loading
# both into one process lets the wrong library's symbols win (an `undefined symbol` crash on
# `import h5py`).

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const REFERENCE_PROJECT = @__DIR__
const OUT_DIR = joinpath(@__DIR__, "output")
mkpath(OUT_DIR)

# GPS span of the GWOSC GW150914 4 kHz strain files (`H-H1_LOSC_4_V1-1126256640-4096.hdf5` and
# its L1 counterpart) that the data-driven comparisons read.
const GW150914_GPS_START = 1126256640
const GW150914_GPS_END = 1126256640 + 4096

_julia() = Base.julia_cmd()

function _ensure_instantiated(project)
    isfile(joinpath(project, "Manifest.toml")) && return
    run(`$(_julia()) --project=$project -e "using Pkg; Pkg.instantiate()"`)
end

"""
    run_gwdata_step(script, args...)

Run `script` (a file in this directory) in the GravitationalWaveData.jl project, passing `args` as `ARGS`.
"""
function run_gwdata_step(script, args...)
    _ensure_instantiated(REPO_ROOT)
    run(`$(_julia()) --project=$REPO_ROOT $(joinpath(@__DIR__, script)) $(string.(args))`)
end

"""
    run_reference_step(script, args...)

Run `script` (a file in this directory) in the reference project, passing `args` as `ARGS`.
The first run instantiates the project and lets CondaPkg create its conda environment, which
downloads gwpy and LALSuite's Python bindings (a few hundred MB).
"""
function run_reference_step(script, args...)
    _ensure_instantiated(REFERENCE_PROJECT)
    run(`$(_julia()) --project=$REFERENCE_PROJECT $(joinpath(@__DIR__, script)) $(string.(args))`)
end

"""
    gw150914_strain_path(ifo) -> String

Local path of the GWOSC GW150914 strain file for detector `ifo` ("H1" or "L1"), downloaded into
the GravitationalWaveData.jl cache on first use. Runs `fetch_open_data` in a GravitationalWaveData.jl subprocess, since the
driver process itself loads neither project.
"""
function gw150914_strain_path(ifo)
    _ensure_instantiated(REPO_ROOT)
    code = "using GravitationalWaveData; print(fetch_open_data(ARGS[1], $GW150914_GPS_START, $GW150914_GPS_END))"
    return readchomp(`$(_julia()) --project=$REPO_ROOT -e $code $ifo`)
end

output_path(name) = joinpath(OUT_DIR, name)
