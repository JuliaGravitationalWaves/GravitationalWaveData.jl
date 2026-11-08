# Locating the real GWOSC data files some tests use. Searched in order: this repo's `data/`
# directory, the directory named by `ENV["GWDATA_TEST_DATA"]` (if set), and the default
# `fetch_open_data` strain cache (`gwdata_cache_dir("strain")`, by default
# `~/.cache/gwdata/strain`). Tests needing a file that isn't found `@test_skip` rather than fail;
# run e.g. `fetch_open_data("H1", 1126256640, 1126260736)` once to populate the cache.

const TEST_DATA_DIRS = filter(!isnothing, [
    joinpath(@__DIR__, "..", "data"),
    get(ENV, "GWDATA_TEST_DATA", nothing),
    gwdata_cache_dir("strain"),
])

"""
    find_test_data(filename) -> Union{String,Nothing}

Path to `filename` in the first of `TEST_DATA_DIRS` that has it, or `nothing`.
"""
function find_test_data(filename::AbstractString)
    for dir in TEST_DATA_DIRS
        path = joinpath(dir, filename)
        isfile(path) && return path
    end
    return nothing
end

const GW150914_H1_HDF5 = "H-H1_LOSC_4_V1-1126256640-4096.hdf5"
const GW150914_L1_HDF5 = "L-L1_LOSC_4_V1-1126256640-4096.hdf5"
const GW150914_H1_GWF = "H-H1_LOSC_4_V1-1126256640-4096.gwf"
