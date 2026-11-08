# Data

Local copies of public GWOSC data files used by the test suite. Everything in this directory
except this file is gitignored. The tests also look in `ENV["GRAVITATIONALWAVEDATA_TEST_DATA"]` and in the
`fetch_open_data` cache (`~/.cache/gravitationalwavedata/strain`); see `test/testdata.jl`. Tests whose files
are missing are skipped.

Files used (all GW150914, GPS `1126256640`-`1126260736`, 4096 s, 4 kHz):

- `H-H1_LOSC_4_V1-1126256640-4096.hdf5`, `L-L1_LOSC_4_V1-1126256640-4096.hdf5`: GWOSC HDF5
  strain files. Obtain with `fetch_open_data("H1", 1126256640, 1126260736)` (and likewise
  `"L1"`), which saves them to the cache.
- `H-H1_LOSC_4_V1-1126256640-4096.gwf`: the same H1 segment in native IGWD frame format, read
  by `load_gwf_strain`/`list_gwf_channels`. Contains channels `H1:LOSC-STRAIN` (REAL_8,
  gzip-compressed), `H1:LOSC-DQMASK`, and `H1:LOSC-INJMASK`. Download it into this directory
  from `https://gwosc.org/archive/data/O1/1126170624/H-H1_LOSC_4_V1-1126256640-4096.gwf`
  (listed by `https://gwosc.org/archive/links/O1/H1/1126256640/1126260736/json/`).
