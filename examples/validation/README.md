# Validation against gwpy and LALSuite

These scripts check GWData.jl numerically against the reference implementations it follows:
[gwpy](https://gwpy.github.io/) (and the scipy routines underneath it) for spectral estimation,
filtering, resampling, spectrograms, coherence and the Q-transform, and LALSuite's `lal` and
`lalsimulation` Python bindings for sidereal time, detector response, time delays and design
noise curves. The agreement they measure is listed in the
[Validation section of the conventions page](../../docs/CONVENTIONS.md#Validation), and the spot
values are kept as regression tests under `test/`, so you only need these scripts to re-check a
change against the references.

## Requirements

- Linux or macOS. LALSuite has no Windows build, so on Windows run these inside WSL2 with a
  Linux install of Julia.
- Julia 1.10 or later.
- Network access on the first run. The data-driven comparisons download the GW150914 H1 and
  L1 strain files from GWOSC (about 130 MB each) into the GWData.jl cache, and the reference
  project creates a conda environment with gwpy and LALSuite (a few hundred MB). Both are
  reused on later runs.

Nothing else needs installing by hand. [CondaPkg.jl](https://github.com/JuliaPy/CondaPkg.jl)
creates the Python environment from [`CondaPkg.toml`](CondaPkg.toml) the first time a reference
step runs. If you have set `JULIA_CONDAPKG_BACKEND=Null` and `JULIA_PYTHONCALL_EXE` to use your
own Python, that Python needs `gwpy`, `lalsuite` and `numpy` installed.

## Running

From the repository root:

```sh
julia examples/validation/run_all.jl              # every comparison
julia examples/validation/run_all.jl psd gmst     # a subset
julia examples/validation/compare_psd.jl          # one comparison
```

Each comparison prints its agreement figures. The intermediate CSV files go to
`examples/validation/output/`, which git ignores.

| Script | GWData.jl | Reference |
|:--|:--|:--|
| `compare_psd.jl` | `psd` | gwpy `TimeSeries.psd` |
| `compare_filtering.jl` | `bandpass`, `highpass`, `lowpass`, `whiten` | gwpy filters and `whiten` |
| `compare_resample.jl` | `resample` | gwpy `TimeSeries.resample` |
| `compare_spectrogram.jl` | `spectrogram` | gwpy `TimeSeries.spectrogram2` |
| `compare_qtransform.jl` | `q_transform` | gwpy `q_scan` |
| `compare_coherence.jl` | `csd`, `coherence` | gwpy `TimeSeries.csd`, `coherence` |
| `compare_gmst.jl` | `gmst` | `lal.GreenwichMeanSiderealTime` |
| `compare_antenna_response.jl` | detector geometry, `compute_antenna_response` | `lal.CachedDetectors`, `lal.ComputeDetAMResponse` |
| `compare_time_delay.jl` | `time_delay_from_geocenter` | `lal.TimeDelayFromEarthCenter` |
| `compare_reference_curves.jl` | `aligo_design_psd`, `avirgo_design_psd` | `lalsimulation.SimNoisePSD*` |

`qtile_reference_step.jl` is not part of `run_all.jl`. It regenerates the single-tile reference
values used by `test/qtransform.jl`; see its header for usage.

## How the scripts are organised

Each `compare_<name>.jl` driver runs two separate Julia processes, which pass their results to
each other through CSV files:

- `compare_<name>_gwdata_step.jl` runs in the GWData.jl project (the repository root).
- `compare_<name>_gwpy_step.jl` or `compare_<name>_lal_step.jl` runs in this directory's project,
  which contains only PythonCall plus the conda environment.

The two sides must not share a process. Julia's `HDF5_jll` (a GWData.jl dependency) and
Python's `h5py` each ship a library named `libhdf5`, and with both loaded the wrong one's
symbols get used, so `import h5py` fails with an `undefined symbol` error. The shared helpers for
starting the two processes and locating the data are in [`common.jl`](common.jl).
