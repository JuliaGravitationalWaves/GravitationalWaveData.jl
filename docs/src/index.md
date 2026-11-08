# GravitationalWaveData.jl

A Julia toolkit for gravitational-wave detector data: `TimeSeries`/`FrequencySeries`/`Spectrogram`
types, GWOSC and frame-file I/O, PSD/ASD estimation, filtering and whitening, resampling,
spectrograms, Q-transforms, GPS time and segments, detector geometry and antenna response, and
coherence, built on `DSP.jl` and `FFTW.jl`.

## Where to start

- **New to the package?** Start with [Getting Started](generated/01_getting_started.md): loading
  real GWOSC strain data into a `TimeSeries` and making a first plot.
- **Looking for a specific topic?** Each tutorial in the sidebar covers one module, and the
  [GW150914 capstone](generated/12_gw150914_capstone.md) combines them to reproduce the
  discovery figure.
- **Want the full function/type list?** See the [API Reference](api/index.md).
- **Want the numerical conventions** (PSD normalization, whitening algorithm, GPS/leap seconds,
  ...) and how they were validated? See [Conventions](conventions.md).
- **Prefer a notebook?** Every tutorial is also available as a Jupyter notebook under
  `notebooks/`, and the Q-transform and filter-design topics have interactive explorers, as
  `Pluto.jl` notebooks under `notebooks/pluto/` and standalone `Bonito.jl` web apps under
  `notebooks/bonito/`.

## Installation

```julia
using Pkg
Pkg.add("GravitationalWaveData")
```
