# Notices

GravitationalWaveData.jl  
Copyright (C) 2026 Paul Hopkins

This program is free software: you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. See [LICENSE](LICENSE).

## Upstream code

Parts of this package are Julia ports of code from the projects below, which are distributed under
the GNU General Public License. The source comments name the routines concerned.

| Upstream | Licence | Used in |
|---|---|---|
| [GWpy](https://github.com/gwpy/gwpy), Copyright (c) 2014-2017 Louisiana State University, 2017-2025 Cardiff University | GPL-3.0-or-later | `Filtering.jl`, `QTransform.jl` |
| [LALSuite](https://git.ligo.org/lscsoft/lalsuite), Copyright (c) the LALSuite authors | GPL-2.0-or-later | `Detectors.jl`, `SiderealTime.jl`, `ReferenceCurves.jl` |

## Earlier snapshot

The tag `paper1-run` records the source tree used for the runs of the accompanying paper, unchanged
so that its hash stays verifiable. Its files say MIT, but it contains the same upstream-derived
code, so that snapshot is also distributed under GPL-3.0-or-later.
