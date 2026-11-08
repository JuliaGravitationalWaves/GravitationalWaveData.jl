#!/usr/bin/env julia
# Reference side of compare_time_delay.jl (see README.md): writes
# lal.TimeDelayFromEarthCenter over the same grid to a CSV. Runs in this
# directory's project and must not load GWData.
#
# Usage: julia --project=examples/validation compare_time_delay_lal_step.jl <delay_csv>

using PythonCall
const lal = pyimport("lal")

delay_csv = ARGS[1]

# lal.cached_detector_by_prefix maps "H1"/"L1"/"V1"/... to the corresponding LALDetector, the
# same objects the lalCachedDetectors C array holds (LALDetectors.h / CreateDetector.c).
const PREFIXES = ["H1", "L1", "V1"]

const POINTS = [
    (0.0, 0.0, 1126259462.4),
    (1.0, 0.5, 1126259462.4),
    (3.0, -1.0, 1187008882.4),
    (0.0, pi / 2 - 0.001, 1126259462.4),
    (2.5, -(pi / 2 - 0.001), 1126259462.4),
    (4.2, 0.05, 1119744016.0),
    (5.9, 1.4, 1440000000.0),
    (0.7, -0.3, 630763213.0),
]

open(delay_csv, "w") do io
    println(io, "det,ra,dec,gps,delay")
    for name in PREFIXES
        det = lal.cached_detector_by_prefix[name]
        loc = det.location
        for (ra, dec, gps) in POINTS
            gpstime = lal.LIGOTimeGPS(gps)
            delay = lal.TimeDelayFromEarthCenter(loc, ra, dec, gpstime)
            d = pyconvert(Float64, delay)
            println(io, "$name,$(repr(ra)),$(repr(dec)),$(repr(gps)),$(repr(d))")
        end
    end
end
