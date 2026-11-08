#!/usr/bin/env julia
# GWData.jl side of compare_time_delay.jl (see README.md): writes
# `time_delay_from_geocenter` over a grid of (det, ra, dec, gps) points to a CSV.
#
# Usage: julia --project=<GWData.jl root> compare_time_delay_gwdata_step.jl <delay_csv>

using GWData: H1, L1, V1, time_delay_from_geocenter

delay_csv = ARGS[1]

const DETECTORS = [("H1", H1), ("L1", L1), ("V1", V1)]

# (ra, dec, gps) test points: a generic spread plus edge cases (near the celestial poles,
# a leap-second-boundary GPS time, the J2000.0 epoch) -- same convention/point set as
# compare_antenna_response's POINTS (minus psi, which this function doesn't take).
const POINTS = [
    (0.0, 0.0, 1126259462.4),
    (1.0, 0.5, 1126259462.4),
    (3.0, -1.0, 1187008882.4),
    (0.0, pi / 2 - 0.001, 1126259462.4),    # near north celestial pole
    (2.5, -(pi / 2 - 0.001), 1126259462.4), # near south celestial pole
    (4.2, 0.05, 1119744016.0),              # leap-second-boundary GPS time
    (5.9, 1.4, 1440000000.0),
    (0.7, -0.3, 630763213.0),               # J2000.0 epoch
]

open(delay_csv, "w") do io
    println(io, "det,ra,dec,gps,delay")
    for (name, det) in DETECTORS
        for (ra, dec, gps) in POINTS
            delay = time_delay_from_geocenter(det, ra, dec, gps)
            println(io, "$name,$(repr(ra)),$(repr(dec)),$(repr(gps)),$(repr(delay))")
        end
    end
end
