#!/usr/bin/env julia
# GravitationalWaveData.jl side of compare_antenna_response.jl (see README.md): writes the H1/L1/V1
# geometry (location + response tensor) and `compute_antenna_response` over a grid of
# (ra, dec, psi, gps) points to CSVs.
#
# Usage: julia --project=<GravitationalWaveData.jl root> compare_antenna_response_gwdata_step.jl <detectors_csv> <antenna_csv>

using GravitationalWaveData: H1, L1, V1, compute_antenna_response, geometry

detectors_csv, antenna_csv = ARGS

const DETECTORS = [("H1", H1), ("L1", L1), ("V1", V1)]

open(detectors_csv, "w") do io
    println(io, "det,loc_x,loc_y,loc_z,r11,r12,r13,r21,r22,r23,r31,r32,r33")
    for (name, det) in DETECTORS
        lx, ly, lz = geometry(det).location
        r = geometry(det).response
        println(io, "$name,$(repr(lx)),$(repr(ly)),$(repr(lz)),",
            "$(repr(r[1][1])),$(repr(r[1][2])),$(repr(r[1][3])),",
            "$(repr(r[2][1])),$(repr(r[2][2])),$(repr(r[2][3])),",
            "$(repr(r[3][1])),$(repr(r[3][2])),$(repr(r[3][3]))")
    end
end

# (ra, dec, psi, gps) test points: a generic spread plus edge cases (near the celestial poles,
# i.e. near a detector's local zenith/nadir depending on latitude; a leap-second-boundary GPS
# time; psi at 0/pi/2 boundary values).
const POINTS = [
    (0.0, 0.0, 0.0, 1126259462.4),
    (1.0, 0.5, 0.3, 1126259462.4),
    (3.0, -1.0, 1.2, 1187008882.4),
    (0.0, pi / 2 - 0.001, 0.0, 1126259462.4),   # near north celestial pole
    (2.5, -(pi / 2 - 0.001), 0.7, 1126259462.4), # near south celestial pole
    (4.2, 0.05, pi / 2, 1119744016.0),          # leap-second-boundary GPS time, psi=pi/2
    (5.9, 1.4, 2.9, 1440000000.0),
    (0.7, -0.3, 0.0, 630763213.0),               # J2000.0 epoch
]

open(antenna_csv, "w") do io
    println(io, "det,ra,dec,psi,gps,fplus,fcross")
    for (name, det) in DETECTORS
        for (ra, dec, psi, gps) in POINTS
            fp, fc = compute_antenna_response(det, ra, dec, psi, gps)
            println(io, "$name,$(repr(ra)),$(repr(dec)),$(repr(psi)),$(repr(gps)),$(repr(fp)),$(repr(fc))")
        end
    end
end
