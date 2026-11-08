#!/usr/bin/env julia
# Reference side of compare_antenna_response.jl (see README.md): writes
# lal.CachedDetectors' H1/L1/V1 geometry and lal.ComputeDetAMResponse over the same grid. Runs in this
# directory's project and must not load GWData.
#
# Usage: julia --project=examples/validation compare_antenna_response_lal_step.jl <detectors_csv> <antenna_csv>

using PythonCall
const lal = pyimport("lal")
np = pyimport("numpy")

detectors_csv, antenna_csv = ARGS

# lal.cached_detector_by_prefix maps "H1"/"L1"/"V1"/... to the corresponding LALDetector, the
# same objects the lalCachedDetectors C array holds (LALDetectors.h / CreateDetector.c).
const PREFIXES = ["H1", "L1", "V1"]

open(detectors_csv, "w") do io
    println(io, "det,loc_x,loc_y,loc_z,r11,r12,r13,r21,r22,r23,r31,r32,r33")
    for name in PREFIXES
        det = lal.cached_detector_by_prefix[name]
        loc = pyconvert(Vector{Float64}, np.asarray(det.location))
        resp = pyconvert(Matrix{Float64}, np.asarray(det.response))
        println(io, "$name,$(repr(loc[1])),$(repr(loc[2])),$(repr(loc[3])),",
            "$(repr(resp[1, 1])),$(repr(resp[1, 2])),$(repr(resp[1, 3])),",
            "$(repr(resp[2, 1])),$(repr(resp[2, 2])),$(repr(resp[2, 3])),",
            "$(repr(resp[3, 1])),$(repr(resp[3, 2])),$(repr(resp[3, 3]))")
    end
end

const POINTS = [
    (0.0, 0.0, 0.0, 1126259462.4),
    (1.0, 0.5, 0.3, 1126259462.4),
    (3.0, -1.0, 1.2, 1187008882.4),
    (0.0, pi / 2 - 0.001, 0.0, 1126259462.4),
    (2.5, -(pi / 2 - 0.001), 0.7, 1126259462.4),
    (4.2, 0.05, pi / 2, 1119744016.0),
    (5.9, 1.4, 2.9, 1440000000.0),
    (0.7, -0.3, 0.0, 630763213.0),
]

open(antenna_csv, "w") do io
    println(io, "det,ra,dec,psi,gps,fplus,fcross")
    for name in PREFIXES
        det = lal.cached_detector_by_prefix[name]
        resp = det.response
        for (ra, dec, psi, gps) in POINTS
            gmst = lal.GreenwichMeanSiderealTime(lal.LIGOTimeGPS(gps))
            result = lal.ComputeDetAMResponse(resp, ra, dec, psi, gmst)
            fp = pyconvert(Float64, result[0])
            fc = pyconvert(Float64, result[1])
            println(io, "$name,$(repr(ra)),$(repr(dec)),$(repr(psi)),$(repr(gps)),$(repr(fp)),$(repr(fc))")
        end
    end
end
