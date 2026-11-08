module GWData

include("GWData/Types.jl")
include("GWData/IO/GWOSC.jl")
include("GWData/IO/Posteriors.jl")
include("GWData/GPSTime.jl")
include("GWData/Segments.jl")
include("GWData/IO/GWF.jl")
include("GWData/PSD.jl")
include("GWData/Coherence.jl")
include("GWData/Resample.jl")
include("GWData/Filtering.jl")
include("GWData/SiderealTime.jl")
include("GWData/Detectors.jl")
include("GWData/QTransform.jl")
include("GWData/Spectrogram.jl")
include("GWData/ReferenceCurves.jl")

using .Types
using .GWOSC
using .Posteriors
using .GPSTime
using .Segments
using .GWF
using .PSD
using .Coherence
using .Resample
using .Filtering
using .SiderealTime
using .Detectors
using .QTransform
using .Spectrograms
using .ReferenceCurves

export TimeSeries, FrequencySeries, times, frequencies, duration, sample_rate, crop
export load_gwosc_strain, condition_segment, fetch_open_data, gwdata_cache_dir
export fetch_posterior_samples, list_posterior_analyses, load_posterior_samples
export GPS_EPOCH, gps_to_utc, utc_to_gps, leap_seconds
export Segment, SegmentList, coalesce_segments, gaps, overlaps
export load_gwf_strain, list_gwf_channels
export psd, asd
export csd, coherence
export resample
export bandpass, highpass, lowpass, notch, whiten
export gmst
export AbstractDetector, Detector, DetectorGeometry, geometry, H1, L1, V1
export compute_antenna_response, antenna_response_gmst, time_delay_from_geocenter, time_delay_from_geocenter_gmst
export normalized_fft, QTileParams, qtile_window, qtile_data_indices, qtile_padding,
       qtile_energy, q_values, tiling_frange, plane_frange, plane_frequencies,
       QRow, QScan, QPeak, energy_series, q_plane_transform, q_scan, q_transform, qscan_peak
export Spectrogram, spectrogram
export aligo_design_psd, aligo_design_asd, avirgo_design_psd, avirgo_design_asd

end
