module Resample

export resample

using ..Types: TimeSeries, sample_rate
import DSP

# ================================================================================================
# Resampling: `resample(ts::TimeSeries, new_rate)` wraps `DSP.resample`, which is itself a
# polyphase FIR resampler (rational-ratio path) / polyphase-with-linear-interpolation resampler
# (arbitrary-ratio path). Both paths in DSP.jl design their own anti-aliasing/anti-imaging FIR
# filter via `DSP.resample_filter`: a Kaiser-windowed lowpass with cutoff at the *new* Nyquist
# frequency (`rel_bw=1.0` by default) and `attenuation=60` dB stopband attenuation -- i.e. GW-
# appropriate anti-aliasing-before-downsampling is DSP.jl's *default* behavior, not something that
# needs to be composed on top: see `DSP.Filters.resample_filter` docstrings and
# `_resample_filter` in `DSP/src/Filters/design.jl` (`cutoff = f_nyq * rel_bw`,
# `trans_width = cutoff * 0.2`, `kaiserord(trans_width, attenuation)`). `DSP.resample` also
# compensates for the FIR filter's own group delay (`undelay!`/`setphase!` in
# `DSP/src/Filters/stream_filt.jl`) so the output lines up in time with the input rather than
# being shifted by the filter's ramp-up -- this is why `resample` below keeps `t0` unchanged
# rather than trying to correct for a filter delay itself.
#
# Rational vs. arbitrary ratio: GW sample rates are essentially always integer Hz related by a
# power-of-two-ish ratio (4096 -> 2048, 16384 -> 4096, 4096 -> 1024, ...). When both the current
# and requested rate are exact integers we build an exact `Rational` ratio and use DSP's rational
# (pure polyphase, no interpolation error) path; otherwise we fall back to DSP's `AbstractFloat`
# arbitrary-rate path (polyphase + linear interpolation between phases), which is approximate but
# handles any ratio.
# ================================================================================================

"""
    resample(ts::TimeSeries, new_rate::Real; attenuation::Real=60, rel_bw::Real=1.0,
              Nphi::Integer=32) -> TimeSeries

Resample `ts` to `new_rate` Hz using `DSP.resample`'s polyphase FIR resampler, with an
automatically-designed Kaiser-window anti-aliasing (when downsampling) / anti-imaging (when
upsampling) lowpass filter at the new Nyquist frequency -- matching `DSP.resample`'s own defaults,
which are already GW-appropriate (see module docs above).

  - `attenuation`: stopband attenuation of the anti-aliasing filter, in dB (`DSP.resample_filter`'s
    default is `60`).
  - `rel_bw`: filter cutoff as a fraction of the new Nyquist frequency (`DSP.resample_filter`'s
    default is `1.0`, i.e. cutoff exactly at the new Nyquist).
  - `Nphi`: number of polyphase filter phases, only used for the arbitrary-ratio (non-integer rate
    ratio) path (`DSP.resample_filter`'s default is `32`).

`t0` is preserved (see module docs: `DSP.resample` already compensates for the anti-aliasing
filter's own group delay), and `dt` is set to exactly `1 / new_rate`.

Returns `ts` unchanged (a copy) if `new_rate == sample_rate(ts)`.
"""
function resample(ts::TimeSeries, new_rate::Real; attenuation::Real=60, rel_bw::Real=1.0,
                   Nphi::Integer=32)
    new_rate > 0 || throw(ArgumentError("new_rate must be positive, got $new_rate"))
    old_rate = sample_rate(ts)

    if new_rate == old_rate
        return TimeSeries(copy(ts.data), ts.t0, ts.dt; name=ts.name)
    end

    x = float.(ts.data)

    if isinteger(new_rate) && isinteger(old_rate)
        ratio = Rational(Int(new_rate), Int(old_rate))
        h = DSP.resample_filter(ratio, rel_bw, attenuation)
        resampled = DSP.resample(x, ratio, h)
    else
        ratio = Float64(new_rate) / Float64(old_rate)
        h = DSP.resample_filter(ratio, Nphi, rel_bw, attenuation)
        resampled = DSP.resample(x, ratio, h, Nphi)
    end

    return TimeSeries(resampled, ts.t0, 1 / Float64(new_rate); name=ts.name)
end

end
