module ReferenceCurves

export aligo_design_psd, aligo_design_asd, avirgo_design_psd, avirgo_design_asd

using ..Types: FrequencySeries
using ..Detectors: _C_SI

# ================================================================================================
# Reference design-sensitivity ASD/PSD curves for aLIGO ("zero-detuning, high-power" configuration)
# and AdvVirgo, ported from LALSimulation's `LALSimNoisePSD.c`.
#
# Source (fetched 2026-08-28 from the canonical LALSuite repository, `master` branch, at
# https://git.ligo.org/lscsoft/lalsuite):
#   - lalsimulation/lib/LALSimNoisePSD.c -- the two top-level functions ported here are
#     `XLALSimNoisePSDaLIGOZeroDetHighPower` and `XLALSimNoisePSDAdvVirgo`.
#   - lal/lib/std/LALConstants.h -- the handful of SI physical constants (`LAL_PI`, `LAL_PI_180`,
#     `LAL_C_SI`, `LAL_K_SI`, `LAL_HBAR_SI`) these formulas need.
#
# IMPORTANT naming note: the well-known
# simple closed-form "Ajith fit" to the aLIGO zero-detuning-high-power curve (a short rational/
# exponential function of `f`, as used by e.g. PyCBC's and bilby's own `aLIGOZeroDetHighPower`
# convenience functions) is NOT what LALSimulation's `XLALSimNoisePSDaLIGOZeroDetHighPower` (and
# hence `lalsimulation.SimNoisePSDaLIGOZeroDetHighPower`, what this module is cross-validated
# against) actually computes. The real LALSimulation function instead composes:
#   - `XLALSimNoisePSDaLIGOQuantumZeroDetHighPower` (quantum shot + radiation-pressure noise),
#     itself `XLALSimNoisePSDQuantum` evaluated at aLIGO's design parameters -- a full
#     Buonanno & Chen, Phys. Rev. D 64, 062006 (2001) input-output relation for a signal-recycled
#     interferometer (their Eq. 5.13, ported from LALSimulation's own adaptation of GWINC's
#     `shotrad.m`), plus
#   - `XLALSimNoisePSDaLIGOThermal` (suspension + mirror-coating thermal noise), a sum of two
#     phenomenological power-law fits (`XLALSimNoisePSDSuspTherm`/`XLALSimNoisePSDMirrorTherm`).
# This module ports that actual composition verbatim (not the simpler Ajith-fit approximation),
# since that is what `lalsimulation.SimNoisePSDaLIGOZeroDetHighPower` returns and what the
# cross-check in `examples/validation/compare_reference_curves.jl` validates against.
#
# AdvVirgo, by contrast, genuinely is a simple closed-form phenomenological log-polynomial/
# exponential fit -- `XLALSimNoisePSDAdvVirgo`, Eq. 6 of Manzotti & Dietz 2012 (arXiv:1202.4031,
# cited in LALSimulation's own bibliography as "md2012"), fit to the AdvVirgo design curve
# published at http://wwwcascina.virgo.infin.it/advirgo.
# ================================================================================================

# --- SI physical constants, ported verbatim from lal/lib/std/LALConstants.h. ---
const _LAL_PI = 3.141592653589793238462643383279502884
const _LAL_PI_180 = 1.745329251994329576923690768488612713e-2
const _LAL_C_SI = _C_SI                                    # speed of light, m/s (Detectors.jl)
const _LAL_K_SI = 1.380649e-23                              # Boltzmann constant, J/K
const _LAL_HBAR_SI = 1.054571817646156391262428003302280745e-34  # reduced Planck constant, J*s

# --- aLIGO design parameters, ported verbatim from LALSimNoisePSD.c's `#define` constants
# --- ("these ones seem vaguely well motivated" / "these ones are not physically motivated but
# --- seem to be phenomenologically rightish", per that file's own comments). ---
const _ALIGO_ARMLENGTH_SI = 3995.0
const _ALIGO_LASER_POWER_HIGH_SI = 125.0
const _ALIGO_LASER_WAVELENGTH_SI = 1.064e-6
const _ALIGO_MIRROR_MASS_SI = 40.0
const _ALIGO_MIRROR_LOSS = 37.5e-6
const _ALIGO_BS_LOSS = 0.002
const _ALIGO_ITM_TRANSMITTANCE = 0.014
const _ALIGO_PRM_TRANSMITTANCE = 0.027
const _ALIGO_SRM_TRANSMITTANCE = 0.2
const _ALIGO_TEMPERATURE_SI = 290.0
const _ALIGO_THERMAL_SUSP_FREQ_SI = 9.0
const _ALIGO_THERMAL_SUSP_QUAL = 6e10
const _ALIGO_THERMAL_COAT_FREQ_SI = 1e4
const _ALIGO_THERMAL_COAT_QUAL = 6e6

# --- XLALSimNoisePSDSuspTherm: crude broadband suspension-thermal-noise power law,
# --- S_h(f) = L^-2 * (2*k*T / (pi^3 * f0^3 * M * Q)) * (f0/f)^5. ---
function _susp_therm(f::Float64, L::Float64, M::Float64, T::Float64, f0::Float64, Q::Float64)
    fac = 2.0 * _LAL_K_SI * T / (L^2 * M * Q * (_LAL_PI * f0)^3)
    return fac * (f0 / f)^5
end

# --- XLALSimNoisePSDMirrorTherm: crude broadband mirror/coating-thermal-noise power law,
# --- S_h(f) = L^-2 * (2*k*T / (pi^3 * f0^3 * M * Q)) * (f0/f). Same prefactor as _susp_therm,
# --- different (and shallower) frequency falloff. ---
function _mirror_therm(f::Float64, L::Float64, M::Float64, T::Float64, f0::Float64, Q::Float64)
    fac = 2.0 * _LAL_K_SI * T / (L^2 * M * Q * (_LAL_PI * f0)^3)
    return fac * (f0 / f)
end

# --- XLALSimNoisePSDaLIGOThermal: suspension + mirror-coating thermal noise at aLIGO's design
# --- parameters. ---
function _aligo_thermal(f::Float64)
    susp = _susp_therm(f, _ALIGO_ARMLENGTH_SI, _ALIGO_MIRROR_MASS_SI, _ALIGO_TEMPERATURE_SI,
                        _ALIGO_THERMAL_SUSP_FREQ_SI, _ALIGO_THERMAL_SUSP_QUAL)
    coat = _mirror_therm(f, _ALIGO_ARMLENGTH_SI, _ALIGO_MIRROR_MASS_SI, _ALIGO_TEMPERATURE_SI,
                          _ALIGO_THERMAL_COAT_FREQ_SI, _ALIGO_THERMAL_COAT_QUAL)
    return susp + coat
end

# --- XLALSimNoisePSDQuantum: quantum noise (shot noise + radiation pressure noise) of a
# --- signal-recycled interferometer, Buonanno & Chen, Phys. Rev. D 64, 062006 (2001) Eq. 5.13,
# --- as adapted in LALSimulation from GWINC's `shotrad.m` (with losses from mirrors only; BS/SRC
# --- coupling losses ignored, matching LALSimulation's own documented simplification). Ported
# --- verbatim, term for term, from the C source -- including its complex (COMPLEX16) intermediate
# --- coefficients C/D/P/Q/N, which become `Complex{Float64}` here. ---
function _quantum(f::Float64, I0::Float64, lambda::Float64, L::Float64, M::Float64, A::Float64,
                   A_BS::Float64, T_ITM::Float64, T_PRM::Float64, T_SRM::Float64, ds::Float64,
                   zeta::Float64, eta::Float64)
    Omega = 2.0 * _LAL_PI * f                         # [BC, table 1] signal angular frequency
    omega_0 = 2.0 * _LAL_PI * _LAL_C_SI / lambda       # [BC, table 1] laser angular frequency
    lambda_SR = A_BS
    lambda_PD = 1.0 - eta
    tau = sqrt(T_SRM)                                  # SRM transmittance (amplitude)
    rho = sqrt(1.0 - tau^2)                            # SRM reflectivity (amplitude)
    phi = (_LAL_PI - ds) / 2.0                         # [BC, between 2.14 & 2.15] SR detuning
    lambda_arm = 2.0 * A                               # [BC, after 5.2] round-trip loss in arm
    gamma_ac = T_ITM * _LAL_C_SI / (4.0 * L)           # [KLMTV-PRD2001] arm-cavity half bandwidth
    epsilon = lambda_arm / (2.0 * gamma_ac * L / _LAL_C_SI)  # [BC, after 5.2] arm-cavity loss coeff.
    r1 = sqrt(1.0 - T_ITM)                             # ITM reflectivity
    rarm = r1 - T_ITM * sqrt(1.0 - 2.0 * A) / (1.0 - r1 * sqrt(1.0 - 2.0 * A))  # arm-cavity reflectivity
    G_PRC = T_PRM / (1.0 + sqrt(1.0 - T_PRM) * rarm * sqrt(1.0 - A_BS))^2       # power-recycling gain
    I_0 = G_PRC * I0                                   # [BC, table 1] power at BS
    I_SQL = M * L^2 * gamma_ac^4 / (4.0 * omega_0)     # [BC, 2.14] power to reach free-mass SQL
    Kappa = 2.0 * ((I_0 / I_SQL) * gamma_ac^4) / (Omega^2 * (gamma_ac^2 + Omega^2))  # [BC, 2.13]
    beta = atan(Omega / gamma_ac)                      # [BnC, after 2.11] GW sideband phase shift
    h_SQL = sqrt(8.0 * _LAL_HBAR_SI / (M * (Omega * L)^2))  # [BnC, 2.12] SQL strain

    e2ib = exp(2.0im * beta)
    e4ib = exp(4.0im * beta)
    eneg2ib = exp(-2.0im * beta)
    eib = exp(1.0im * beta)
    enegib = exp(-1.0im * beta)

    # Coefficients [BC, Equations 5.8 to 5.12].
    C11_L = sqrt(1.0 - lambda_PD) * (
        (1.0 + rho^2) * (cos(2 * phi) + Kappa / 2.0 * sin(2 * phi)) - 2.0 * rho * cos(2 * beta) -
        1.0 / 4.0 * epsilon * (-2.0 * (1.0 + e2ib)^2 * rho +
                               4.0 * (1.0 + rho^2) * cos(beta)^2 * cos(2 * phi) +
                               (3.0 + e2ib) * Kappa * (1.0 + rho^2) * sin(2 * phi)) +
        lambda_SR * (e2ib * rho - 1.0 / 2.0 * (1.0 + rho^2) * (cos(2 * phi) + Kappa / 2.0 * sin(2 * phi))))
    C22_L = C11_L
    C12_L = sqrt(1.0 - lambda_PD) * tau^2 * (
        -(sin(2 * phi) + Kappa * sin(phi)^2) +
        1.0 / 2.0 * epsilon * sin(phi) * ((3.0 + e2ib) * Kappa * sin(phi) + 4.0 * cos(beta)^2 * cos(phi)) +
        1.0 / 2.0 * lambda_SR * (sin(2 * phi) + Kappa * sin(phi)^2))
    C21_L = sqrt(1.0 - lambda_PD) * tau^2 * (
        (sin(2 * phi) - Kappa * cos(phi)^2) +
        1.0 / 2.0 * epsilon * cos(phi) * ((3.0 + e2ib) * Kappa * sin(phi) - 4.0 * cos(beta)^2 * sin(phi)) +
        1.0 / 2.0 * lambda_SR * (-sin(2 * phi) + Kappa * cos(phi)^2))

    D1_L = sqrt(1.0 - lambda_PD) * (
        -(1.0 + rho * e2ib) * sin(phi) +
        1.0 / 4.0 * epsilon * (3.0 + rho + 2.0 * rho * e4ib + e2ib * (1.0 + 5.0 * rho)) * sin(phi) +
        1.0 / 2.0 * lambda_SR * e2ib * rho * sin(phi))
    D2_L = sqrt(1.0 - lambda_PD) * (
        -(-1.0 + rho * e2ib) * cos(phi) +
        1.0 / 4.0 * epsilon * (-3.0 + rho + 2.0 * rho * e4ib + e2ib * (-1.0 + 5.0 * rho)) * cos(phi) +
        1.0 / 2.0 * lambda_SR * e2ib * rho * cos(phi))

    P11 = 1.0 / 2.0 * sqrt(1.0 - lambda_PD) * sqrt(lambda_SR) * tau *
          (-2.0 * rho * e2ib + 2.0 * cos(2 * phi) + Kappa * sin(2 * phi))
    P22 = P11
    P12 = -sqrt(1.0 - lambda_PD) * sqrt(lambda_SR) * tau * sin(phi) * (2.0 * cos(phi) + Kappa * sin(phi))
    P21 = sqrt(1.0 - lambda_PD) * sqrt(lambda_SR) * tau * cos(phi) * (2.0 * sin(phi) - Kappa * cos(phi))

    Q11 = sqrt(lambda_PD) * (
        eneg2ib + rho^2 * e2ib - rho * (2.0 * cos(2 * phi) + Kappa * sin(2 * phi)) +
        1.0 / 2.0 * epsilon * rho * (eneg2ib * cos(2 * phi) +
                                     e2ib * (-2.0 * rho - 2.0 * rho * cos(2 * beta) + cos(2 * phi) +
                                             Kappa * sin(2 * phi)) +
                                     2.0 * cos(2 * phi) + 3.0 * Kappa * sin(2 * phi)) -
        1.0 / 2.0 * lambda_SR * rho * (2.0 * rho * e2ib - 2.0 * cos(2 * phi) - Kappa * sin(2 * phi)))
    Q22 = Q11
    Q12 = 0.0 + 0.0im
    Q21 = 0.0 + 0.0im

    N11 = sqrt(1.0 - lambda_PD) * sqrt(epsilon / 2.0) * tau * (
        Kappa * (1.0 + rho * e2ib) * sin(phi) +
        2.0 * cos(beta) * (enegib * cos(phi) - rho * eib * (cos(phi) + Kappa * sin(phi))))
    N22 = -sqrt(1.0 - lambda_PD) * sqrt(2.0 * epsilon) * tau * (-enegib + rho * eib) * cos(beta) * cos(phi)
    N12 = -sqrt(1.0 - lambda_PD) * sqrt(2.0 * epsilon) * tau * (enegib + rho * eib) * cos(beta) * sin(phi)
    N21 = sqrt(1.0 - lambda_PD) * sqrt(2.0 * epsilon) * tau * (
        -Kappa * (1.0 + rho) * cos(phi) + 2.0 * cos(beta) * (enegib + rho * eib) * cos(beta) * sin(phi))

    # >>>>>>>> QUANTUM NOISE POWER SPECTRAL DENSITY [BC, 5.13] <<<<<<<<<<<<<<<<<
    n = h_SQL^2 / (2.0 * Kappa * tau^2 * abs2(D1_L * sin(zeta) + D2_L * cos(zeta))) * (
        abs2(C11_L * sin(zeta) + C21_L * cos(zeta)) +
        abs2(C12_L * sin(zeta) + C22_L * cos(zeta)) +
        abs2(P11 * sin(zeta) + P21 * cos(zeta)) +
        abs2(P12 * sin(zeta) + P22 * cos(zeta)) +
        abs2(Q11 * sin(zeta) + Q21 * cos(zeta)) +
        abs2(Q12 * sin(zeta) + Q22 * cos(zeta)) +
        abs2(N11 * sin(zeta) + N21 * cos(zeta)) +
        abs2(N12 * sin(zeta) + N22 * cos(zeta)))

    return n
end

# --- XLALSimNoisePSDaLIGOQuantumZeroDetHighPower: the quantum-noise piece of the aLIGO
# --- "zero-detuning, high-power" configuration (no signal-recycling-cavity detuning, 125W input
# --- laser power). ---
function _aligo_quantum_zero_det_high_power(f::Float64)
    eta = 0.9                       # quantum efficiency of photodiode
    ds = 0.0                        # detuning phase (rad) -- "zero detuning"
    zeta = 116.0 * _LAL_PI_180      # homodyne detection phase (rad)
    return _quantum(f, _ALIGO_LASER_POWER_HIGH_SI, _ALIGO_LASER_WAVELENGTH_SI, _ALIGO_ARMLENGTH_SI,
                     _ALIGO_MIRROR_MASS_SI, _ALIGO_MIRROR_LOSS, _ALIGO_BS_LOSS,
                     _ALIGO_ITM_TRANSMITTANCE, _ALIGO_PRM_TRANSMITTANCE, _ALIGO_SRM_TRANSMITTANCE,
                     ds, zeta, eta)
end

# --- XLALSimNoisePSDaLIGOZeroDetHighPower: quantum + thermal noise, i.e. the full aLIGO
# --- zero-detuning-high-power design PSD (valid above ~9 Hz per LALSimulation's own documented
# --- caveat -- it does not include seismic/Newtonian-noise contributions). ---
function _aligo_zero_det_high_power_psd(f::Float64)
    return _aligo_quantum_zero_det_high_power(f) + _aligo_thermal(f)
end

# --- XLALSimNoisePSDAdvVirgo: phenomenological log-polynomial/exponential fit to the AdvVirgo
# --- design ASD, Eq. 6 of Manzotti & Dietz 2012 (arXiv:1202.4031). ---
function _avirgo_design_psd(f::Float64)
    x = log(f / 300.0)
    x2 = x^2
    asd = 1.259e-24 * (0.07 * exp(-0.142 - 1.437 * x + 0.407 * x2) +
                        3.1 * exp(-0.466 - 1.043 * x - 0.548 * x2) +
                        0.4 * exp(-0.304 + 2.896 * x - 0.293 * x2) +
                        0.09 * exp(1.466 + 3.722 * x - 0.984 * x2))
    return asd^2
end

# ================================================================================================
# Public API. Each curve is available as either:
#   - `<curve>_design_psd(f::AbstractVector{<:Real}) -> Vector{Float64}` -- evaluate at arbitrary
#     frequencies, returning a plain vector; or
#   - `<curve>_design_psd([f0::Real,] df::Real, n::Integer) -> FrequencySeries` -- evaluate on a
#     uniform grid of `n` bins (matching `FrequencySeries`'s own `(data, [f0,] df)` constructor
#     convention in src/GWData/Types.jl), returning a `FrequencySeries` directly.
#
# In the `FrequencySeries` form, the DC bin (f=0, only reachable when `f0=0.0`, the default) is
# set to exactly 0.0 rather than evaluated: both underlying formulas are singular or ill-defined
# at f=0 (the aLIGO quantum-noise formula divides by `Omega^2`; the AdvVirgo fit takes `log(0)`),
# and LALSimulation's own general-purpose PSD-array routine (`XLALSimNoisePSD` in
# LALSimNoisePSD.c) makes the identical choice ("kmin = 1; /* will set DC to zero */") rather than
# evaluate its `psdfunc` there. The plain-vector form does not add this guard -- a caller who
# explicitly asks for `f=0.0` gets whatever IEEE754 NaN/Inf the raw formula produces there, same
# as calling `lalsimulation.SimNoisePSDaLIGOZeroDetHighPower(0.0)` directly would.
# ================================================================================================

"""
    aligo_design_psd(f::AbstractVector{<:Real}) -> Vector{Float64}
    aligo_design_psd([f0::Real,] df::Real, n::Integer) -> FrequencySeries

Advanced LIGO "zero-detuning, high-power" design-sensitivity power spectral density (strain^2/Hz),
ported from LALSimulation's `XLALSimNoisePSDaLIGOZeroDetHighPower` (`lalsimulation/lib/
LALSimNoisePSD.c`) -- see this file's module docstring for the exact formula, its provenance
(Buonanno & Chen, Phys. Rev. D 64, 062006 (2001), quantum noise; plus a suspension/coating thermal
-noise fit), and an important naming caveat (this is NOT the simpler closed-form "Ajith fit" some
other codebases use under a similar name).

Valid above ~9 Hz (LALSimulation's own documented range for this model -- it omits seismic/
Newtonian-noise contributions that dominate the design curve below that).

See the "Public API" note above this docstring for the DC-bin (f=0) convention used by the
`FrequencySeries`-returning form.
"""
aligo_design_psd(f::AbstractVector{<:Real}) = Float64[_aligo_zero_det_high_power_psd(Float64(fi)) for fi in f]

function aligo_design_psd(f0::Real, df::Real, n::Integer)
    f0f, dff = Float64(f0), Float64(df)
    freqs = f0f .+ dff .* (0:(n - 1))
    data = [fi == 0.0 ? 0.0 : _aligo_zero_det_high_power_psd(fi) for fi in freqs]
    return FrequencySeries(data, f0f, dff; name="aLIGO zero-detuning high-power design PSD")
end
aligo_design_psd(df::Real, n::Integer) = aligo_design_psd(0.0, df, n)

"""
    aligo_design_asd(f::AbstractVector{<:Real}) -> Vector{Float64}
    aligo_design_asd([f0::Real,] df::Real, n::Integer) -> FrequencySeries

Advanced LIGO "zero-detuning, high-power" design-sensitivity amplitude spectral density
(strain/sqrt(Hz)), `sqrt.(aligo_design_psd(...))`. See [`aligo_design_psd`](@ref) for the source
and caveats.
"""
aligo_design_asd(f::AbstractVector{<:Real}) = sqrt.(aligo_design_psd(f))

function aligo_design_asd(f0::Real, df::Real, n::Integer)
    p = aligo_design_psd(f0, df, n)
    return FrequencySeries(sqrt.(p.data), p.f0, p.df; name="aLIGO zero-detuning high-power design ASD")
end
aligo_design_asd(df::Real, n::Integer) = aligo_design_asd(0.0, df, n)

"""
    avirgo_design_psd(f::AbstractVector{<:Real}) -> Vector{Float64}
    avirgo_design_psd([f0::Real,] df::Real, n::Integer) -> FrequencySeries

Advanced Virgo design-sensitivity power spectral density (strain^2/Hz), ported from
LALSimulation's `XLALSimNoisePSDAdvVirgo` (`lalsimulation/lib/LALSimNoisePSD.c`): a
phenomenological log-polynomial/exponential fit, Eq. 6 of Manzotti & Dietz 2012
(arXiv:1202.4031), to the AdvVirgo design curve published at
http://wwwcascina.virgo.infin.it/advirgo. See this file's module docstring for further
provenance detail.

See the "Public API" note above [`aligo_design_psd`](@ref) for the DC-bin (f=0) convention used
by the `FrequencySeries`-returning form.
"""
avirgo_design_psd(f::AbstractVector{<:Real}) = Float64[_avirgo_design_psd(Float64(fi)) for fi in f]

function avirgo_design_psd(f0::Real, df::Real, n::Integer)
    f0f, dff = Float64(f0), Float64(df)
    freqs = f0f .+ dff .* (0:(n - 1))
    data = [fi == 0.0 ? 0.0 : _avirgo_design_psd(fi) for fi in freqs]
    return FrequencySeries(data, f0f, dff; name="AdvVirgo design PSD")
end
avirgo_design_psd(df::Real, n::Integer) = avirgo_design_psd(0.0, df, n)

"""
    avirgo_design_asd(f::AbstractVector{<:Real}) -> Vector{Float64}
    avirgo_design_asd([f0::Real,] df::Real, n::Integer) -> FrequencySeries

Advanced Virgo design-sensitivity amplitude spectral density (strain/sqrt(Hz)),
`sqrt.(avirgo_design_psd(...))`. See [`avirgo_design_psd`](@ref) for the source.
"""
avirgo_design_asd(f::AbstractVector{<:Real}) = sqrt.(avirgo_design_psd(f))

function avirgo_design_asd(f0::Real, df::Real, n::Integer)
    p = avirgo_design_psd(f0, df, n)
    return FrequencySeries(sqrt.(p.data), p.f0, p.df; name="AdvVirgo design ASD")
end
avirgo_design_asd(df::Real, n::Integer) = avirgo_design_asd(0.0, df, n)

end # module ReferenceCurves
