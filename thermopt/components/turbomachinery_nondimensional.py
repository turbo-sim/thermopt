"""Design-point models for centrifugal compressors and radial inflow turbines.

Reference: Parisi et al., Advancing non-dimensional models for radial
turbomachinery applied to pumped thermal energy storage, Energy Conversion
and Management 342 (2025), 120085, doi:10.1016/j.enconman.2025.120085.

States use jaxprop names (T, p, h, s, rho, mu, a) and SI units. Performance
functions evaluate trial states without a machine object. They use NumPy by
default; pass xp=jax.numpy to trace their algebra. Turbines additionally need
an h,s property flash, so tracing requires a compatible jaxprop fluid. The
SciPy outlet solvers are separate host functions, not differentiable kernels.

These are single-phase design correlations, not off-design performance maps.
Compressor correlation efficiency is polytropic; turbine efficiency is overall
isentropic. The models are not validated for wet flow.
"""

import numpy as np
import jaxprop as props
from scipy.optimize import brentq

__all__ = [
    "compressor_defaults", "turbine_defaults", "friction_factor",
    "polytropic_head", "turbine_baseline_efficiency",
    "compressor_performance", "turbine_performance",
    "compressor_reference_coefficients", "compressor_mach_correction",
    "compressor_relative_velocity_ratio", "compressor_diameter_to_chord_ratio",
    "turbine_kinetic_energy_fraction",
    "solve_compressor", "solve_turbine",
]


def compressor_defaults():
    """Return fresh compressor design parameters (Section 3.2).

    Flow coefficient uses Eq. (3), based on total inlet density. Backsweep is
    in degrees, roughness in m, and metal density in kg/m³. Clearance is the
    tip gap divided by exit blade height. d_eta_other is an additive efficiency
    correction (fraction, not percent). is_shrouded selects the impeller;
    has_vaned selects a vaned rather than vaneless diffuser. The baseline is
    an unshrouded, 40-degree backswept impeller with a vaned diffuser.
    """
    return dict(flow_coefficient=0.07, backsweep=40.0, is_shrouded=False,
                has_vaned=True, roughness=3.6e-6, clearance_height_ratio=0.02,
                metal_density=7800.0, d_eta_other=0.0)


def turbine_defaults():
    """Return fresh turbine design parameters (Section 3.3).

    specific_speed is dimensionless and uses rotor-exit static density in
    Eq. (4), not inlet or overall-outlet density; the study spans 0.35–0.85.
    is_scalloped selects a scalloped disk. diffuser_loss_fraction is the
    fraction of rotor-exit kinetic energy lost downstream. Clearance is the
    tip gap/blade-height ratio, metal_density is in kg/m³, and d_eta_other is
    an additive rotor-efficiency correction expressed as a fraction.
    """
    return dict(specific_speed=0.6, is_scalloped=True,
                diffuser_loss_fraction=0.35, clearance_height_ratio=0.02,
                metal_density=7800.0, d_eta_other=0.0)


def friction_factor(reynolds, relative_roughness, *, xp=np):
    """Return the Darcy friction factor for the compressor Reynolds correction.

    reynolds is based on blade chord and relative inlet velocity; roughness
    is sand-grain height/chord, approximated by average roughness/chord here.
    The laminar and turbulent expressions are blended smoothly following the
    Casey–Robinson procedure cited beside Eqs. (21)–(22). The factor of four
    converts the blended skin-friction coefficient to a Darcy factor.
    Inputs need positive Reynolds number and nonnegative roughness.
    xp may be numpy or jax.numpy; broadcastable array inputs also work.
    """
    laminar = 2.656 / xp.sqrt(reynolds)
    turbulent = 0.136 / (-xp.log10(0.2 * relative_roughness + 12.5 / reynolds))**2.15
    laminar_weight = 1 / (1 + xp.exp(-5 * (laminar / turbulent - 1)))
    return 4 * (laminar_weight * laminar + (1 - laminar_weight) * turbulent)


def polytropic_head(inlet, outlet, *, xp=np):
    """Return signed Mallen–Saville polytropic head in J/kg (Eq. 9).

    States or mappings contain canonical temperature [K], enthalpy [J/kg] and
    entropy [J/(kg K)] fields. Compression
    gives positive head and expansion negative head. The entropy correction
    uses the logarithmic mean of endpoint temperatures. log1p and its series
    limit keep values and derivatives finite when temperatures coincide.
    This approximation is exact for a calorically perfect gas; it is not a
    path integration. xp selects NumPy or JAX without converting tracers.
    """
    change = (outlet["temperature"] - inlet["temperature"]) / inlet["temperature"]
    near_equal = xp.abs(change) < 1e-4
    # Prevent 0/0 even in the inactive branch: JAX differentiates both branches.
    safe_change = xp.where(near_equal, 1.0, change)
    series = 1 + change / 2 - change**2 / 12 + change**3 / 24
    log_mean = inlet["temperature"] * xp.where(near_equal, series, safe_change / xp.log1p(safe_change))
    return outlet["enthalpy"] - inlet["enthalpy"] - (outlet["entropy"] - inlet["entropy"]) * log_mean


def compressor_reference_coefficients(flow_coefficient, *, is_shrouded=False, xp=np):
    """Return reference work and vaned polytropic-head coefficients (Eqs. 12-13).

    These dimensionless coefficients describe a 40-degree backswept impeller
    before Mach, Reynolds, clearance or diffuser corrections. Table 4 selects
    shrouded/unshrouded coefficients. Their ratio is the baseline polytropic
    efficiency in Eq. (14), as plotted in Figure 1. The work coefficient is
    the reference curve in Figure 2. flow_coefficient follows Eq. (3).
    Scalar and array inputs are supported; xp selects NumPy or JAX.
    """
    flow = flow_coefficient
    work = (xp.where(is_shrouded, 0.62, 0.68)
            + (flow / (-np.pi / 4 * xp.where(is_shrouded, 0.4, 0.37)))**3
            + np.pi / 4 * xp.where(is_shrouded, 0.0014, 0.002) / flow)
    head = (xp.where(is_shrouded, 0.51, 0.59)
            + 4 / np.pi * xp.where(is_shrouded, 1.0, 0.7) * flow
            - (4 / np.pi)**2 * xp.where(is_shrouded, 7.6, 7.5) * flow**2
            - np.pi / 4 * 0.00025 / flow)
    return work, head


def compressor_mach_correction(flow_coefficient, peripheral_mach, *, xp=np):
    """Return the additive efficiency penalty for high tip Mach (Eqs. 18-20).

    peripheral_mach is tip speed divided by the total-inlet speed of sound.
    The correction is zero up to Mach 0.8 and negative above it. Figure 3
    illustrates this penalty at different pressure ratios. Both arguments
    are dimensionless and may be broadcastable NumPy or JAX arrays.
    """
    excess = xp.maximum(0.0, flow_coefficient * (peripheral_mach - 0.8))
    return -(0.05 + 3 * excess) * excess


def compressor_relative_velocity_ratio(flow_coefficient):
    """Return relative inlet velocity / tip speed using the Figure 4 correlation.

    This dimensionless ratio sets the chord Reynolds number in Eq. (22).
    Parisi et al. (2025), doi:10.1016/j.enconman.2025.120085, contains an
    error in the printed Eq. (23): 0.42 + 5*phi - 14*phi**2 disagrees with
    Figure 4 and its mean-line points. Use 0.38 + 3*phi - 7.8*phi**2, which
    reproduces the plotted correlation within PDF coordinate precision.
    The fitted flow-coefficient range is 0.01-0.15. The original thermopt
    implementation followed the erroneous printed equation.

    flow_coefficient follows Eq. (3); the returned w/U ratio is dimensionless.
    NumPy and JAX scalar/array inputs work directly through arithmetic only.
    Both actual and reference Reynolds numbers must use this same ratio.
    """
    return 0.38 + 3 * flow_coefficient - 7.8 * flow_coefficient**2


def compressor_diameter_to_chord_ratio(flow_coefficient, backsweep=40.0, *, xp=np):
    """Return impeller diameter divided by blade chord (Eqs. 24-25).

    The reference ratio is 1.1 + 22*phi at 40-degree backsweep. Smaller
    backsweep angles [degrees] increase the ratio; no further correction is
    applied above 40 degrees. Figure 5 compares mean-line designs corrected
    to this reference. The dimensionless ratio sets the chord used in the
    Reynolds correction. xp supports NumPy arrays or JAX transformations.
    """
    return 1.1 + 22 * flow_coefficient + 0.01 * xp.maximum(0.0, 40.0 - backsweep)


def turbine_kinetic_energy_fraction(specific_speed, *, xp=np):
    """Return rotor-exit kinetic energy divided by isentropic rotor work (Eq. 36).

    specific_speed is based on the static rotor outlet. This dimensionless
    kinetic fraction relates baseline total-to-total and total-to-static
    efficiencies through eta_ts = eta_tt / (1 + kappa), Eq. (37), shown in
    Figure 7. With a diffuser, only its loss fraction times kappa contributes
    to the overall efficiency penalty (Eq. 39). xp selects NumPy or JAX.
    """
    return xp.maximum(0.02 + 0.08 * specific_speed**2, 0.2 * specific_speed**4)


def compressor_performance(inlet, outlet, outlet_isentropic, mass_flow, design=None, *, xp=np):
    """Evaluate compressor efficiency, size, speed and stress at a trial outlet.

    Implements Section 3.2, Eqs. (12)–(27) and Table 4. inlet/outlet are total
    states with canonical JaxProp property names; outlet_isentropic is at the same outlet pressure
    and inlet entropy. mass_flow is kg/s. design partially overrides
    compressor_defaults(); inputs are not modified. The paper's
    flow-coefficient data span approximately 0.01–0.15.

    Returns correlation/polytropic/isentropic efficiencies as fractions,
    positive consumed power [W], impeller diameter [m], peripheral speed
    [m/s], angular speed [rad/s] and disk stress [Pa]. The solved residual is
    correlation_efficiency minus polytropic_efficiency (Eq. 11). Appendix B
    supplies stress factors: shrouded stress interpolates between flow
    coefficients 0.05 and 0.12 and clips outside those endpoints.

    No property flashes or scalar validation occur here. Pass xp=jax.numpy
    to jit/grad/vmap the algebra, including continuous design parameters.
    """
    design = compressor_defaults() | ({} if design is None else design)
    flow = design["flow_coefficient"]
    shrouded = design["is_shrouded"]
    # Table 4: work and polytropic-head coefficients at 40-degree backsweep.
    work_reference, head_reference = compressor_reference_coefficients(
        flow, is_shrouded=shrouded, xp=xp,
    )
    vaned_efficiency = head_reference / work_reference
    vaneless_penalty = 0.017 / (0.04 + 20 / np.pi * flow + vaned_efficiency**3)
    baseline = vaned_efficiency - xp.where(design["has_vaned"], 0.0, vaneless_penalty)

    backsweep_deficit = xp.maximum(0.0, 40.0 - design["backsweep"])
    work_coefficient = work_reference + (0.004 + 0.2 * flow**2) * backsweep_deficit
    work = outlet["enthalpy"] - inlet["enthalpy"]
    tip_speed = xp.sqrt(work / work_coefficient)
    diameter = xp.sqrt(mass_flow / (inlet["density"] * flow * tip_speed))
    mach_correction = compressor_mach_correction(flow, tip_speed / inlet["speed_of_sound"], xp=xp)

    # Eqs. (22)–(25): chord and relative inlet speed set Reynolds number.
    diameter_over_chord = compressor_diameter_to_chord_ratio(flow, design["backsweep"], xp=xp)
    chord = diameter / diameter_over_chord
    velocity_ratio = compressor_relative_velocity_ratio(flow)
    relative_speed = velocity_ratio * tip_speed
    viscosity = inlet["viscosity"] / inlet["density"]
    reynolds = relative_speed * chord / viscosity
    reference_chord = 0.4 / diameter_over_chord
    reference_reynolds = velocity_ratio * 277 * reference_chord / 15.7e-6
    friction = friction_factor(reynolds, design["roughness"] / chord, xp=xp)
    reference_friction = friction_factor(reference_reynolds, 3.6e-6 / reference_chord, xp=xp)
    reynolds_correction = (
        -(0.05 + 0.002 / (flow + 0.0025))
        * (friction - reference_friction) / reference_friction
    )

    correlation_efficiency = (baseline - 0.0006 * backsweep_deficit + mach_correction
                              + reynolds_correction - 0.5 * (design["clearance_height_ratio"] - 0.02)
                              + design["d_eta_other"])
    shrouded_stress_factor = 0.67 + (0.75 - 0.67) * (xp.clip(flow, 0.05, 0.12) - 0.05) / (0.12 - 0.05)
    stress_factor = xp.where(shrouded, shrouded_stress_factor, 0.41)
    return dict(correlation_efficiency=correlation_efficiency,
                polytropic_efficiency=polytropic_head(inlet, outlet, xp=xp) / work,
                isentropic_efficiency=(outlet_isentropic["enthalpy"] - inlet["enthalpy"]) / work,
                power=mass_flow * work, diameter=diameter, peripheral_speed=tip_speed,
                angular_speed=2 * tip_speed / diameter,
                mechanical_stress=stress_factor * design["metal_density"] * tip_speed**2)


def turbine_baseline_efficiency(specific_speed, *, xp=np):
    """Return baseline rotor isentropic efficiency (Eq. 29 and Figure 6).

    specific_speed is dimensionless and defined at the rotor static outlet
    (Eq. 4). The extra quadratic penalty below 0.45 models low-speed
    departures. This baseline precedes scalloping, clearance, diffuser and
    Reynolds corrections; it is not final machine efficiency. The fitted
    study range is 0.35–0.85. xp supports NumPy arrays or JAX tracing.
    """
    deviation = specific_speed - 0.6
    # Keep 0.5 to match the published Figure 6 and the authors' implementation.
    # Table 6 contains a typo: it lists 0.30 for this cubic coefficient.
    # Parisi et al. (2025), doi:10.1016/j.enconman.2025.120085.
    baseline = 0.93 - 0.3 * deviation**2 - 0.5 * deviation**3
    return baseline - 4.2 * xp.minimum(0.0, specific_speed - 0.45)**2


def turbine_performance(fluid, inlet, outlet, outlet_isentropic, mass_flow, design=None, *, xp=np, errors="raise"):
    """Evaluate a radial turbine at a trial total outlet (Section 3.3).

    Implements Eqs. (29)–(41), (46)–(48) and Appendix B. States use jaxprop
    names and SI units; mass_flow is kg/s. design partially overrides
    turbine_defaults(). The supplied fluid performs the rotor static h,s
    flash: this density is essential to specific-speed sizing in Eq. (4).
    Neither state mappings nor design parameters are modified.

    Returns the same eight performance quantities/units as the compressor,
    with positive produced power [W]. Here correlation_efficiency is overall
    isentropic efficiency; its difference from isentropic_efficiency is the
    solved residual (Eq. 10). The rotor entropy allocation follows the
    paper's diffuser-loss model, not an additional empirical loss.

    xp=jax.numpy enables tracing if fluid.get_state supports h,s under JAX
    (e.g. FluidPerfectGas). A HEOS fluid still requires a host flash; changing
    xp alone cannot make that backend differentiable. For a native model under
    JAX transformations, select errors="mask" explicitly for the rotor flash.
    """
    design = turbine_defaults() | ({} if design is None else design)
    speed = design["specific_speed"]
    baseline = turbine_baseline_efficiency(speed, xp=xp)
    kinetic_fraction = turbine_kinetic_energy_fraction(speed, xp=xp)
    clearance_correction = -1.15 * baseline * (design["clearance_height_ratio"] - 0.02)
    scalloping_correction = xp.where(design["is_scalloped"], -0.03, 0.0)
    rotor_efficiency = baseline + clearance_correction + scalloping_correction + design["d_eta_other"]
    diffuser_factor = 1 + design["diffuser_loss_fraction"] * kinetic_fraction
    overall_efficiency = rotor_efficiency / diffuser_factor

    work = inlet["enthalpy"] - outlet["enthalpy"]
    isentropic_work = inlet["enthalpy"] - outlet_isentropic["enthalpy"]
    rotor_isentropic_work = isentropic_work / diffuser_factor
    rotor_entropy = (
        inlet["entropy"] + (outlet["entropy"] - inlet["entropy"])
        * (1 - rotor_efficiency) / (diffuser_factor - rotor_efficiency)
    )
    rotor_enthalpy = outlet["enthalpy"] - kinetic_fraction * rotor_isentropic_work
    rotor_options = {} if errors == "raise" else {"errors": errors}
    rotor_static = fluid.get_state(props.HmassSmass_INPUTS, rotor_enthalpy, rotor_entropy, **rotor_options)
    tip_speed = 0.7 * xp.sqrt(2 * rotor_isentropic_work)
    angular_speed = speed * rotor_isentropic_work**0.75 / xp.sqrt(mass_flow / rotor_static["density"])
    diameter = 2 * tip_speed / angular_speed
    reynolds = 2 * mass_flow / (inlet["viscosity"] * diameter)
    correlation_efficiency = 1 - (1 - overall_efficiency) * (
        0.3 + 0.7 * xp.minimum(1.0, reynolds / 352000.0)**-0.2
    )
    stress_factor = xp.where(design["is_scalloped"], 0.15, 0.35)
    return dict(correlation_efficiency=correlation_efficiency,
                polytropic_efficiency=-work / polytropic_head(inlet, outlet, xp=xp),
                isentropic_efficiency=work / isentropic_work,
                power=mass_flow * work, diameter=diameter, peripheral_speed=tip_speed,
                angular_speed=angular_speed,
                mechanical_stress=stress_factor * design["metal_density"] * tip_speed**2)


def _validate_inputs(inlet, p_out, mass_flow, design, defaults, *, compressor):
    """Reject invalid scalar host inputs before starting an outlet root solve.

    Validate finite states, positive pressure/temperature/flow, pressure
    direction and design parameter names/domains. Published study ranges are
    documented, not clamped: extrapolation remains the caller's responsibility.
    Validation stays outside arithmetic functions so it never inspects JAX
    tracers or adds Python control flow to a differentiated kernel.
    """
    for key in ("temperature", "pressure", "enthalpy", "entropy", "density", "viscosity", "speed_of_sound"):
        value = float(inlet[key])
        if not np.isfinite(value) or (key not in ("enthalpy", "entropy") and value <= 0):
            raise ValueError(f"Inlet {key} must be finite" + (" and positive." if key not in ("enthalpy", "entropy") else "."))
    for key, value in (("mass_flow", mass_flow), ("p_out", p_out)):
        if value is None or not np.isfinite(value) or value <= 0:
            raise ValueError(f"{key} must be finite and positive.")
    equal_pressure = np.isclose(p_out, inlet["pressure"], rtol=1e-12, atol=0)
    if equal_pressure or (compressor and p_out < inlet["pressure"]) or (not compressor and p_out > inlet["pressure"]):
        direction = "greater" if compressor else "less"
        raise ValueError(f"Outlet pressure must be {direction} than inlet pressure.")
    unknown = design.keys() - defaults.keys()
    if unknown:
        raise ValueError(f"Unknown design parameters: {', '.join(sorted(unknown))}")
    for key, value in design.items():
        if isinstance(defaults[key], bool):
            if not isinstance(value, (bool, np.bool_)):
                raise ValueError(f"{key} must be a boolean.")
        elif not np.isfinite(value):
            raise ValueError(f"{key} must be finite.")
    for key in ("roughness", "clearance_height_ratio", "backsweep"):
        if key in design and design[key] < 0:
            raise ValueError(f"{key} must be nonnegative.")
    for key in ("metal_density", "specific_speed"):
        if key in design and design[key] <= 0:
            raise ValueError(f"{key} must be positive.")
    if compressor and design["flow_coefficient"] < 1e-4:
        raise ValueError("flow_coefficient must be at least 1e-4.")
    if not compressor and not 0 <= design["diffuser_loss_fraction"] <= 1:
        raise ValueError("diffuser_loss_fraction must be between 0 and 1.")


def _solve_outlet(fluid, inlet, p_out, mass_flow, design, *, compressor):
    """Match correlation and state efficiencies in outlet enthalpy (Section 4).

    Cache the inlet and isentropic outlet, then use a bracketed Brent solve
    over outlet enthalpy [J/kg] at fixed pressure [Pa]. Original brackets
    correspond to isentropic efficiencies of 1 to 0.5 for compression and
    1 to 0.1 for expansion. Failed bracketing raises a contextual ValueError;
    an unconverged or nonphysical root is never returned as a valid design.

    Returns jaxprop states, performance, efficiency residual, the completed
    design and SciPy iteration/function-call diagnostics. All flashes use the
    explicit fluid argument. This host routine is not JIT-able or
    differentiable: that requires performance kernels and a differentiable
    root solver. It does not mutate or retain intermediate machine states.
    """
    defaults = compressor_defaults() if compressor else turbine_defaults()
    design = defaults | ({} if design is None else design)
    _validate_inputs(inlet, p_out, mass_flow, design, defaults, compressor=compressor)
    isentropic = fluid.get_state(props.PSmass_INPUTS, p_out, inlet["entropy"])
    isentropic_work = float(isentropic["enthalpy"] - inlet["enthalpy"])
    low_efficiency = 0.5 if compressor else 0.1
    other_bound = float(inlet["enthalpy"]) + isentropic_work * (1 / low_efficiency if compressor else low_efficiency)
    efficiency_key = "polytropic_efficiency" if compressor else "isentropic_efficiency"

    def evaluate(enthalpy):
        """Flash one trial h,p outlet and return its state and performance."""
        outlet = fluid.get_state(props.HmassP_INPUTS, enthalpy, p_out)
        if compressor:
            performance = compressor_performance(inlet, outlet, isentropic, mass_flow, design)
        else:
            performance = turbine_performance(fluid, inlet, outlet, isentropic, mass_flow, design)
        return outlet, performance

    def residual(enthalpy):
        """Return the scalar correlation-minus-thermodynamic efficiency mismatch."""
        _, performance = evaluate(enthalpy)
        return float(performance["correlation_efficiency"] - performance[efficiency_key])

    kind = "compressor" if compressor else "turbine"
    try:
        enthalpy, status = brentq(residual, float(isentropic["enthalpy"]), other_bound, full_output=True)
    except ValueError as error:
        raise ValueError(
            f"Could not solve {kind} at p_out={p_out:g} Pa, "
            f"mass_flow={mass_flow:g} kg/s: {error}"
        ) from error
    outlet, performance = evaluate(enthalpy)
    performance = {key: float(value) for key, value in performance.items()}
    final_residual = performance["correlation_efficiency"] - performance[efficiency_key]
    if (not status.converged or abs(final_residual) > 1e-8
            or not all(np.isfinite(value) for value in performance.values())
            or not 0 < performance[efficiency_key] <= 1):
        raise RuntimeError(f"Invalid {kind} solution: efficiency residual={final_residual:g}.")
    return dict(state_in=inlet, state_out=outlet, state_out_isentropic=isentropic,
                performance=performance, design=design, residual=final_residual,
                solver=dict(converged=status.converged, iterations=status.iterations,
                            function_calls=status.function_calls))


def solve_compressor(fluid, inlet, p_out, mass_flow, design=None):
    """Solve a centrifugal-compressor design point from an existing inlet state.

    fluid is a jaxprop model providing p,s and h,p flashes; inlet is its total
    inlet FluidState. p_out [Pa] must exceed inlet pressure and mass_flow
    [kg/s] must be positive. design partially overrides compressor_defaults().
    Section 4's efficiency matching uses only outlet enthalpy; there is no
    temperature-iteration option or stored machine state.

    Returns a dictionary with state_in, state_out, state_out_isentropic,
    performance (units in compressor_performance), completed design, residual
    and solver diagnostics. This is a scalar SciPy host solve; use
    compressor_performance for JAX transformations of the correlation itself.
    """
    return _solve_outlet(fluid, inlet, p_out, mass_flow, design, compressor=True)


def solve_turbine(fluid, inlet, p_out, mass_flow, design=None):
    """Solve a radial-inflow-turbine design point from an existing inlet state.

    fluid is a jaxprop model providing p,s, h,p and h,s flashes; inlet is its
    total inlet FluidState. p_out [Pa] must be below inlet pressure and
    mass_flow [kg/s] must be positive. design partially overrides
    turbine_defaults(). The rotor static flash sets size/speed; state_out is
    the overall total outlet downstream of the diffuser.

    Returns the same fields as solve_compressor, with positive produced
    power and overall isentropic efficiency. Section 4's scalar efficiency
    matching uses outlet enthalpy and SciPy, outside the JAX kernel boundary.
    """
    return _solve_outlet(fluid, inlet, p_out, mass_flow, design, compressor=False)
