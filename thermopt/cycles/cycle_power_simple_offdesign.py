r"""Simple power cycle with partial evaporation, off design on fixed hardware.

States: 1 pump inlet, 2 pump outlet, 3 turbine inlet, 4 turbine outlet.

Variables: p1, p3, h3, T_h,out, plus the controls ``pump_speed`` and
``heat_sink_mass_flow`` (optimized if in ``design_variables``, fixed if in
``fixed_parameters.controls``).

Equalities (set to 0 in the YAML): ``heater.UA_residual``,
``cooler.UA_residual``, ``compressor.head_residual``,
``heater.energy_balance_residual``.
"""

import copy
import math

import jaxprop as jxp
import numpy as np
import pysolver_view as psv

from .. import utilities
from ..components import compression_process, compute_component_energy_flows, expansion_process, heat_exchanger
from ..components import components_offdesign as cod

COLORS_MATLAB = utilities.COLORS_MATLAB
CONTROL_NAMES = ("heat_sink_mass_flow", "pump_speed")
Z_NAMES = ("compressor_inlet_pressure", "expander_inlet_pressure", "expander_inlet_enthalpy",
           "heat_source_exit_temperature")


def _f(x):
    return float(np.asarray(x))


def saturation_quality(fluid, p, h):
    """Quality from saturation enthalpies (also outside the dome); returns (Q, h_l, h_v)."""
    if p >= _f(fluid.critical_point.p):
        raise psv.DomainError(f"Pressure {p:.4g} Pa at or above critical")
    hl = _f(fluid.get_state(jxp.PQ_INPUTS, p, 0.0).h)
    hv = _f(fluid.get_state(jxp.PQ_INPUTS, p, 1.0).h)
    return (h - hl) / (hv - hl), hl, hv


def _exchanger(fluid_hot, h_in_hot, h_out_hot, p_in_hot, p_out_hot,
               fluid_cold, h_in_cold, h_out_cold, p_in_cold, p_out_cold, spec):
    return heat_exchanger(fluid_hot, h_in_hot, h_out_hot, p_in_hot, p_out_hot,
                          fluid_cold, h_in_cold, h_out_cold, p_in_cold, p_out_cold,
                          counter_current=True, num_steps=int(spec["num_elements"]),
                          include_saturation_nodes=bool(spec.get("include_saturation_nodes", False)))


def evaluate_cycle(variables, parameters, constraints, objective_function):
    variables = copy.deepcopy(variables)
    parameters = copy.deepcopy(parameters)

    # Controls: optimized or fixed
    controls = dict(parameters.get("controls") or {})
    for name in CONTROL_NAMES:
        if name in variables:
            controls[name] = variables.pop(name)
    z = [variables.pop(name) for name in Z_NAMES]
    utilities.check_for_unused_keys(variables, "variables", raise_error=True)

    output = evaluate_offdesign(parameters, z, controls)
    f = utilities.evaluate_objective_function(output, objective_function)
    c_eq, c_ineq, constraint_report = utilities.evaluate_constraints(output, constraints)

    # Set colors for plotting
    orange, blue, red = COLORS_MATLAB[1], COLORS_MATLAB[0], COLORS_MATLAB[6]
    comps = output["components"]
    comps["heater"]["hot_side"]["plot_params"] = {"color": red, "linestyle": "-"}
    comps["heater"]["cold_side"]["plot_params"] = {"color": orange, "linestyle": "-"}
    comps["cooler"]["hot_side"]["plot_params"] = {"color": orange, "linestyle": "-"}
    comps["cooler"]["cold_side"]["plot_params"] = {"color": blue, "linestyle": "-"}
    comps["expander"]["plot_params"] = {"color": orange, "linestyle": "-"}
    comps["compressor"]["plot_params"] = {"color": orange, "linestyle": "-"}

    return {
        **output,
        "objective_function": f,
        "equality_constraints": c_eq,
        "inequality_constraints": c_ineq,
        "constraints_report": constraint_report,
    }


def evaluate_offdesign(parameters, z, controls, sizing=None):
    """Evaluate the cycle at z = [p1, p3, h3, T_h,out] on ``parameters["hardware"]``.

    With ``sizing={"working_fluid_mass_flow": m}`` (design point), z is
    [p1, p3, h3], the hardware is sized and returned under ``"hardware"``.
    """
    wf = jxp.Fluid(**parameters["working_fluid"], identifier="working_fluid")
    src = jxp.Fluid(**parameters["heating_fluid"], identifier="heating_fluid")
    snk = jxp.Fluid(**parameters["cooling_fluid"], identifier="cooling_fluid")
    source, sink = parameters["heat_source"], parameters["heat_sink"]
    heater_par, cooler_par, noz = parameters["heater"], parameters["cooler"], parameters["nozzle"]
    design = parameters.get("design", {})
    hw = parameters.get("hardware")

    z = [float(v) for v in z]
    if sizing:
        p1, p3, h3 = z[:3]
    else:
        p1, p3, h3, T_h_out = z
    m_h = float(source["mass_flow"])
    m_k, omega = float(controls["heat_sink_mass_flow"]), float(controls["pump_speed"])

    # Pressures
    p2 = p3 / (1.0 - float(heater_par["pressure_drop_cold_side"]))
    p4 = p1 / (1.0 - float(cooler_par["pressure_drop_hot_side"]))

    # State 1 (fixed subcooling) and state 3
    T_sat1 = _f(wf.get_state(jxp.PQ_INPUTS, p1, 0.0).T)
    st1 = wf.get_state(jxp.PT_INPUTS, p1, T_sat1 - float(parameters["closure"]["subcooling"]))
    h1 = _f(st1.h)
    st3 = wf.get_state(jxp.HmassP_INPUTS, h3, p3, generalize_quality=True)
    q3, _, _ = saturation_quality(wf, p3, h3)

    # Turbine map: velocity ratio relative to design x = nu / nu_d = c0_d / c0, with the
    # spouting velocity c0 = sqrt(2 dh_s) (the blade speed is fixed, so it cancels)
    st4s = wf.get_state(jxp.PSmass_INPUTS, p4, st3.s, supersaturation=True, generalize_quality=True)
    dh_s = h3 - _f(st4s.h)
    if dh_s <= 0.0:
        raise psv.DomainError(f"Nonpositive turbine isentropic enthalpy drop ({dh_s:.4g} J/kg)")
    c0 = math.sqrt(2.0 * dh_s)
    if sizing:
        c0_d = c0
        eta_td = float(design["turbine_efficiency"])
    else:
        c0_d = float(hw["turbine"]["design_spouting_velocity"])
        eta_td = float(hw["turbine"]["design_efficiency"])
    x = c0_d / c0
    eta_t = cod.turbine_efficiency(x, eta_td)
    if not 0.0 < x < 2.0:
        raise psv.DomainError(f"Turbine velocity ratio relative to design {x:.4g} outside (0, 2)")
    expander = expansion_process(wf, h3, p3, p4, eta_t, "isentropic")
    h4 = _f(expander["state_out"].h)

    # Equivalent nozzle
    flow = cod.nozzle_flow(wf, h3, _f(st3.s), p4, float(noz["search_minimum_temperature"]))
    Cd = float(noz["discharge_coefficient"])
    if sizing:
        if not flow["choked"]:
            raise psv.DomainError("Design point is not choked; the throat area cannot be sized")
        m_wf = float(sizing["working_fluid_mass_flow"])
        A_t = m_wf / (Cd * flow["G"])
    else:
        A_t = float(hw["nozzle"]["throat_area"])
        m_wf = Cd * A_t * flow["G"]

    # Pump map
    rho1 = _f(st1.d)
    Qv = m_wf / rho1
    H_req = (_f(wf.get_state(jxp.PSmass_INPUTS, p2, st1.s).h) - h1) / cod.GRAVITY
    if sizing:
        pd = design["pump"]
        pump = cod.PumpMap(float(pd["design_speed"]), Qv, H_req, float(pd["shutoff_head_ratio"]),
                           float(pd["design_efficiency"]))
    else:
        pump = cod.PumpMap.from_dict(hw["pump"])
    phi = pump.flow_coefficient(Qv, omega)
    if not 0.0 < phi < pump.phi0 or omega <= 0.0:
        raise psv.DomainError(f"Pump flow coefficient {phi:.4g} outside (0, phi0)")
    eta_p = pump.efficiency(phi)
    compressor = compression_process(wf, h1, p1, p2, eta_p, "isentropic")
    h2 = _f(compressor["state_out"].h)

    # Evaporator
    p_h_in = float(source["inlet_pressure"])
    p_h_out = p_h_in * (1.0 - float(heater_par["pressure_drop_hot_side"]))
    h_h_in = _f(src.get_state(jxp.PT_INPUTS, p_h_in, float(source["inlet_temperature"])).h)
    Q_e = m_wf * (h3 - h2)
    if Q_e <= 0.0:
        raise psv.DomainError("Nonpositive evaporator duty")
    if sizing:
        h_h_out = h_h_in - Q_e / m_h
    else:
        h_h_out = _f(src.get_state(jxp.PT_INPUTS, p_h_out, T_h_out).h)
    heater = _exchanger(src, h_h_in, h_h_out, p_h_in, p_h_out, wf, h2, h3, p2, p3, heater_par)
    ua_e = cod.integrate_ua(heater, m_wf)

    # Sink pump and condenser
    p_k_in = float(sink["inlet_pressure"])
    p_k_out = float(sink["exit_pressure"]) / (1.0 - float(cooler_par["pressure_drop_cold_side"]))
    h_k_in = _f(snk.get_state(jxp.PT_INPUTS, p_k_in, float(sink["inlet_temperature"])).h)
    heat_sink_pump = compression_process(snk, h_k_in, p_k_in, p_k_out,
                                         float(parameters["heat_sink_pump"]["efficiency"]), "isentropic")
    p_k_c_in, h_k_c_in = _f(heat_sink_pump["state_out"].p), _f(heat_sink_pump["state_out"].h)
    p_k_c_out = p_k_c_in * (1.0 - float(cooler_par["pressure_drop_cold_side"]))
    Q_c = m_wf * (h4 - h1)
    cooler = _exchanger(wf, h4, h1, p4, p1, snk, h_k_c_in, h_k_c_in + Q_c / m_k, p_k_c_in, p_k_c_out,
                        cooler_par)
    ua_c = cod.integrate_ua(cooler, m_k)

    # Source pump downstream of the evaporator
    heat_source_pump = compression_process(src, h_h_out, p_h_out, float(source["exit_pressure"]),
                                           float(parameters["heat_source_pump"]["efficiency"]), "isentropic")

    # Hardware (sized at design) and available conductance
    if sizing:
        hw = {
            "evaporator": {"UA_design": ua_e["UA"], "ua_exponent": float(design["evaporator"]["ua_exponent"]),
                           "source_mass_flow_design": m_h},
            "condenser": {"UA_design": ua_c["UA"], "ua_exponent": float(design["condenser"]["ua_exponent"]),
                          "sink_mass_flow_design": m_k},
            "nozzle": {"throat_area": A_t},
            "turbine": {"design_spouting_velocity": c0_d, "design_efficiency": eta_td},
            "pump": pump.to_dict(),
        }
    UA_e_av = float(hw["evaporator"]["UA_design"]) * (
        m_h / float(hw["evaporator"]["source_mass_flow_design"])) ** float(hw["evaporator"]["ua_exponent"])
    UA_c_av = float(hw["condenser"]["UA_design"]) * (
        m_k / float(hw["condenser"]["sink_mass_flow_design"])) ** float(hw["condenser"]["ua_exponent"])
    Q_ref = float(hw.get("reference_duty", Q_e))
    dh_ref = Q_ref / float(hw.get("reference_working_fluid_mass_flow", m_wf))

    # Components
    heater["hot_side"]["mass_flow"], heater["cold_side"]["mass_flow"] = m_h, m_wf
    cooler["hot_side"]["mass_flow"], cooler["cold_side"]["mass_flow"] = m_wf, m_k
    expander["mass_flow"] = compressor["mass_flow"] = m_wf
    heat_source_pump["mass_flow"], heat_sink_pump["mass_flow"] = m_h, m_k
    heater.update({
        "UA_required": ua_e["UA"], "UA_available": UA_e_av,
        "UA_residual": ua_e["ln_UA"] - math.log(UA_e_av),
        "energy_balance_residual": (m_h * (h_h_in - h_h_out) - Q_e) / Q_ref,
        "specific_duty_ratio": (h3 - h2) / dh_ref,
        "min_approach": ua_e["min_approach"],
    })
    cooler.update({
        "UA_required": ua_c["UA"], "UA_available": UA_c_av,
        "UA_residual": ua_c["ln_UA"] - math.log(UA_c_av),
        "min_approach": ua_c["min_approach"],
    })
    expander.update({
        "inlet_quality": q3, "velocity_ratio_relative": x, "isentropic_efficiency": eta_t,
        "nozzle": {"throat_area": A_t, "mass_flux": flow["G"], "choked": flow["choked"],
                   "critical_pressure": flow["critical_pressure"]},
    })
    compressor.update({
        "speed": omega, "flow_coefficient": phi, "isentropic_efficiency": eta_p,
        "head_residual": (H_req - pump.head(Qv, omega)) / pump.H_d,
    })
    components = {
        "expander": expander, "compressor": compressor, "heater": heater, "cooler": cooler,
        "heat_source_pump": heat_source_pump, "heat_sink_pump": heat_sink_pump,
    }
    compute_component_energy_flows(components)

    # Energy analysis
    W_t, W_p = expander["power"], compressor["power"]
    W_aux = heat_source_pump["power"] + heat_sink_pump["power"]
    W_net = W_t - W_p - W_aux
    h_min = _f(src.get_state(jxp.PT_INPUTS, float(source["exit_pressure"]), float(source["minimum_temperature"])).h)
    Q_max = m_h * (h_h_in - h_min)
    energy_analysis = {
        "heater_heat_flow": Q_e,
        "heater_heat_flow_max": Q_max,
        "cooler_heat_flow": Q_c,
        "expander_power": W_t,
        "compressor_power": W_p,
        "heat_source_pump_power": heat_source_pump["power"],
        "heat_sink_pump_power": heat_sink_pump["power"],
        "net_cycle_power": W_net,
        "net_system_power": W_net,
        "net_system_power_normalized": W_net / float(hw.get("reference_net_power", W_net)),
        "mass_flow_heating_fluid": m_h,
        "mass_flow_working_fluid": m_wf,
        "mass_flow_cooling_fluid": m_k,
        "cycle_efficiency": W_net / Q_e,
        "system_efficiency": W_net / Q_max,
        "backwork_ratio": (W_p + W_aux) / W_t,
        "energy_balance": (Q_e + W_p) - (W_t + Q_c),
    }

    output = {
        "components": components,
        "energy_analysis": energy_analysis,
        "working_fluid": wf,
        "heating_fluid": src,
        "cooling_fluid": snk,
    }
    if sizing:
        hw["reference_duty"] = Q_e
        hw["reference_working_fluid_mass_flow"] = m_wf
        hw["reference_net_power"] = W_net
        output["hardware"] = hw
    return output


def size_hardware(design_output, parameters):
    """Size the hardware from a converged ``simple`` design output.

    Flows are scaled to ``heat_source.mass_flow``. Returns
    (hardware, design_output, design_controls).
    """
    comps = design_output["components"]
    ea = design_output["energy_analysis"]
    lam = float(parameters["heat_source"]["mass_flow"]) / float(ea["mass_flow_heating_fluid"])
    controls = {"heat_sink_mass_flow": lam * float(ea["mass_flow_cooling_fluid"]),
                "pump_speed": float(parameters["design"]["pump"]["design_speed"])}
    z = [_f(comps["compressor"]["state_in"].p), _f(comps["expander"]["state_in"].p),
         _f(comps["expander"]["state_in"].h)]
    sizing = {"working_fluid_mass_flow": lam * float(ea["mass_flow_working_fluid"])}
    output = evaluate_offdesign({k: v for k, v in parameters.items() if k != "hardware"}, z, controls, sizing)
    return output["hardware"], output, controls
