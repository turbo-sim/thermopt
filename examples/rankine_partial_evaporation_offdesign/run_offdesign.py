"""Partial-evaporation ORC: design and sizing once, then off-design points on demand.

Run from this folder. In IPython:  %run run_offdesign.py   then   f(0.9, 5.0, 4.0)
"""

import copy

import numpy as np
import thermopt as th
from thermopt.cycles import cycle_power_simple_offdesign as od

# Step 1. Design
design = th.ThermodynamicCycleOptimization("case_novec649_design.yaml", out_dir="results")
design.config["solver_options"]["max_iterations"] = int(design.config["solver_options"]["max_iterations"])  # pysolver requires a Python int
design.load_config(design.config)
design.run_optimization()

# Step 2. Sizing
config = th.read_configuration_file("case_novec649_offdesign.yaml")
parameters = config["problem_formulation"]["fixed_parameters"]
hardware, design_point, design_controls = od.size_hardware(design.problem.cycle_data, parameters)
parameters["hardware"] = hardware
config["solver_options"]["max_iterations"] = int(config["solver_options"]["max_iterations"])

# Step 3. Variables: start from the design point, with bounds around it.
# The controls (pump speed, sink flow) are optimized; put them in parameters["controls"] to fix them.
c = design_point["components"]
p1, p3, h3 = float(c["compressor"]["state_in"].p), float(c["expander"]["state_in"].p), float(c["expander"]["state_in"].h)
_, h_liq, h_vap = od.saturation_quality(design_point["working_fluid"], p3, h3)
T_source_out_design = float(c["heater"]["hot_side"]["state_out"].T)
config["problem_formulation"]["design_variables"] = {
    "compressor_inlet_pressure": {"min": 0.4 * p1, "max": 2.5 * p1, "value": p1},
    "expander_inlet_pressure": {"min": 0.25 * p3, "max": 1.6 * p3, "value": p3},
    "expander_inlet_enthalpy": {"min": h3 - 1.2 * (h_vap - h_liq), "max": h3 + 1.2 * (h_vap - h_liq), "value": h3},
    "heat_source_exit_temperature": None,  # bounds depend on the operating point, set in offdesign_config
    "heat_sink_mass_flow": {"min": 0.3 * design_controls["heat_sink_mass_flow"],
                            "max": design_controls["heat_sink_mass_flow"],
                            "value": design_controls["heat_sink_mass_flow"]},
    "pump_speed": {"min": 0.4 * design_controls["pump_speed"], "max": 1.3 * design_controls["pump_speed"],
                   "value": design_controls["pump_speed"]},
}


def offdesign_config(mass_flow_ratio, delta_T_source, delta_T_sink):
    """Off-design configuration for the given changes relative to design."""
    cfg = copy.deepcopy(config)
    fp = cfg["problem_formulation"]["fixed_parameters"]
    fp["heat_source"]["mass_flow"] *= mass_flow_ratio
    fp["heat_source"]["inlet_temperature"] += delta_T_source
    fp["heat_sink"]["inlet_temperature"] += delta_T_sink
    T_source_in, T_sink_in = fp["heat_source"]["inlet_temperature"], fp["heat_sink"]["inlet_temperature"]
    fp["special_points"]["maximum_temperature"] = T_source_in
    fp["special_points"]["ambient_temperature"] = T_sink_in
    cfg["problem_formulation"]["design_variables"]["heat_source_exit_temperature"] = {
        "min": T_sink_in, "max": T_source_in - 0.1,
        "value": float(np.clip(T_source_out_design, T_sink_in, T_source_in - 0.1)),
    }
    return cfg


# Step 4. Off-design point (maximum net power), always solved from the design point.
# The design ThermoOpt object is reused: load_config replaces its problem with the off-design one.
offdesign = design


def f(mass_flow_ratio=1.0, delta_T_source=0.0, delta_T_sink=0.0):
    """Solve one off-design point and return the expander component.

    mass_flow_ratio: heat-source mass flow / design value
    delta_T_source:  heat-source inlet temperature - design value [K]
    delta_T_sink:    heat-sink (cooling water) inlet temperature - design value [K]
    """
    offdesign.load_config(offdesign_config(mass_flow_ratio, delta_T_source, delta_T_sink))
    offdesign.run_optimization()
    if not offdesign.solver.success:
        print(f"Warning: solver did not converge ({offdesign.solver.message})")
    return offdesign.problem.cycle_data["components"]["expander"]


if __name__ == "__main__":
    for name, args in (("design conditions", ()), ("off design", (0.9, 5.0, 4.0))):
        expander = f(*args)
        print(f"{name:18s} turbine inlet {float(expander['state_in'].p) / 1e5:.3f} bar, "
              f"quality {expander['inlet_quality']:.4f}, outlet {float(expander['state_out'].p) / 1e5:.4f} bar")
