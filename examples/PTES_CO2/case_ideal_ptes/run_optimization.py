import os
import pathlib
import thermopt as th

th.set_plot_options(fontsize=16)

CONFIG_FILE = "../case_PTES_CO2.yaml"
OUT_DIR_BASE = "results"
os.makedirs(OUT_DIR_BASE, exist_ok=True)

# Set up cycle and config
cycle = th.ThermodynamicCycleOptimization(CONFIG_FILE, out_dir=OUT_DIR_BASE)
cycle.set_config_value("solver_options.callbacks.plot_cycle", True)
# cycle.set_config_value("solver_options.max_iterations", 5)

# Set turbine efficiencies to 1.0
for var in [
    "problem_formulation.fixed_parameters.expander_charge.efficiency",
    "problem_formulation.fixed_parameters.expander_discharge.efficiency",
]:
    cycle.set_config_value(var, 1.0)

# Set compressor efficiencies to 1.0
for var in [
    "problem_formulation.fixed_parameters.compressor_charge.efficiency",
    "problem_formulation.fixed_parameters.compressor_discharge.efficiency",
]:
    cycle.set_config_value(var, 1.0)

# Set HX pinch point constraints to 0 °C
hx_constraints = [
    "$components.heater_charge.temperature_difference",
    "$components.cooler_charge.temperature_difference",
    "$components.recuperator_charge.temperature_difference",
    "$components.heater_discharge.temperature_difference",
    "$components.cooler_discharge.temperature_difference",
    "$components.recuperator_discharge.temperature_difference",
]

for var in hx_constraints:
    cycle.set_constraint(variable=var, type=">", value=0.0, normalize=10.0)

cycle.set_constraint(
    variable="$energy_analysis.hot_storage_upper_temperature - $energy_analysis.hot_storage_lower_temperature",
    type=">",
    value=0.0,
    normalize=False,
)

cycle.set_constraint(
    variable="$energy_analysis.cold_storage_upper_temperature - $energy_analysis.cold_storage_lower_temperature",
    type=">",
    value=0.0,
    normalize=False,
)

# Run optimization
cycle.run_optimization()
cycle.save_results()
cycle.problem.plot_cycle()
th.savefig_in_formats(cycle.problem.figure, "ideal_ptes_cycle")
