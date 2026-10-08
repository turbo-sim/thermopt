"""Compare fixed superheat and partial evaporation using ThermoOpt's native workflow."""

from copy import deepcopy
from functools import lru_cache
import json
from pathlib import Path
import shutil
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import jaxprop as props
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
import thermopt as th

# Run configuration: edit these globals, then run this script without arguments.
RUN_FAMILIES = ("cyclopentane", "butane", "novec")
RUN_MODES = ("superheated", "partial_evaporation")
OPTIMIZE = False  # Comparison exports existing solves; no duplicate optimization.
IMPORT_SENSITIVITY_RESULTS = True  # Export the matching 70 C / eta=0.90 solved cases.
COMPARISON_SOURCE_MINIMUM_C = 70
INLET_SUPERHEAT_K = 2.0

FAMILIES = {
    "butane": ("n-Butane", "simple", "Butane / simple"),
    "novec": ("Novec649", "recuperated", "Novec 649 / recuperated"),
    "cyclopentane": ("Cyclopentane", "simple", "Cyclopentane / simple"),
}
RESULTS = HERE / "results" / "fluid_comparison"


class Cycle(th.ThermodynamicCycleOptimization):
    """Use the existing manager with an already loaded configuration dictionary."""

    def read_config(self, config):
        return self.load_config(config)


@lru_cache(maxsize=1)
def fixed_reference_heat_J_per_kg():
    """Source-water enthalpy drop from its inlet to cooling-water inlet T.

    Both endpoint enthalpies use the source inlet pressure (10 bar). This
    reference is independent of the reinjection limit and working fluid.
    """
    config = th.read_configuration_file(HERE / "case_template.yaml")
    fixed = config["problem_formulation"]["fixed_parameters"]
    water = props.Fluid(**fixed["heating_fluid"])
    pressure = fixed["heat_source"]["inlet_pressure"]
    hot = float(water.get_state(props.PT_INPUTS, pressure,
                               fixed["heat_source"]["inlet_temperature"]).h)
    cold = float(water.get_state(props.PT_INPUTS, pressure,
                                fixed["heat_sink"]["inlet_temperature"]).h)
    return hot - cold




def violation(problem):
    return float(max(np.max(np.abs(problem.c_eq), initial=0), np.max(problem.c_ineq, initial=0)))


def physical_cycle_valid(problem):
    energy = problem.cycle_data["energy_analysis"]
    expander = problem.cycle_data["components"]["expander"]
    values = [float(energy[key]) for key in (
        "mass_flow_working_fluid", "mass_flow_heating_fluid", "mass_flow_cooling_fluid",
        "heater_heat_flow", "cooler_heat_flow", "net_system_power")]
    return (all(np.isfinite(value) and value > 0 for value in values)
            and float(expander["state_in"].p) > float(expander["state_out"].p))






def plot_grid(problems, titles, filename):
    """Plot fresh evaluated problems in a 2-by-N grid, save PNG/PDF, and return it."""
    if not problems or len(problems) != len(titles):
        raise ValueError("Provide at least one cycle problem and one title per problem.")

    figure, axes = plt.subplots(2, len(problems), figsize=(5.2 * len(problems), 9.6),
                                squeeze=False)
    for column, (problem, title) in enumerate(zip(problems, titles)):
        # Each problem owns its native artists on its two supplied axes.
        problem.figure, problem.axes = figure, axes[:, column]
        problem.plot_cycle()
        axes[0, column].set_title(title)

    figure.legend(handles=[
        Line2D([], [], color=th.COLORS_MATLAB[index], label=label)
        for index, label in ((1, "Working fluid"), (6, "Heat source"), (0, "Cooling water"))
    ], loc="lower center", ncol=3)
    figure.tight_layout(rect=(0, 0.06, 1, 1))
    th.savefig_in_formats(figure, filename, formats=[".png", ".pdf"])
    return figure


def compare(families):
    summaries = []
    for family in families:
        problems, titles = [], []
        for mode in RUN_MODES:
            directory = RESULTS / f"{family}_{mode}"
            saved = th.read_configuration_file(directory / "optimal_solution.yaml")
            problem = th.ThermodynamicCycleProblem(saved, out_dir=str(directory))
            problem.fitness(problem.x0)
            problems.append(problem)
            summaries.append(pd.read_csv(directory / "summary.csv"))
            inlet_label = f"{INLET_SUPERHEAT_K:g} K superheat" if mode == "superheated" else "Partial evaporation"
            titles.append(f"{FAMILIES[family][2]}\n{inlet_label}")
        figure = plot_grid(problems, titles, RESULTS / f"comparison_grid_{family}")
        plt.close(figure)
    table = pd.concat(summaries, ignore_index=True)
    baseline = table.loc[table["case"].eq("cyclopentane_superheated"), "system_efficiency_pct"]
    if not baseline.empty:
        table["relative_gain_pct"] = 100 * (table["system_efficiency_pct"] / baseline.iloc[0] - 1)
    # Keep system efficiency last in the exported comparison table.
    table = table[[column for column in table if column != "system_efficiency_pct"] + ["system_efficiency_pct"]]
    basename = "comparison" if len(families) == len(FAMILIES) else f"{families[0]}_comparison"
    table.to_csv(RESULTS / f"{basename}.csv", index=False)
    columns = ["case", "inlet_quality", "solver_success", "max_constraint_violation", "function_evaluations"]
    if "relative_gain_pct" in table:
        columns.append("relative_gain_pct")
    columns.append("system_efficiency_pct")
    text = table[columns].to_string(index=False, float_format=lambda x: f"{x:.6g}")
    (RESULTS / f"{basename}.txt").write_text(text + "\n", encoding="utf-8")
    print("\n" + text)


def export_sensitivity_case(family, mode):
    """Export the selected matching sweep solve, preserving its solver records.

    The comparison is a slice of the sensitivity study at 70 C and eta=0.90.
    Reconstructing its native cycle and convergence plots does not optimize it.
    """
    import pysolver_view as psv

    case = f"{mode}_Tmin{COMPARISON_SOURCE_MINIMUM_C:03d}_eta90"
    source = HERE / "results" / f"sensitivity_{family}" / case
    destination = RESULTS / f"{family}_{mode}"
    selected = json.loads((source / "result.json").read_text(encoding="utf-8"))["row"]
    if not selected["success"]:
        raise RuntimeError(f"No accepted sensitivity solution to export: {source}")
    config = th.read_configuration_file(source / "optimal_solution.yaml")
    problem = th.ThermodynamicCycleProblem(config, out_dir=str(destination))
    problem.fitness(problem.x0)
    if not physical_cycle_valid(problem) or violation(problem) > 1e-5:
        raise RuntimeError(f"Saved solution failed reevaluation: {case}")
    efficiency = 100 * float(problem.cycle_data["energy_analysis"]["system_efficiency"])
    if not np.isclose(efficiency, selected["system_efficiency_pct"], rtol=0, atol=1e-7):
        raise RuntimeError(f"Saved efficiency does not match checkpoint: {case}")
    destination.mkdir(parents=True, exist_ok=True)
    for filename in ("optimal_solution.yaml", "convergence_history.txt", "optimization_report.txt"):
        shutil.copy2(source / filename, destination / filename)
    shutil.copy2(source / "result.json", destination / "source_result.json")
    shutil.copy2(source / "attempts.csv", destination / "source_attempts.csv")

    bounds = pd.DataFrame({"variable": problem.variable_names, "lower": problem.lb,
                           "value": problem.x0, "upper": problem.ub})
    fraction = (bounds["value"] - bounds["lower"]) / (bounds["upper"] - bounds["lower"])
    bounds["distance_to_nearest_bound_fraction"] = np.minimum(fraction, 1 - fraction)
    bounds["active_bound"] = np.where(fraction < 1e-4, "lower", np.where(fraction > 1 - 1e-4, "upper", ""))
    bounds.to_csv(destination / "design_variable_bounds.csv", index=False)
    state = bounds[bounds["variable"].str.startswith(("compressor_", "expander_"))]
    summary = dict(selected, case=f"{family}_{mode}", source_case=case,
                   source_study=f"sensitivity_{family}",
                   minimum_state_bound_distance=float(state["distance_to_nearest_bound_fraction"].min()),
                   active_bounds="; ".join(bounds.loc[bounds["active_bound"] != "", "variable"]))
    pd.DataFrame([summary]).to_csv(destination / "summary.csv", index=False)
    problem.save_data_to_excel(filename=str(destination / "optimal_solution.xlsx"))
    problem.plot_cycle()
    problem.figure.tight_layout(pad=1)
    problem.figure.suptitle(None)
    th.savefig_in_formats(problem.figure, destination / "optimal_solution", formats=[".png", ".pdf"])
    plt.close(problem.figure)

    # Plot the recorded history with the native solver plotting function. No
    # counters, KKT values, or optimization statistics are recalculated here.
    columns = ("grad_count", "func_count", "objective_value", "constraint_violation", "norm_step")
    history = {name: [] for name in columns}
    for line in (source / "convergence_history.txt").read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) != 6:
            continue
        try:
            values = [float(value) for value in fields]
        except ValueError:
            continue
        if values[0].is_integer() and values[1].is_integer():
            for name, value in zip(columns, values[:5]):
                history[name].append(int(value) if name.endswith("count") else value)
    if not history["objective_value"]:
        raise ValueError(f"Missing native solver history in {source}")
    solver = psv.OptimizationSolver(problem, update_on="function", print_convergence=False,
                                    plot_convergence=False, plot_scale_constraints="log")
    solver.x_final = problem.x0.copy()
    solver.success = bool(selected["solver_success"])
    solver.message = selected["solver_message"]
    solver.elapsed_time = selected["solver_seconds"]
    solver.func_count_tot = selected["function_evaluations"]
    solver.grad_count = selected["gradient_evaluations"]
    solver.convergence_history.update(history)
    solver.plot_convergence_history(savefile=True, filename="convergence_history",
                                    output_dir=str(destination), showfig=False)
    plt.close("all")
    print(f"Exported {family}/{case}: system efficiency {efficiency:.6f}%", flush=True)


def main():
    if OPTIMIZE:
        raise RuntimeError("Run the sensitivity script for the single optimization of each case.")
    if IMPORT_SENSITIVITY_RESULTS:
        for family in RUN_FAMILIES:
            for mode in RUN_MODES:
                export_sensitivity_case(family, mode)
    compare(RUN_FAMILIES)


if __name__ == "__main__":
    main()
