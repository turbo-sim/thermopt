"""Optimize two simple Rankine cycles and compare uniform and snapped HX grids."""

from contextlib import redirect_stderr, redirect_stdout
from copy import deepcopy
from importlib.metadata import version
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
import thermopt as th

# Each grid starts from the same YAML guess for its inlet formulation.
UNIFORM_NODES = (10, 50, 100, 200)
SNAP_NODES = 10
CHECK_NODES = 1000
CASES = {"superheated": "20 K superheat", "partial_evaporation": "Partial evaporation"}
RESULTS = HERE / "results"
FIGURES = HERE / "figures"


def profiles(problem):
    """Both streams share the cold-side heat-duty coordinate, including snapped nodes."""
    data = {}
    for name in ("heater", "cooler"):
        hx = problem.cycle_data["components"][name]
        hot, cold = hx["hot_side"], hx["cold_side"]
        x = (cold["states"].h - cold["state_in"].h) / (
            cold["state_out"].h - cold["state_in"].h)
        data[name] = pd.DataFrame({
            "duty_fraction": x,
            "hot_temperature_C": hot["states"].T - 273.15,
            "cold_temperature_C": cold["states"].T - 273.15,
            "approach_K": hx["temperature_difference"],
        })
    return data


def constraint_violation(problem):
    return float(max(np.max(np.abs(problem.c_eq), initial=0),
                     np.max(problem.c_ineq, initial=0)))


def solve_case(base, mode, nodes, snap):
    config = deepcopy(base)
    inlet = config.pop("inlet_cases")[mode]
    formulation = config["problem_formulation"]
    for name, value in inlet["initial_guess"].items():
        formulation["design_variables"][name]["value"] = value
    formulation["constraints"].extend(inlet["constraints"])
    for name in ("heater", "cooler"):
        formulation["fixed_parameters"][name].update(
            num_elements=nodes, include_saturation_nodes=snap)

    name = f"{mode}_{'snap' if snap else 'uniform'}_{nodes:03d}"
    folder = RESULTS / name
    folder.mkdir(parents=True, exist_ok=True)
    input_file = folder / "input.yaml"
    input_file.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    with (folder / "solver.log").open("w", encoding="utf-8") as log:
        with redirect_stdout(log), redirect_stderr(log):
            cycle = th.ThermodynamicCycleOptimization(str(input_file), out_dir=str(folder))
            cycle.run_optimization()
            problem, solver = cycle.problem, cycle.solver
            function_evaluations = int(solver.func_count_tot)
            problem.fitness(solver.x_final)
            problem.save_current_configuration(str(folder / "solution.yaml"))
            cycle.print_convergence_history(savefile=True)
            cycle.print_optimization_report(savefile=True)
            # Reporting may evaluate perturbed states; restore the actual solution.
            problem.fitness(solver.x_final)

    energy = problem.cycle_data["energy_analysis"]
    turbine_inlet = problem.cycle_data["components"]["expander"]["state_in"]
    valid = all(np.isfinite(float(energy[key])) and float(energy[key]) > 0 for key in (
        "mass_flow_working_fluid", "mass_flow_heating_fluid", "mass_flow_cooling_fluid",
        "heater_heat_flow", "cooler_heat_flow", "net_system_power"))
    valid = valid and abs(float(energy["energy_balance"])) < 0.1
    violation = constraint_violation(problem)
    row = dict(case=mode, nodes=nodes, snapping=snap,
               success=bool(solver.success and violation < 1e-5 and valid),
               solver_success=bool(solver.success), solver_message=str(solver.message),
               system_efficiency_pct=100 * float(energy["system_efficiency"]),
               turbine_inlet_superheat_K=(0.0 if turbine_inlet.is_two_phase
                                         else float(turbine_inlet.superheating)),
               turbine_inlet_quality=float(turbine_inlet.Q),
               turbine_inlet_is_two_phase=bool(turbine_inlet.is_two_phase),
               constraint_violation=violation,
               energy_balance_W=float(energy["energy_balance"]),
               function_evaluations=function_evaluations,
               solve_seconds=float(solver.elapsed_time))
    coarse = profiles(problem)

    # Reevaluate this same design, without optimizing, on the common check grid.
    dense_config = th.read_configuration_file(folder / "solution.yaml")
    for name in ("heater", "cooler"):
        dense_config["fixed_parameters"][name].update(
            num_elements=CHECK_NODES, include_saturation_nodes=False)
    dense_problem = th.ThermodynamicCycleProblem(dense_config, out_dir=str(folder))
    dense_problem.fitness(dense_problem.x0)
    dense = profiles(dense_problem)
    row["dense_constraint_violation"] = constraint_violation(dense_problem)
    for name in ("heater", "cooler"):
        coarse[name].to_csv(folder / f"{name}_model.csv", index=False)
        dense[name].to_csv(folder / f"{name}_check.csv", index=False)
        row[f"{name}_model_min_K"] = float(coarse[name].approach_K.min())
        row[f"{name}_check_min_K"] = float(dense[name].approach_K.min())
    (folder / "summary.json").write_text(json.dumps(row, indent=2), encoding="utf-8")
    plt.close("all")
    return row, coarse, dense


def save_figure(fig, name):
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES / f"{name}.{suffix}", dpi=180, bbox_inches="tight")
    plt.close(fig)


def make_figures(table, snapped_profiles):
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 10,
                         "axes.labelsize": 10, "axes.titlesize": 12,
                         "xtick.labelsize": 9, "ytick.labelsize": 9,
                         "legend.fontsize": 9, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex="col")
    for col, (mode, title) in enumerate(CASES.items()):
        rows = table.loc[table.case.eq(mode)]
        uniform = rows.loc[~rows.snapping].sort_values("nodes")
        snap = rows.loc[rows.snapping].iloc[0]
        axes[0, col].plot(uniform.nodes, uniform.system_efficiency_pct, "o-",
                          color="#4477AA", label="Uniform")
        axes[0, col].plot(snap.nodes, snap.system_efficiency_pct, "D",
                          color="#EE7733", ms=7, label=f"{SNAP_NODES} + snapping")
        axes[0, col].set(title=title, ylabel="System efficiency [%]")
        axes[0, col].legend()
        for hx, label, color in (("heater", "Evaporator", "#4477AA"),
                                  ("cooler", "Condenser", "#228833")):
            axes[1, col].plot(uniform.nodes, uniform[f"{hx}_check_min_K"], "o-",
                              color=color, label=f"{label}: uniform")
            axes[1, col].plot(snap.nodes, snap[f"{hx}_check_min_K"], "D", ms=7,
                              color=color, markerfacecolor=color,
                              markeredgecolor="#EE7733", markeredgewidth=2,
                              label=f"{label}: snapping")
        axes[1, col].axhline(5, color="black", ls="--", lw=1, label="5 K requirement")
        axes[1, col].set(xlabel="Nodes per heat exchanger",
                          ylabel=f"Minimum approach on {CHECK_NODES}-node check [K]")
        axes[1, col].legend(fontsize=8)
        for ax in axes[:, col]:
            ax.set_xscale("log")
            ax.set_xticks(UNIFORM_NODES, labels=[str(n) for n in UNIFORM_NODES])
            ax.grid(alpha=.2)
    save_figure(fig, "grid_sensitivity")

    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    for col, (mode, title) in enumerate(CASES.items()):
        coarse, dense = snapped_profiles[mode]
        for r, name in enumerate(("heater", "cooler")):
            ax = axes[r, col]
            for side, color in (("hot", "#CC3311"), ("cold", "#0077BB")):
                column = f"{side}_temperature_C"
                ax.plot(100 * dense[name].duty_fraction, dense[name][column], color=color,
                        label=f"{side.capitalize()}: {CHECK_NODES}-node check")
                ax.plot(100 * coarse[name].duty_fraction, coarse[name][column], "o",
                        ms=5, fillstyle="none", color=color, label=f"{side.capitalize()}: snapped nodes")
            ax.set(title=f"{title} / {'evaporator' if name == 'heater' else 'condenser'}",
                   xlabel="Heat duty from cold inlet [%]", ylabel="Temperature [°C]")
            ax.grid(alpha=.2)
            ax.legend(fontsize=8)
    save_figure(fig, "temperature_profiles")


def write_readme(table):
    lines = [
        "# Heat-exchanger node snapping\n",
        "Optimize a simple cyclopentane Rankine cycle with either **20 K turbine-inlet superheat** "
        "or **partial evaporation**. Compare uniform heat-exchanger grids with a coarse grid that "
        "snaps nodes to the liquid/vapor saturation boundaries.\n",
        "## Run\n",
        "From the repository root, with ThermoOpt and its dependencies installed:\n",
        "```console\npython examples/demo_hx_node_snap/run_demo.py\n```\n",
        "`case.yaml` contains the cycle, solver settings, initial guesses and inlet constraints. "
        "Edit `UNIFORM_NODES`, `SNAP_NODES` or `CHECK_NODES` at the top of `run_demo.py` to change "
        "the comparison. Each run performs all optimizations again and updates this README and the figures.\n",
        "Enable snapping on both heat exchangers using:\n",
        "```yaml\nheater:\n  num_elements: 10\n  include_saturation_nodes: true\n"
        "cooler:\n  num_elements: 10\n  include_saturation_nodes: true\n```\n",
        "These entries belong under `problem_formulation.fixed_parameters`. Set "
        "`include_saturation_nodes: false` for the unchanged uniform calculation. "
        "Snapping moves the interior two-phase nodes next to saturation boundaries and adjusts "
        "the opposite stream at the same heat duty. It preserves endpoints and node counts, "
        "and uses a linear pressure estimate without an iterative root finder.\n",
        "## Results\n",
        "The heat source is water at 140 °C and 10 bar, with a 70 °C minimum outlet; cooling "
        "water enters at 20 °C and leaves at approximately 25 °C. Turbine efficiency is 90%, "
        "pump efficiencies are 70%, pressure losses are 1% per HX stream, and minimum approach "
        "is 5 K. Net cycle output is 1 MW before the external-water pumps. System efficiency "
        "includes those pumps and references source heat available down to 70 °C.\n",
        f"All {len(table)} cases start from the same YAML guess within each inlet formulation "
        "and use the same SLSQP settings. Each final design is independently reevaluated with "
        f"**{CHECK_NODES} uniform nodes and snapping disabled**, without reoptimization. "
        f"{int(table.success.sum())}/{len(table)} optimizations converged and satisfied their sampled constraints.\n",
        f"| Inlet case | HX grid | System efficiency [%] | Evaporator check [K] | Condenser check [K] | Solve time [s] |\n"
        "|---|---|---:|---:|---:|---:|",
    ]
    for r in table.itertuples():
        grid = f"{r.nodes} + snap" if r.snapping else f"{r.nodes} uniform"
        lines.append(f"| {CASES[r.case]} | {grid} | {r.system_efficiency_pct:.4f} | "
                     f"{r.heater_check_min_K:.4f} | {r.cooler_check_min_K:.4f} | {r.solve_seconds:.2f} |")
    lines.append("\n![Grid sensitivity](figures/grid_sensitivity.png)\n")
    for mode, label in CASES.items():
        rows = table.loc[table.case.eq(mode)]
        snap = rows.loc[rows.snapping].iloc[0]
        fine = rows.loc[~rows.snapping].sort_values("nodes").iloc[-1]
        quality = (f" Its optimized turbine-inlet vapor fraction is {100*snap.turbine_inlet_quality:.1f}%."
                   if mode == "partial_evaporation" else "")
        lines.append(f"- **{label}:** {SNAP_NODES} snapped nodes give {snap.system_efficiency_pct:.4f}% efficiency "
                     f"versus {fine.system_efficiency_pct:.4f}% with {fine.nodes} uniform nodes; "
                     f"the snapped design's evaporator/condenser checks are "
                     f"{snap.heater_check_min_K:.4f}/{snap.cooler_check_min_K:.4f} K.{quality}")
    lines.extend([
        "",
        "A coarse uniform grid can miss a sharp pinch and therefore report optimistic efficiency. "
        "Judge efficiency together with the independent approach-temperature checks, not by objective "
        "agreement alone. Snapping places the few available nodes at the phase boundaries where "
        "these pinches occur.\n",
        "![Snapped nodes and independent temperature profiles](figures/temperature_profiles.png)\n",
        "The circles are the snapped model nodes; the lines reevaluate the same optimized designs "
        "on the uniform check grid. A finite check grid can still miss an exact corner. Pressure "
        "location is approximate with pressure loss, and a grid must contain interior two-phase "
        "nodes to snap them. These are single-start local optimizations; solve times exclude checks "
        "and plotting and are not a controlled speed benchmark.\n",
        "## Files\n",
        "- `case.yaml`: all cycle and solver inputs, including both inlet formulations.\n"
        "- `run_demo.py`: optimization, node sensitivity, independent checks and report generation.\n"
        "- `results/`: per-case input/solution YAML, solver records, model/check profiles, and `summary.csv` "
        "(ignored by Git).\n"
        "- `figures/`: the PNG figures above and matching PDFs, retained in Git.\n",
    ])
    (HERE / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    base = yaml.safe_load((HERE / "case.yaml").read_text(encoding="utf-8"))
    rows, snapped_profiles = [], {}
    grids = [(n, False) for n in UNIFORM_NODES] + [(SNAP_NODES, True)]
    for mode in CASES:
        for nodes, snap in grids:
            print(f"{CASES[mode]}: {nodes} nodes, snapping={snap}", flush=True)
            row, coarse, dense = solve_case(base, mode, nodes, snap)
            rows.append(row)
            if snap:
                snapped_profiles[mode] = coarse, dense
            pd.DataFrame(rows).to_csv(RESULTS / "summary.csv", index=False)
            print(f"  accepted={row['success']}, efficiency={row['system_efficiency_pct']:.4f}%, "
                  f"check={row['heater_check_min_K']:.4f}/{row['cooler_check_min_K']:.4f} K", flush=True)
    table = pd.DataFrame(rows)
    make_figures(table, snapped_profiles)
    write_readme(table)
    versions = {name: version(name) for name in ("numpy", "scipy", "jaxprop", "CoolProp")}
    (RESULTS / "environment.json").write_text(json.dumps(versions, indent=2), encoding="utf-8")
    if not table.success.all():
        raise RuntimeError("Some optimizations failed; see results/summary.csv and per-case solver logs.")
    print(f"Finished. Results and instructions: {HERE / 'README.md'}", flush=True)


if __name__ == "__main__":
    main()
