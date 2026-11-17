"""Rerun the recuperated sCO2 example and export backend comparisons.

Run with the repository's virtual environment from any working directory.
Generated configurations, histories, and summaries stay beside this script.
"""

import argparse
import copy
import csv
import json
import os
from contextlib import nullcontext
from importlib.metadata import version
from pathlib import Path
from time import perf_counter

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import thermopt as th


WORK_DIR = Path(__file__).resolve().parent
SOURCE_CONFIG = WORK_DIR.parents[1] / "examples/sCO2_recuperated/case_sCO2_recuperated.yaml"


def json_value(value):
    """Convert NumPy arrays and scalars for the result files."""
    if hasattr(value, "__array__"):
        return np.asarray(value).tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, default=json_value), encoding="utf-8")


def cycle_snapshot(problem, point):
    """Reevaluate a physical design and record feasibility and cycle performance."""
    problem.fitness(np.asarray(point, dtype=float))
    data = problem.cycle_data
    states = {}
    for component in ("compressor", "expander"):
        states[component] = {
            end: {
                prop: float(getattr(data["components"][component][end], prop))
                for prop in ("p", "T", "h", "s", "rho")
            }
            for end in ("state_in", "state_out")
        }
    lower, upper = map(np.asarray, problem.get_bounds())
    bound_violation = np.maximum(lower - point, point - upper) / (upper - lower)
    checked_values = np.r_[
        problem.f,
        problem.c_eq,
        problem.c_ineq,
        bound_violation,
        list(data["energy_analysis"].values()),
        [value for component in states.values() for state in component.values() for value in state.values()],
    ]
    if not np.all(np.isfinite(checked_values)):
        raise ValueError("Cycle evaluation contains nonfinite performance or constraint values.")
    return {
        "design_variables": dict(zip(problem.variable_names, point)),
        "energy_analysis": copy.deepcopy(data["energy_analysis"]),
        "states": states,
        "minimum_temperature_difference_K": {
            name: float(np.min(data["components"][name]["temperature_difference"]))
            for name in ("heater", "recuperator", "cooler")
        },
        "maximum_normalized_constraint_violation": float(
            max(0.0, np.max(problem.c_ineq), np.max(bound_violation))
        ),
    }


def run_case(backend, output_dir, grid_size=100, all_fluids=False):
    """Use the example's normal optimization manager with plotting disabled."""
    run_name = backend if backend == "coolprop" else f"bicubic_{grid_size}"
    if backend == "bicubic" and all_fluids:
        run_name = f"bicubic_all_{grid_size}"
    run_dir = output_dir / run_name
    run_dir.mkdir(parents=True, exist_ok=True)
    config = th.read_configuration_file(str(SOURCE_CONFIG))
    config["solver_options"]["callbacks"] = {}
    config["solver_options"]["plot_convergence"] = False
    config["solver_options"]["max_iterations"] = int(
        config["solver_options"]["max_iterations"]
    )
    config_path = run_dir / "config.yaml"
    import yaml

    config_path.write_text(
        yaml.safe_dump(json.loads(json.dumps(config, default=json_value)), sort_keys=False),
        encoding="utf-8",
    )
    context = nullcontext({})
    if backend == "bicubic":
        from bicubic_adapter import bicubic_cycle

        context = bicubic_cycle(grid_size, output_dir / "tables", all_fluids=all_fluids)
    setup_start = perf_counter()
    with context as audit:
        table_seconds = perf_counter() - setup_start
        start = perf_counter()
        cycle = th.ThermodynamicCycleOptimization(str(config_path), out_dir=str(run_dir))
        initial = cycle_snapshot(cycle.problem, cycle.problem.x0.copy())
        cycle.run_optimization()
        solve_seconds = perf_counter() - start
        final_point = cycle.solver.x_final.copy()
        final = cycle_snapshot(cycle.problem, final_point)
        cross_evaluations = {}
        reference_path = output_dir / "coolprop/summary.json"
        if backend == "bicubic" and reference_path.exists():
            reference = json.loads(reference_path.read_text(encoding="utf-8"))
            reference_point = np.array(
                [reference["final"]["design_variables"][name] for name in cycle.problem.variable_names]
            )
            cross_evaluations["coolprop_optimum_with_bicubic"] = cycle_snapshot(
                cycle.problem, reference_point
            )
    if backend == "bicubic":
        reference_problem = th.ThermodynamicCycleProblem(
            copy.deepcopy(config["problem_formulation"]), out_dir=str(run_dir)
        )
        cross_evaluations["bicubic_optimum_with_coolprop"] = cycle_snapshot(
            reference_problem, final_point
        )
    solver = cycle.solver
    history = solver.convergence_history
    write_json(run_dir / "history.json", history)
    fields = ("iteration", "event_kind", "func_count_total", "objective_value", "constraint_violation")
    with (run_dir / "history.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        writer.writerows(zip(*(history[field] for field in fields)))
    summary = {
        "backend": backend,
        "source_config": str(SOURCE_CONFIG),
        "versions": {
            name: version(name)
            for name in ("thermopt", "jaxprop", "pysolver_view", "CoolProp", "scipy", "numpy")
        },
        "success": bool(solver.success),
        "message": str(solver.message),
        "iterations": int(solver.iteration_count),
        "model_calls": int(solver.func_count_tot),
        "final_optimality_solver_coordinates": float(history["optimality"][-1]),
        "setup_and_solve_seconds": solve_seconds,
        "table_setup_seconds": table_seconds,
        "bicubic_audit": audit,
        "initial": initial,
        "final": final,
        "cross_evaluations": cross_evaluations,
    }
    write_json(run_dir / "summary.json", summary)
    print(
        f"{run_name}: success={summary['success']}, "
        f"efficiency={100 * final['energy_analysis']['cycle_efficiency']:.9f}%, "
        f"iterations={summary['iterations']}, solve={solve_seconds:.3f}s",
        flush=True,
    )
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=("coolprop", "bicubic"), default="coolprop")
    parser.add_argument("--grid-size", type=int, default=100)
    parser.add_argument("--all-fluids", action="store_true")
    parser.add_argument("--output", type=Path, default=WORK_DIR / "results")
    args = parser.parse_args()
    run_case(args.backend, args.output.resolve(), args.grid_size, args.all_fluids)


if __name__ == "__main__":
    main()
