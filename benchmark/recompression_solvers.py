"""Compare three solvers with CoolProp and bicubic on recompression CO2.

Run from any directory after installing ThermOpt. The working-fluid table is
128 by 128 for bicubic, independently of the YAML's normal table resolution.
The input YAML is read only. Each solver receives the same physical initial
design, bounds, constraints, finite differences and problem scaling.
"""

import argparse
import copy
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import time

import jax
import matplotlib.pyplot as plt
import numpy as np
import pysolver_view as psv
import thermopt as th

EXAMPLE_DIR = Path(__file__).resolve().parents[1] / "case_studies/sCO2_recompression"
GRID_SIZE = 128
FEASIBILITY_TOL = 1e-6
OPTIMALITY_TOL = 1e-4
BACKENDS = ("coolprop", "bicubic")
SOLVERS = (
    ("slsqp", "Controlled SLSQP", "pysolver", "slsqp"),
    ("pounce_ipm", "Pounce interior point", "pounce", "ipm"),
    ("pounce_sqp", "Pounce SQP", "pounce", "sqp"),
)


def comparison_configuration(config, table_directory, backend="bicubic"):
    """Select a native working-fluid evaluator without modifying the YAML.

    Auxiliary fluids keep their configured CoolProp models in both runs.
    Removing bicubic options is required for an ordinary CoolProp constructor.
    """
    if backend not in BACKENDS:
        raise ValueError(f"Property backend must be one of {BACKENDS}.")
    configured = copy.deepcopy(config)
    fluid = configured["problem_formulation"]["fixed_parameters"]["working_fluid"]
    fluid["method"] = backend
    if backend == "bicubic":
        fluid["table_options"].update(
            N_h=GRID_SIZE, N_p=GRID_SIZE, table_dir=str(table_directory.resolve())
        )
    else:
        fluid.pop("table_options", None)
    return configured


def solver_options(config, library, method, max_iterations, verbose):
    """Keep common numerics while passing each backend only its own controls."""
    options = copy.deepcopy(config["solver_options"])
    options.pop("callbacks", None)
    options.update(
        library=library,
        method=method,
        max_iterations=int(max_iterations or options["max_iterations"]),
        update_on="iteration",
        print_convergence=verbose,
        plot_convergence=False,
    )
    common = dict(feasibility_tol=FEASIBILITY_TOL, optimality_tol=OPTIMALITY_TOL)
    if library == "pysolver":
        common.update(
            first_step_fraction=float(options["extra_options"]["first_step_fraction"]),
            max_step_fraction=1.0,
            step_growth_factor=1.5,
        )
    else:
        common["stop_on_kkt"] = True
        if method == "sqp":
            # The default 1e-9 inner-QP tolerances stalled at the first QP for
            # this discretized HX problem. These remain stricter than the
            # outer feasibility threshold and work with numerical derivatives.
            common.update(sqp_qp_feas_tol=1e-7, sqp_qp_opt_tol=1e-7)
    # Pounce has its own line search; it does not expose SLSQP step fractions.
    options["extra_options"] = common
    return options


def plain(value):
    """Serialize native numerical storage, retaining unavailable data as null."""
    if isinstance(value, dict):
        return {
            key.item() if isinstance(key, np.generic) else key: plain(item)
            for key, item in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [plain(item) for item in value]
    if hasattr(value, "__array__"):
        return plain(np.asarray(value).tolist())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def save_json(path, value):
    """Save actual observations and settings without nonstandard JSON NaNs."""
    path.write_text(
        json.dumps(plain(value), indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


def assess_optimality(history, config, directory):
    """Use the same independent KKT assessment at each recorded solver point.

    These postprocessing probes do not enter solve time, model-call totals or
    trajectories. Pounce native observations are not certified accepted steps.
    Restoration observations without a user-space point remain in raw events.
    """
    problem = th.ThermodynamicCycleProblem(
        copy.deepcopy(config["problem_formulation"]), out_dir=str(directory)
    )
    problem.fitness(problem.x0)
    diagnostics = psv.OptimizationSolver(
        problem,
        library="pysolver",
        method="slsqp",
        problem_scale=config["solver_options"]["problem_scale"],
        derivative_method=config["solver_options"]["derivative_method"],
        derivative_abs_step=config["solver_options"]["derivative_abs_step"],
        print_convergence=False,
    )
    values, failures = [], []
    for index, point in enumerate(history["x"]):
        try:
            assessment = diagnostics.evaluate_kkt_conditions(point, tol=FEASIBILITY_TOL)
            values.append(float(assessment["first_order_optimality"]))
        except psv.DomainError as error:
            values.append(float("nan"))
            failures.append(dict(observation=index, message=str(error)))
    history["comparison_optimality"] = values
    return 1 + diagnostics.func_count_tot, failures


def run_case(case, config, output, max_iterations=None, verbose=False):
    """Solve from the shared initial design and preserve honest stopped results."""
    name, label, library, method = case
    backend = config["problem_formulation"]["fixed_parameters"]["working_fluid"][
        "method"
    ]
    directory = output / backend / name
    directory.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    problem = th.ThermodynamicCycleProblem(
        copy.deepcopy(config["problem_formulation"]), out_dir=str(directory)
    )
    initial = problem.x0.copy()
    jax.block_until_ready(problem.fitness(initial))
    options = solver_options(config, library, method, max_iterations, verbose)
    raw_events = []
    solver = psv.OptimizationSolver(
        problem, **options, backend_callback=raw_events.append
    )
    setup_seconds = time.perf_counter() - started
    properties = (
        f"bicubic {GRID_SIZE}x{GRID_SIZE}" if backend == "bicubic" else "CoolProp HEOS"
    )
    print(f"\nRunning {label}, {properties}", flush=True)
    started = time.perf_counter()
    error = None
    try:
        solver.solve(initial.copy())
    except psv.DomainError as exception:
        error = str(exception)
    # Include pending work in the complete returned cycle before stopping time.
    jax.block_until_ready(problem.cycle_data)
    solve_seconds = time.perf_counter() - started
    history = copy.deepcopy(solver.convergence_history)
    diagnostic_calls, diagnostic_failures = assess_optimality(
        history, config, directory
    )
    final_point = getattr(solver, "x_final", None)
    final_fitness = None
    if final_point is not None:
        try:
            final_fitness = problem.fitness(np.asarray(final_point))
        except psv.DomainError as exception:
            error = str(exception)
    summary = dict(
        name=name,
        label=label,
        property_backend=backend,
        grid=[GRID_SIZE, GRID_SIZE] if backend == "bicubic" else None,
        success=bool(solver.success) and error is None,
        message=error or str(solver.message),
        options=options,
        initial_design=initial,
        final_design=final_point,
        final_fitness=final_fitness,
        setup_seconds=setup_seconds,
        solve_seconds=solve_seconds,
        iterations=solver.iteration_count,
        model_calls=solver.func_count_tot,
        postprocessing_model_calls=diagnostic_calls,
        postprocessing_failures=diagnostic_failures,
        final_kkt=(
            getattr(solver, "result", {}).get("kkt")
            if solver.result is not None
            else None
        ),
        backend_result=solver.backend_result,
    )
    save_json(directory / "summary.json", summary)
    save_json(directory / "history.json", history)
    save_json(directory / "backend_events.json", raw_events)
    fields = (
        "iteration",
        "event_kind",
        "func_count_total",
        "objective_value",
        "constraint_violation",
        "comparison_optimality",
    )
    with (directory / "history.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        writer.writerows(zip(*(history[field] for field in fields)))
    print(
        f"{'Converged' if summary['success'] else 'Stopped'}: {summary['message']} ({solve_seconds:.2f} s)",
        flush=True,
    )
    return summary, history


def plot_comparison(results, output):
    """Use solver colors and backend line styles for the six actual histories."""
    th.set_plot_options(fontsize=11, grid=False)
    figure, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    colors = dict(
        zip((case[0] for case in SOLVERS), ("tab:blue", "tab:orange", "tab:green"))
    )
    for summary, history in results:
        backend = summary["property_backend"]
        label = f"{summary['label']} / {backend}" + (
            "" if summary["success"] else " (stopped)"
        )
        for ax, values in zip(
            axes,
            (
                -np.asarray(history["objective_value"]),
                history["constraint_violation"],
                history["comparison_optimality"],
            ),
        ):
            ax.plot(
                history["iteration"],
                values,
                label=label,
                color=colors[summary["name"]],
                linestyle="--" if backend == "coolprop" else "-",
                marker=".",
                markersize=4,
            )
    axes[0].set(title="Cycle efficiency", ylabel="Efficiency [-]")
    axes[1].set(title="Constraint violation", ylabel="Maximum normalized violation [-]")
    axes[2].set(title="First-order optimality", ylabel=r"$\|\nabla_x L\|_\infty$")
    for ax, tolerance in zip(axes[1:], (FEASIBILITY_TOL, OPTIMALITY_TOL)):
        ax.set_yscale("symlog", linthresh=1e-8)
        ax.axhline(tolerance, color="0.5", linestyle=":", linewidth=1)
        ax.set_ylim(bottom=0)
    for ax in axes:
        ax.set_xlabel("Solver iteration / native observation")
        ax.set_xlim(left=0)
        ax.grid(axis="y", alpha=0.2)
    figure.legend(
        *axes[0].get_legend_handles_labels(),
        loc="lower center",
        ncol=3,
        frameon=False,
        fontsize=9,
    )
    figure.tight_layout(rect=(0, 0.18, 1, 1))
    for suffix in ("png", "svg", "pdf"):
        figure.savefig(
            output / f"solver_comparison.{suffix}", dpi=200, bbox_inches="tight"
        )
    return figure


def timing_table(results):
    """Format one observed complete solve per backend/solver, including status.

    Setup includes native model construction/table loading and initial cycle
    evaluation. Solve time includes the backend's own stopping/final KKT work;
    extra history diagnostics, output files and plots are excluded.
    """
    headers = (
        "Solver",
        "Properties",
        "Setup [s]",
        "Solve [s]",
        "Iter",
        "Calls",
        "Eff. [%]",
        "Status",
    )
    rows = []
    for summary, _ in results:
        fitness = summary["final_fitness"]
        efficiency = "n/a" if fitness is None else f"{-float(fitness[0])*100:.5f}"
        rows.append(
            (
                summary["label"],
                (
                    "Bicubic 128x128"
                    if summary["property_backend"] == "bicubic"
                    else "CoolProp HEOS"
                ),
                f"{summary['setup_seconds']:.3f}",
                f"{summary['solve_seconds']:.3f}",
                str(summary["iterations"]),
                str(summary["model_calls"]),
                efficiency,
                "Converged" if summary["success"] else "Stopped",
            )
        )
    widths = [
        max(len(row[index]) for row in (headers, *rows))
        for index in range(len(headers))
    ]
    border = "+" + "+".join("-" * (width + 2) for width in widths) + "+"

    def row_text(row):
        return (
            "| "
            + " | ".join(
                value.ljust(width) if index in (0, 1, 7) else value.rjust(width)
                for index, (value, width) in enumerate(zip(row, widths))
            )
            + " |"
        )

    lines = [
        "One complete solve per solver/backend (seconds)",
        "Setup includes table loading/generation and initial-cycle warmup.",
        "Solve excludes setup, extra history diagnostics and plots.",
        border,
        row_text(headers),
        border,
    ]
    lines.extend(row_text(row) for row in rows)
    lines.append(border)
    return "\n".join(lines)


def main():
    """Run the six backend/solver combinations and print their observed times."""
    # Thermodynamic finite differences need full precision in native JAX tables.
    jax.config.update("jax_enable_x64", True)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=EXAMPLE_DIR / "case_sCO2_recompression.yaml"
    )
    parser.add_argument(
        "--output", type=Path, default=EXAMPLE_DIR / "results" / "solver_comparison"
    )
    parser.add_argument("--max-iterations", type=int)
    parser.add_argument("--no-show", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    if args.max_iterations is not None and args.max_iterations < 1:
        parser.error("--max-iterations must be positive")
    source = args.config.resolve()
    before = hashlib.sha256(source.read_bytes()).hexdigest()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    base_config = th.read_configuration_file(str(source))
    configurations = {
        backend: comparison_configuration(base_config, output / "tables", backend)
        for backend in BACKENDS
    }
    for backend, config in configurations.items():
        (output / backend).mkdir(exist_ok=True)
        save_json(output / backend / "configuration.json", config)
    results = [
        run_case(
            case, configurations[backend], output, args.max_iterations, args.verbose
        )
        for case in SOLVERS
        for backend in BACKENDS
    ]
    assert all(
        np.array_equal(results[0][0]["initial_design"], item[0]["initial_design"])
        for item in results
    )
    unchanged = hashlib.sha256(source.read_bytes()).hexdigest() == before
    if not unchanged:
        raise RuntimeError(
            "The input YAML changed during the comparison; repeat with a stable configuration."
        )
    save_json(
        output / "comparison.json",
        dict(
            yaml=str(source),
            yaml_sha256=before,
            yaml_unchanged=unchanged,
            grid=[GRID_SIZE, GRID_SIZE],
            property_backends=list(BACKENDS),
            timing_scope="One sequential solve per combination. Setup includes model/table construction and initial complete-cycle warmup. Solve includes native/frontend stopping and final KKT work, excludes extra history diagnostics, output and plots.",
            versions={
                name: importlib.metadata.version(name)
                for name in ("thermopt", "jaxprop", "pysolver_view", "pounce-solver")
            },
            observation_semantics="SLSQP: accepted steps; Pounce: native observations using the last objective point, not guaranteed accepted steps. Raw restoration events are retained separately.",
            results=[item[0] for item in results],
        ),
    )
    table = timing_table(results)
    (output / "timing_comparison.txt").write_text(table + "\n", encoding="utf-8")
    plot_comparison(results, output)
    print(f"\nComparison saved to {output}", flush=True)
    print(f"\n{table}", flush=True)
    if args.no_show:
        plt.close("all")
    else:
        plt.show()


if __name__ == "__main__":
    main()
