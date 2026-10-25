"""Compare native and controlled SLSQP on the same recuperated CO2 cycle.

Run from the repository root with ``poetry run python
examples/sCO2_recuperated/compare_slsqp.py``. Add ``--no-show`` for a headless
run. The YAML configuration and physical initial guess are shared by all runs.
"""

import argparse
import copy
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pysolver_view as psv
import thermopt as th

EXAMPLE_DIR = Path(__file__).resolve().parent
FIRST_STEP_FRACTIONS = (1.0, 0.3, 0.2, 0.1)
FEASIBILITY_TOL = 1e-6
OPTIMALITY_TOL = 1e-4


class ComparisonCycleProblem(th.ThermodynamicCycleProblem):
    """Expose property-domain failures to the custom solver's backtracking."""

    def fitness(self, x):
        """Evaluate the cycle, translating JaxProp's property-failure wrapper."""
        try:
            return super().fitness(x)
        except RuntimeError as error:
            if not str(error).startswith("Thermodynamic property calculations failed:"):
                raise
            raise psv.DomainError(str(error)) from error


def comparison_cases(current_scale):
    """Return the two native runs and four custom runs with their plot styles."""
    cases = [
        ("scipy_scale_1", "SciPy, scale = 1", "scipy", 1.0, None, "black", "-"),
        (
            "scipy_current_scale",
            f"SciPy, scale = {current_scale:g}",
            "scipy",
            current_scale,
            None,
            "black",
            "--",
        ),
    ]
    colors = plt.get_cmap("plasma")(np.linspace(0.2, 0.8, 4))
    for fraction, color in zip(FIRST_STEP_FRACTIONS, colors):
        cases.append(
            (
                f"custom_step_{fraction:g}",
                f"Custom, initial fraction = {fraction:.1f}",
                "pysolver",
                1.0,
                fraction,
                color,
                "-",
            )
        )
    return cases


def evaluate_optimality_history(history, config, scale, run_dir):
    """Assess every recorded point independently, in that run's solver coordinates.

    A separate cycle and solver supply bounded numerical derivatives and the
    shared KKT multiplier fit. These postprocessing calls do not enter the
    optimization histories or their model-call totals.
    """
    problem = ComparisonCycleProblem(
        copy.deepcopy(config["problem_formulation"]), out_dir=str(run_dir)
    )
    problem.fitness(problem.x0)
    settings = config["solver_options"]
    diagnostics = psv.OptimizationSolver(
        problem,
        library="pysolver",
        method="slsqp",
        problem_scale=scale,
        derivative_method=settings["derivative_method"],
        derivative_abs_step=settings["derivative_abs_step"],
        print_convergence=False,
    )
    print(f"Assessing optimality at {len(history['x'])} recorded points", flush=True)
    values = []
    for point in history["x"]:
        try:
            assessment = diagnostics.evaluate_kkt_conditions(point, tol=FEASIBILITY_TOL)
            values.append(assessment["first_order_optimality"])
        except psv.DomainError as error:
            print(f"Optimality unavailable at one recorded point: {error}", flush=True)
            values.append(float("nan"))
    history["optimality"] = values
    # Include the initialization call made before the diagnostic solver was created.
    return 1 + diagnostics.func_count_tot


def run_case(case, config, output_dir, max_iterations=None, verbose=False):
    """Solve a fresh problem and save actual recorded observations and status."""
    name, label, library, scale, fraction, _, _ = case
    run_dir = output_dir / name
    problem = ComparisonCycleProblem(
        copy.deepcopy(config["problem_formulation"]), out_dir=str(run_dir)
    )
    # Initialize the constraint counts before constructing the solver.
    problem.fitness(problem.x0)
    options = copy.deepcopy(config["solver_options"])
    options.pop("callbacks", None)
    options.update(
        library=library,
        method="slsqp",
        problem_scale=scale,
        update_on="iteration",
        print_convergence=verbose,
        plot_convergence=False,
    )
    if max_iterations is not None:
        options["max_iterations"] = max_iterations
    # The YAML reader uses NumPy integers; custom SLSQP expects a Python int.
    options["max_iterations"] = int(options["max_iterations"])
    if library == "pysolver":
        options["extra_options"].update(
            first_step_fraction=fraction,
            max_step_fraction=1.0,
            step_growth_factor=1.5,
            feasibility_tol=FEASIBILITY_TOL,
            optimality_tol=OPTIMALITY_TOL,
        )
    solver = psv.OptimizationSolver(problem, **options)
    print(f"\nRunning {label} (problem scale = {scale:g})", flush=True)
    error_message = None
    try:
        solver.solve(problem.x0.copy())
    except psv.DomainError as error:
        # Native SciPy has no domain backtracking; retain its partial history.
        error_message = str(error)

    history = copy.deepcopy(solver.convergence_history)
    optimality_model_calls = evaluate_optimality_history(
        history, config, scale, run_dir
    )
    fields = (
        "iteration",
        "event_kind",
        "func_count_total",
        "objective_value",
        "constraint_violation",
        "optimality",
    )
    with (run_dir / "history.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        writer.writerows(zip(*(history[field] for field in fields)))
    summary = dict(
        name=name,
        label=label,
        library=library,
        problem_scale=float(scale),
        first_step_fraction=fraction,
        success=bool(solver.success) if error_message is None else False,
        message=error_message or str(solver.message),
        objective=getattr(solver, "f_final", None),
        model_calls=int(solver.func_count_tot),
        optimality_model_calls=int(optimality_model_calls),
    )
    (run_dir / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    status = "converged" if summary["success"] else "stopped"
    print(
        f"{status}: {summary['message']}\nModel calls: {summary['model_calls']}",
        flush=True,
    )
    return case, history, summary


def plot_comparison(results, output_dir):
    """Plot objective, feasibility, and optimality in one row with full axes boxes."""
    th.set_plot_options(fontsize=11, grid=False, linewidth=1.5)
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.8))
    for case, history, summary in results:
        _, label, _, _, _, color, linestyle = case
        if not summary["success"]:
            label += " (stopped)"
        for ax, field in zip(
            axes, ("objective_value", "constraint_violation", "optimality")
        ):
            ax.plot(
                history["iteration"],
                history[field],
                label=label,
                color=color,
                linestyle=linestyle,
                marker="o",
                markersize=3.2,
                markerfacecolor="white",
                markeredgewidth=0.9,
            )
    axes[0].set(
        title="(a) Objective function",
        ylabel=r"Objective, $f=-\eta_{\mathrm{cycle}}$ [-]",
    )
    axes[1].set(
        title="(b) Constraint violation", ylabel="Maximum normalized violation [-]"
    )
    axes[2].set(
        title="(c) Optimality", ylabel=r"Lagrangian gradient, $\|\nabla_x L\|_\infty$"
    )
    # Symlog preserves exact zeros without replacing them by invented positive values.
    for ax, tolerance, label in zip(
        axes[1:],
        (FEASIBILITY_TOL, OPTIMALITY_TOL),
        (r"$\varepsilon_c=10^{-6}$", r"$\varepsilon_{\mathrm{opt}}=10^{-4}$"),
    ):
        ax.set_yscale("symlog", linthresh=1e-8, linscale=0.5)
        ax.set_ylim(bottom=0)
        ax.axhline(tolerance, color="0.6", linestyle=":", linewidth=1)
        ax.annotate(
            label,
            xy=(1, tolerance),
            xycoords=("axes fraction", "data"),
            xytext=(-6, 4),
            textcoords="offset points",
            ha="right",
            fontsize=9,
            color="0.4",
        )
    for ax in axes:
        ax.set_xlabel("Solver iteration")
        ax.set_xlim(left=0)
        for spine in ax.spines.values():
            spine.set_visible(True)
        ax.tick_params(top=True, right=True, which="both")
        ax.grid(axis="y", color="0.9", linewidth=0.6)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False)
    fig.tight_layout(rect=(0, 0.17, 1, 1))
    for suffix in ("png", "svg", "pdf"):
        fig.savefig(
            output_dir / f"slsqp_comparison.{suffix}", dpi=300, bbox_inches="tight"
        )
    return fig


def main():
    """Run the six configurations and export their convergence comparison."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=EXAMPLE_DIR / "case_sCO2_recuperated.yaml"
    )
    parser.add_argument(
        "--output", type=Path, default=EXAMPLE_DIR / "results_slsqp_comparison"
    )
    parser.add_argument("--max-iterations", type=int, default=None)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()
    if args.max_iterations is not None and args.max_iterations < 1:
        parser.error("--max-iterations must be positive")
    config = th.read_configuration_file(str(args.config.resolve()))
    output_dir = args.output.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cases = comparison_cases(float(config["solver_options"]["problem_scale"]))
    results = [
        run_case(case, config, output_dir, args.max_iterations, args.verbose)
        for case in cases
    ]
    plot_comparison(results, output_dir)
    print(f"\nComparison saved to {output_dir}")
    if not args.no_show:
        plt.show()
    else:
        plt.close("all")


if __name__ == "__main__":
    main()
