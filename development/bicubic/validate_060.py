"""Compare native JaxProp evaluators through ThermOpt's eager cycle code.

Run with the desired JaxProp wheel installed. No model adapters, monkey patches,
outer JIT, altered cycle equations or silently clipped property calls are used.
Generated YAML files are ordinary ThermOpt inputs, usable without this runner.
"""

import argparse
import copy
import hashlib
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import statistics
import sys
from time import perf_counter

os.environ.setdefault("MPLBACKEND", "Agg")

import jax
import jaxprop as jp
import numpy as np
import yaml

import thermopt as th

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests"))
from regression_support import EXAMPLE_CASES, snapshot

# Separate interpolation allowances; the existing HEOS suite stays unchanged.
BUDGETS = dict(
    temperature_K=0.3,
    efficiency_absolute=5e-4,
    normalized_residual_absolute=0.02,
    density_relative=0.02,
)
CO2 = dict(h_min=110e3, h_max=1.3e6, p_min=8e6, p_max=39e6)
AUXILIARY = {
    "nitrogen": dict(h_min=150e3, h_max=1e6, p_min=90e3, p_max=115e3),
    "water": dict(h_min=60e3, h_max=180e3, p_min=90e3, p_max=115e3),
}


def plain(value):
    if hasattr(value, "__array__"):
        return np.asarray(value).tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, default=plain, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def configuration(case, method, grid, all_fluids, directory, table_directory=None):
    path = ROOT / "examples" / EXAMPLE_CASES[case]
    config = th.read_configuration_file(str(path))
    config["solver_options"].update(
        callbacks={}, plot_convergence=False, print_convergence=False
    )
    # Preserve the example's solver equations/options/tolerance/iteration budget.
    config["solver_options"]["max_iterations"] = int(
        config["solver_options"]["max_iterations"]
    )
    fixed = config["problem_formulation"]["fixed_parameters"]
    for role in ("working_fluid", "heating_fluid", "cooling_fluid"):
        settings = fixed[role]
        settings.pop("exceptions", None)
        settings["method"] = "coolprop"
        if method != "bicubic" or (role != "working_fluid" and not all_fluids):
            continue
        name = settings["name"].casefold()
        if role == "working_fluid":
            if name == "co2":
                domain = CO2
            elif name in ("nitrogen", "air"):
                domain = dict(h_min=100e3, h_max=700e3, p_min=40e3, p_max=1.2e6)
            else:
                raise ValueError(
                    "Declare a validated single-phase table domain for this working fluid first"
                )
            shape = (grid, grid)
        else:
            domain = AUXILIARY[name]
            shape = (60, 8)
        settings.update(
            method="bicubic",
            table_options=dict(
                **domain,
                N_h=shape[0],
                N_p=shape[1],
                N_p_sat=8,
                gradient_method="central",
                coarse_step=1,
                table_dir=str((table_directory or directory / "tables").resolve()),
            ),
        )
    return config, path


def make_manager(case, method, grid, all_fluids, directory, table_directory=None):
    config, source = configuration(
        case, method, grid, all_fluids, directory, table_directory
    )
    label = (
        method
        if method == "coolprop"
        else f"bicubic_{grid}" + ("_all" if all_fluids else "")
    )
    target = directory / case / label
    target.mkdir(parents=True, exist_ok=True)
    config_path = target / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(json.loads(json.dumps(config, default=plain)), sort_keys=False),
        encoding="utf-8",
    )
    start = perf_counter()
    manager = th.ThermodynamicCycleOptimization(str(config_path), out_dir=str(target))
    setup_seconds = perf_counter() - start
    for role, fluid in manager.problem.fluids.items():
        assert isinstance(fluid, jp.Fluid)
        if (
            config["problem_formulation"]["fixed_parameters"][role]["method"]
            == "bicubic"
        ):
            assert type(fluid) is jp.FluidBicubic
    return manager, dict(
        configuration=str(config_path.relative_to(ROOT)),
        original_configuration=str(source.relative_to(ROOT)),
        setup_seconds=setup_seconds,
        evaluator_classes={
            role: type(f).__name__ for role, f in manager.problem.fluids.items()
        },
    )


def temperatures(data):
    result = []
    density = []
    for component in data["components"].values():
        streams = (
            [component[s] for s in ("hot_side", "cold_side")]
            if "hot_side" in component
            else [component]
        )
        for stream in streams:
            result.extend(
                [float(stream[end]["T"]) for end in ("state_in", "state_out")]
            )
            density.extend(
                [float(stream[end]["rho"]) for end in ("state_in", "state_out")]
            )
            if "temperature_profile" in stream:
                result.extend(
                    np.asarray(stream["temperature_profile"]).ravel().tolist()
                )
        if "temperature_difference" in component:
            result.extend(component["temperature_difference"])
    return np.asarray(result), np.asarray(density)


def compare(actual, reference):
    temp, density = temperatures(actual)
    ref_temp, ref_density = temperatures(reference)
    residual_delta = np.abs(np.asarray(actual["residuals"]) - reference["residuals"])
    fields = [
        name
        for name in reference["energy_analysis"]
        if name in ("cycle_efficiency", "COP_refrigeration", "COP_heat_pump")
    ]
    efficiency_error = max(
        abs(float(actual["energy_analysis"][name]) - reference["energy_analysis"][name])
        for name in fields
    )
    values = dict(
        maximum_temperature_error_K=float(np.max(np.abs(temp - ref_temp))),
        maximum_density_relative_error=float(np.max(np.abs(density / ref_density - 1))),
        maximum_residual_error=float(np.max(residual_delta)),
        efficiency_absolute_error=efficiency_error,
    )
    assert values["maximum_temperature_error_K"] <= BUDGETS["temperature_K"], values
    assert (
        values["maximum_density_relative_error"] <= BUDGETS["density_relative"]
    ), values
    assert (
        values["maximum_residual_error"] <= BUDGETS["normalized_residual_absolute"]
    ), values
    assert efficiency_error <= BUDGETS["efficiency_absolute"], values
    return values


def timed(problem, point):
    start = perf_counter()
    result = problem.fitness(point)
    jax.block_until_ready(problem.cycle_data)
    assert np.all(np.isfinite(result))
    return perf_counter() - start


def feasibility(problem, point):
    lower, upper = map(np.asarray, problem.get_bounds())
    bound_error = np.maximum(lower - point, point - upper) / (upper - lower)
    return float(
        max(
            0,
            np.max(bound_error, initial=0),
            np.max(np.abs(problem.c_eq), initial=0),
            np.max(problem.c_ineq, initial=0),
        )
    )


def run(case, grids, all_fluids, optimize, samples, directory, table_directory=None):
    reference, setup = make_manager(case, "coolprop", 0, False, directory)
    point = reference.problem.x0.copy()
    source_initial = snapshot(reference.problem, point)
    result = dict(
        case=case,
        source_setup=setup,
        interpolation_budgets=BUDGETS,
        initial_source=source_initial,
        comparisons={},
    )
    managers = []
    for grid in grids:
        candidate, info = make_manager(
            case, "bicubic", grid, all_fluids, directory, table_directory
        )
        np.testing.assert_array_equal(candidate.problem.x0, point)
        initial = snapshot(candidate.problem, point)
        record = dict(
            setup=info, initial=initial, agreement=compare(initial, source_initial)
        )
        # Interleave complete, synchronized eager cycle evaluations after warmup.
        timed(reference.problem, point)
        timed(candidate.problem, point)
        observations = {"source": [], "bicubic": []}
        for sample in range(samples):
            order = ("source", "bicubic") if sample % 2 == 0 else ("bicubic", "source")
            for label in order:
                model = reference if label == "source" else candidate
                observations[label].append(timed(model.problem, point))
        record["warm_samples_seconds"] = observations
        record["warm_median_seconds"] = {
            k: statistics.median(v) for k, v in observations.items()
        }
        record["bicubic_to_source_time_ratio"] = (
            record["warm_median_seconds"]["bicubic"]
            / record["warm_median_seconds"]["source"]
        )
        result["comparisons"][str(grid)] = record
        managers.append((grid, candidate))
    if optimize:
        start = perf_counter()
        reference.run_optimization()
        source_seconds = perf_counter() - start
        assert reference.solver.success, reference.solver.message
        source_point = np.asarray(reference.solver.x_final).copy()
        source_final = snapshot(reference.problem, source_point)
        source_violation = feasibility(reference.problem, source_point)
        assert source_violation < 2e-6, source_violation
        result["source_optimization"] = dict(
            seconds=source_seconds,
            success=True,
            iterations=int(reference.solver.iteration_count),
            calls=int(reference.solver.func_count_tot),
            point=source_point,
            final=source_final,
            maximum_constraint_violation=source_violation,
        )
        for grid, candidate in managers:
            # First compare exactly the same source-optimum design.
            record = result["comparisons"][str(grid)]
            record["source_optimum_agreement"] = compare(
                snapshot(candidate.problem, source_point), source_final
            )
            start = perf_counter()
            candidate.run_optimization()
            seconds = perf_counter() - start
            assert candidate.solver.success, candidate.solver.message
            final_point = np.asarray(candidate.solver.x_final).copy()
            final = snapshot(candidate.problem, final_point)
            oracle = snapshot(reference.problem, final_point)
            agreement = compare(final, oracle)
            violation = feasibility(candidate.problem, final_point)
            source_violation = feasibility(reference.problem, final_point)
            assert violation < 2e-6, violation
            assert (
                source_violation < BUDGETS["normalized_residual_absolute"]
            ), source_violation
            objective_error = abs(
                float(final["residuals"][0]) - source_final["residuals"][0]
            )
            assert objective_error < BUDGETS["efficiency_absolute"], objective_error
            record["optimization"] = dict(
                seconds=seconds,
                success=True,
                iterations=int(candidate.solver.iteration_count),
                calls=int(candidate.solver.func_count_tot),
                point=final_point,
                final=final,
                independently_evaluated_source=oracle,
                same_design_agreement=agreement,
                maximum_constraint_violation=violation,
                source_constraint_violation=source_violation,
                objective_difference_from_source_optimum=objective_error,
            )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        choices=(
            "sco2_recuperated",
            "sco2_simple_polytropic",
            "sco2_recompression",
            "refrigeration_air",
        ),
        default="sco2_recuperated",
    )
    parser.add_argument("--grids", nargs="+", type=int, default=[100, 200])
    parser.add_argument("--all-fluids", action="store_true")
    parser.add_argument("--optimize", action="store_true")
    parser.add_argument("--samples", type=int, default=12)
    parser.add_argument(
        "--table-directory",
        type=Path,
        help="Reuse validated tables across separate output directories",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).parent / "results/compatibility_060/native",
    )
    args = parser.parse_args()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    context = dict(
        jaxprop_version=jp.__version__,
        jaxprop_path=str(Path(jp.__file__).resolve()),
        thermopt_path=str(Path(th.__file__).resolve()),
        dependencies={
            name: metadata.version(name)
            for name in (
                "jaxprop",
                "jax",
                "jaxlib",
                "CoolProp",
                "numpy",
                "scipy",
                "pysolver_view",
            )
        },
        numerical_source_sha256=hashlib.sha256(
            b"".join(p.read_bytes() for p in sorted((ROOT / "thermopt").rglob("*.py")))
        ).hexdigest(),
        execution="eager ThermOpt/NumPy/SciPy; native JaxProp objects; no application JIT or adapter",
        shared_machine_timings=True,
    )
    arguments = vars(args) | {
        "output": str(args.output),
        "table_directory": str(args.table_directory) if args.table_directory else None,
    }
    save(
        args.output / f"{args.case}_manifest.json",
        dict(context=context, arguments=arguments, budgets=BUDGETS),
    )
    result = run(
        args.case,
        args.grids,
        args.all_fluids,
        args.optimize,
        args.samples,
        args.output,
        args.table_directory,
    )
    save(
        args.output / f"{args.case}_summary.json",
        dict(context=context, **result, passed=True),
    )
    for grid, record in result["comparisons"].items():
        print(
            args.case,
            grid,
            record["agreement"],
            record["warm_median_seconds"],
            flush=True,
        )


if __name__ == "__main__":
    main()
