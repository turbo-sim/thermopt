"""Time complete solver-boundary evaluations separately from reporting I/O.

The base route is the inherited PySolverView method used before Phase 6. Both
routes evaluate the same complete cycle with one configured problem, avoiding
independent-model JIT-cache-order differences. This measures only the added
exception boundary, not a comparison with a different JaxProp release. Output
timings describe the current, numerically verified writers, not the old writers
which could emit text-valued energy cells or incomplete manager configurations.
"""

import argparse
import hashlib
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import platform
import statistics
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=12)
    parser.add_argument("--report-samples", type=int, default=4)
    parser.add_argument("--grid", type=int, default=32)
    args = parser.parse_args()
    if (
        args.samples < 4
        or args.samples % 2
        or args.report_samples < 1
        or args.grid < 16
    ):
        parser.error("Use an even sample count >=4, report samples >=1 and grid >=16")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    os.environ.update(
        MPLBACKEND="Agg",
        DISABLE_PLOTS="1",
        JAX_ENABLE_X64="true",
        JAX_PLATFORM_NAME="cpu",
        OMP_NUM_THREADS="1",
    )
    sys.path[:0] = [str(ROOT), str(ROOT / "tests")]
    import jax
    import jaxprop as jp
    import matplotlib.pyplot as plt
    import numpy as np
    import pysolver_view as psv
    import thermopt as th
    import yaml
    from regression_support import load_case
    from run import cpu_affinity

    assert jax.config.x64_enabled and jax.default_backend() == "cpu"
    assert Path(jp.__file__).resolve().is_relative_to(Path(sys.prefix).resolve())
    report = dict(
        scope="Phase 6 whole-cycle normalized-input exception boundary and reporting",
        baseline="Inherited PySolverView OptimizationProblem.fitness_normalized_input",
        control="One problem/model instance per provider; alternate route order; full cycle evaluated",
        packages={
            name: metadata.version(name)
            for name in (
                "jaxprop",
                "pysolver-view",
                "pounce-solver",
                "jax",
                "jaxlib",
                "numpy",
                "scipy",
            )
        },
        source_sha256={
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted((ROOT / "thermopt").rglob("*.py"))
        },
        harness_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        samples=args.samples,
        reporting_samples=args.report_samples,
        grid=args.grid,
        warm_ratio_limit=1.10,
        environment=dict(
            python=platform.python_version(),
            platform=platform.platform(),
            cpu_affinity=cpu_affinity(),
            hash_seed=os.environ.get("PYTHONHASHSEED"),
        ),
        workloads={},
    )
    for method in ("coolprop", "bicubic"):
        config = load_case("sco2_recuperated")
        config["solver_options"].update(
            callbacks={}, print_convergence=False, plot_convergence=False
        )
        formulation = config["problem_formulation"]
        formulation["plot_settings"].update(
            fluid=dict(plot_saturation_line=False, plot_critical_point=False)
        )
        if method == "bicubic":
            formulation["fixed_parameters"]["working_fluid"].update(
                method=method,
                table_options=dict(
                    h_min=110e3,
                    h_max=1.3e6,
                    p_min=8e6,
                    p_max=39e6,
                    N_h=args.grid,
                    N_p=args.grid,
                    N_p_sat=4,
                    coarse_step=1,
                    table_dir=str(output / "tables"),
                ),
            )
        folder = output / method
        folder.mkdir()
        source = folder / "input.yaml"
        source.write_text(
            yaml.safe_dump(
                th.convert_numpy_to_python(config, precision=None), sort_keys=False
            ),
            encoding="utf-8",
        )
        manager = th.ThermodynamicCycleOptimization(
            source, out_dir=str(folder / "results")
        )
        problem = manager.problem
        point = problem.scale_physical_to_normalized(problem.x0)
        routes = dict(
            inherited=lambda: psv.OptimizationProblem.fitness_normalized_input(
                problem, point
            ),
            current=lambda: problem.fitness_normalized_input(point),
        )
        baseline = routes["inherited"]()
        np.testing.assert_array_equal(routes["current"](), baseline)
        for _ in range(3):
            for route in routes.values():
                route()
        samples = {name: [] for name in routes}
        for index in range(args.samples):
            order = list(routes) if index % 2 == 0 else list(reversed(routes))
            for name in order:
                start = perf_counter()
                actual = routes[name]()
                samples[name].append(perf_counter() - start)
                np.testing.assert_array_equal(actual, baseline)
        medians = {name: statistics.median(values) for name, values in samples.items()}
        manager.solver.fitness(point)
        calls_before_report = manager.solver.func_count_tot
        start = perf_counter()
        manager.save_report_callback(point, 0)
        first_report = perf_counter() - start
        first_report_calls = manager.solver.func_count_tot - calls_before_report
        reporting = {
            name: [] for name in ("yaml", "excel", "plot_update", "cached_report")
        }
        problem.plot_cycle()
        for index in range(args.report_samples):
            operations = dict(
                yaml=lambda: manager.save_current_configuration(
                    folder / "snapshot.yaml"
                ),
                excel=lambda: problem.save_data_to_excel(folder / "cycle.xlsx"),
                plot_update=problem.plot_cycle,
                cached_report=lambda: manager.save_report_callback(point, index + 1),
            )
            for name, operation in operations.items():
                start = perf_counter()
                operation()
                reporting[name].append(perf_counter() - start)
        plt.close("all")
        ratio = medians["current"] / medians["inherited"]
        report["workloads"][method] = dict(
            native_provider=problem.fluid.provider,
            samples_seconds=samples,
            medians_seconds=medians,
            ratio=ratio,
            timing_flag=ratio > report["warm_ratio_limit"],
            numerical_equal=True,
            reporting_samples_seconds=reporting,
            first_report_seconds=first_report,
            first_report_solver_tracked_calls=first_report_calls,
            reporting_note="Elapsed report time includes I/O-only design synchronization; solver counters exclude those direct refreshes",
            reporting_medians_seconds={
                name: statistics.median(values) for name, values in reporting.items()
            },
        )
        print(
            f"{method}: inherited {medians['inherited']*1000:.3f} ms, "
            f"current {medians['current']*1000:.3f} ms, ratio {ratio:.3f}",
            flush=True,
        )
        (output / "report.json").write_text(
            json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )


if __name__ == "__main__":
    main()
