"""Compare configured cycle construction and complete reload/evaluation costs.

Use independent native models per version, explicit configuration inputs and
the installed JaxProp package. Model identities and complete physical cycle
snapshots are checked outside timings. No table metadata work is omitted from
reload timings, and no source-identity checks enter fitness timing itself.
"""

import argparse
from copy import deepcopy
import importlib.metadata as metadata
import importlib.util
import json
import os
from pathlib import Path
import shutil
import statistics
import sys
from time import perf_counter

from run import ROOT, canonical, cpu_affinity, save, sha
from compare_versions import source_hash
from first_use import compare, plain


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=12)
    parser.add_argument("--cohorts", type=int, default=3)
    parser.add_argument("--cpu-affinity", type=lambda v: int(v, 0))
    args = parser.parse_args()
    if args.samples < 4 or args.samples % 2 or args.cohorts < 1:
        parser.error("Use even samples >=4 and at least one cohort")
    reference_root = args.reference_root.resolve()
    if reference_root == ROOT:
        parser.error("Specify a separate reference tree")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(__file__, output / "harness.py")
    affinity = cpu_affinity(args.cpu_affinity)
    os.environ.update(
        JAX_ENABLE_X64="true",
        JAX_PLATFORM_NAME="cpu",
        MPLBACKEND="Agg",
        OMP_NUM_THREADS="1",
    )
    import jax
    import jaxprop as jp
    import numpy as np
    import thermopt

    assert Path(jp.__file__).resolve().is_relative_to(Path(sys.prefix).resolve())
    assert jax.default_backend() == "cpu" and jax.config.x64_enabled
    spec = importlib.util.spec_from_file_location(
        "thermopt_lifecycle_reference",
        reference_root / "thermopt/__init__.py",
        submodule_search_locations=[str(reference_root / "thermopt")],
    )
    reference = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = reference
    spec.loader.exec_module(reference)
    packages = dict(reference=reference, current=thermopt)
    identities = dict(reference=source_hash(reference_root), current=source_hash(ROOT))
    sys.path.insert(0, str(ROOT))
    sys.path.insert(0, str(ROOT / "tests"))
    from benchmark.native_backends import configuration
    from regression_support import load_case

    configs = {}
    for method in ("coolprop", "bicubic"):
        config, _ = configuration(
            "sco2_recuperated",
            method,
            100,
            False,
            output,
            ROOT / "benchmark/results/tables",
        )
        configs[method] = config["problem_formulation"]
    ideal = deepcopy(load_case("sco2_simple_polytropic")["problem_formulation"])
    ideal["fixed_parameters"]["working_fluid"] = dict(
        name="Air", method="perfect_gas", model_options=dict(T_ref=600.0, p_ref=1e5)
    )
    model = jp.Fluid(**ideal["fixed_parameters"]["working_fluid"])
    values = dict(
        compressor_inlet_pressure=1e5,
        expander_inlet_pressure=2e5,
        compressor_inlet_enthalpy=float(
            model.get_state(jp.PT_INPUTS, 1e5, 300.0).enthalpy
        ),
        expander_inlet_enthalpy=float(
            model.get_state(jp.PT_INPUTS, 2e5, 750.0).enthalpy
        ),
    )
    for key, value in values.items():
        ideal["design_variables"][key] = dict(
            min=0.8 * value, value=value, max=1.2 * value
        )
    configs["perfect_gas"] = ideal
    jp_root = Path(jp.__file__).resolve().parent
    jp_hashes = {
        str(p.relative_to(jp_root)).replace("\\", "/"): sha(p)
        for p in sorted(jp_root.rglob("*.py"))
    }
    report = dict(
        schema=1,
        harness_sha256=sha(__file__),
        lock_sha256=sha(ROOT / "poetry.lock"),
        source_sha256=identities,
        jaxprop_source_sha256=jp_hashes,
        cpu_affinity=affinity,
        samples=args.samples,
        cohorts=args.cohorts,
        inputs=plain({key: canonical(value) for key, value in configs.items()}),
        environment=dict(
            python=sys.version,
            executable=sys.executable,
            jaxprop_path=jp.__file__,
            dependencies={
                name: metadata.version(name)
                for name in (
                    "jaxprop",
                    "jax",
                    "jaxlib",
                    "CoolProp",
                    "numpy",
                    "scipy",
                    "equinox",
                )
            },
        ),
        interpretation="Synchronized complete cycle returns. Reload setup includes all provenance/rendering work; constructor compilation/warmup and fitness are recorded separately. Shared-machine observations, no new stopwatch test assertions.",
        results={},
    )
    destination = output / "comparison.json"
    save(destination, report)
    for method, config in configs.items():
        problems = {}
        initial = {}
        cold = {}
        for version, package in packages.items():
            start = perf_counter()
            problems[version] = package.ThermodynamicCycleProblem(
                deepcopy(config), out_dir=str(output / method / version)
            )
            setup = perf_counter() - start
            start = perf_counter()
            values = problems[version].fitness(problems[version].x0)
            jax.block_until_ready(problems[version].cycle_data)
            cold[version] = dict(
                setup_seconds=setup, first_cycle_seconds=perf_counter() - start
            )
            assert np.all(np.isfinite(values))
            initial[version] = plain(problems[version].cycle_data)
        compare(initial["current"], initial["reference"])
        for mode in ("unchanged", "nonfluid"):
            records = []
            for cohort in range(args.cohorts):
                observations = {
                    name: dict(reload_seconds=[], fitness_seconds=[], total_seconds=[])
                    for name in packages
                }
                for index in range(args.samples):
                    candidate = deepcopy(config)
                    if mode == "nonfluid":
                        candidate["fixed_parameters"]["net_power"] = 1e6 * (
                            1 + 0.1 * (index % 2)
                        )
                    for version in (
                        ("reference", "current")
                        if (index + cohort) % 2 == 0
                        else ("current", "reference")
                    ):
                        problem = problems[version]
                        previous = dict(problem.fluids)
                        start = perf_counter()
                        problem.update_problem(candidate)
                        reload_seconds = perf_counter() - start
                        start = perf_counter()
                        values = problem.fitness(problem.x0)
                        jax.block_until_ready(problem.cycle_data)
                        fitness_seconds = perf_counter() - start
                        assert np.all(np.isfinite(values))
                        if version == "current":
                            assert all(
                                problem.fluids[key] is model
                                for key, model in previous.items()
                            )
                        for key, value in (
                            ("reload_seconds", reload_seconds),
                            ("fitness_seconds", fitness_seconds),
                            ("total_seconds", reload_seconds + fitness_seconds),
                        ):
                            observations[version][key].append(value)
                    # Check full physical outputs after both evaluations.
                    compare(
                        plain(problems["current"].cycle_data),
                        plain(problems["reference"].cycle_data),
                    )
                medians = {
                    name: {
                        key: statistics.median(samples)
                        for key, samples in measurements.items()
                    }
                    for name, measurements in observations.items()
                }
                records.append(
                    dict(raw_samples_seconds=observations, medians_seconds=medians)
                )
            medians = {
                name: {
                    key: statistics.median(
                        cohort["medians_seconds"][name][key] for cohort in records
                    )
                    for key in ("reload_seconds", "fitness_seconds", "total_seconds")
                }
                for name in packages
            }
            ratios = {
                key: medians["current"][key] / medians["reference"][key]
                for key in medians["reference"]
            }
            report["results"][f"{method}/{mode}"] = dict(
                observations=records,
                medians=medians,
                ratios=ratios,
                construction_and_first_cycle=cold,
                physical_snapshots_agree=True,
                models_reused=True,
            )
            save(destination, report)
            print(method, mode, ratios, flush=True)
    assert identities == dict(
        reference=source_hash(reference_root), current=source_hash(ROOT)
    )
    assert (
        sha(__file__) == report["harness_sha256"]
        and sha(ROOT / "poetry.lock") == report["lock_sha256"]
    )
    assert jp_hashes == {
        str(p.relative_to(jp_root)).replace("\\", "/"): sha(p)
        for p in sorted(jp_root.rglob("*.py"))
    }
    assert cpu_affinity()["actual_mask"] == affinity["actual_mask"]
    report["source_and_protocol_unchanged"] = True
    report["numerical_checks_passed"] = True
    save(destination, report)


if __name__ == "__main__":
    main()
