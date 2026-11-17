"""Alternate complete cycle evaluations from two ThermOpt source trees.

Both implementations run in independent package namespaces in the same process,
using installed JaxProp and identical ordinary configurations. No numerical
function or backend is patched. This complements fresh-process checkpoints when
shared-machine drift obscures a small interface change's warm evaluation cost.
"""

import argparse
import copy
import hashlib
import importlib.metadata as metadata
import importlib.util
import json
import os
from pathlib import Path
import platform
import statistics
import shutil
import sys
from time import perf_counter

from run import ROOT, canonical, cpu_affinity, save, sha


def source_hash(root):
    return hashlib.sha256(
        b"".join(
            path.read_bytes() for path in sorted((root / "thermopt").rglob("*.py"))
        )
    ).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=12)
    parser.add_argument("--cohorts", type=int, default=3)
    parser.add_argument("--cpu-affinity", type=lambda value: int(value, 0))
    parser.add_argument(
        "--workload",
        action="append",
        help="Restrict to named report workloads, e.g. simple_polytropic/bicubic",
    )
    parser.add_argument(
        "--first-version",
        choices=("reference", "current"),
        default="reference",
        help="Counterbalance construction and initial cache warmup",
    )
    parser.add_argument(
        "--identical-source-control",
        action="store_true",
        help="Use reference code for both variants, with separate model instances",
    )
    parser.add_argument(
        "--shared-models",
        action="store_true",
        help="Compare public functional cycle calls with the same explicitly supplied native models",
    )
    args = parser.parse_args()
    if args.samples < 4 or args.samples % 2 or args.cohorts < 1:
        parser.error("Use an even sample count >=4 and at least one cohort")
    reference_root = args.reference_root.resolve()
    if (
        reference_root == ROOT
        or not (reference_root / "thermopt/__init__.py").is_file()
    ):
        parser.error("Reference must be a separate ThermOpt source tree")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(__file__, args.output / "harness.py")
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
    import thermopt as current

    if not Path(jp.__file__).resolve().is_relative_to(Path(sys.prefix).resolve()):
        raise RuntimeError("Use installed JaxProp without source-checkout shadowing")
    if not jax.config.x64_enabled or jax.default_backend() != "cpu":
        raise RuntimeError("This benchmark requires CPU/x64 execution")
    assert Path(current.__file__).resolve().parent == ROOT / "thermopt"
    spec = importlib.util.spec_from_file_location(
        "thermopt_benchmark_reference",
        reference_root / "thermopt/__init__.py",
        submodule_search_locations=[str(reference_root / "thermopt")],
    )
    reference = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = reference
    spec.loader.exec_module(reference)
    current_package = reference if args.identical_source_control else current
    sys.path.insert(0, str(ROOT))
    sys.path.insert(0, str(ROOT / "tests"))
    from benchmark.native_backends import configuration, plain
    from regression_support import load_case, snapshot

    def output_signature(value):
        """Compare complete returned data, normalizing only the path container.

        Raw optional NaNs, status and identity are retained through the public
        state serialization interface. Projection work stays outside timings.
        """
        from collections.abc import Mapping

        def project(item, key=None):
            if key == "states" and isinstance(item, list):
                item = jp.FluidState.stack(item)
            if isinstance(item, jp.FluidState):
                return project(item.to_data())
            if isinstance(item, jp.Fluid):
                return {
                    name: repr(getattr(item, name))
                    for name in (
                        "name",
                        "provider",
                        "model",
                        "identifier",
                        "reference_id",
                    )
                }
            if isinstance(item, Mapping):
                return {name: project(child, name) for name, child in item.items()}
            if isinstance(item, (tuple, list)):
                return [project(child) for child in item]
            if isinstance(item, (set, frozenset)):
                return sorted(project(child) for child in item)
            if hasattr(item, "__array__"):
                return np.asarray(item).tolist()
            return item

        data = json.dumps(project(value), sort_keys=True, allow_nan=True)
        return hashlib.sha256(data.encode()).hexdigest()

    current_source_root = reference_root if args.identical_source_control else ROOT
    identities = dict(
        reference=source_hash(reference_root), current=source_hash(current_source_root)
    )
    report = dict(
        reference_root=str(reference_root),
        current_root=str(current_source_root),
        source_sha256=identities,
        harness_sha256=sha(__file__),
        protocol_sha256=sha(ROOT / "benchmark/run.py"),
        workload_sha256=sha(ROOT / "benchmark/native_backends.py"),
        helper_sha256=sha(ROOT / "tests/regression_support.py"),
        lock_sha256=sha(ROOT / "poetry.lock"),
        environment=dict(
            python=platform.python_version(),
            platform=platform.platform(),
            processor=platform.processor(),
            cpu_affinity=affinity,
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
        ),
        samples=args.samples,
        cohorts=args.cohorts,
        first_version=args.first_version,
        identical_source_control=args.identical_source_control,
        shared_models=args.shared_models,
        warm_ratio_limit=1.10,
        interpretation="Same-process alternating complete evaluations; separate native models, shared installed JaxProp and CPU pool. Warm only; compilation/setup recorded separately. Fresh-process evidence remains necessary context.",
        results={},
    )
    destination = args.output / "comparison.json"
    save(destination, report)
    jobs = [
        ("recuperated_working", "sco2_recuperated", False),
        ("recuperated_all", "sco2_recuperated", True),
        ("simple_polytropic", "sco2_simple_polytropic", False),
        ("recompression", "sco2_recompression", False),
        ("refrigeration", "refrigeration_air", True),
    ]
    configurations = []
    for label, case, all_streams in jobs:
        for method in ("coolprop", "bicubic"):
            config, _ = configuration(
                case,
                method,
                100,
                all_streams,
                args.output,
                ROOT / "benchmark/results/tables",
            )
            configurations.append((f"{label}/{method}", config["problem_formulation"]))
    configurations += [
        (case, load_case(case)["problem_formulation"])
        for case in ("ptes", "ptes_turbo")
    ]
    if args.workload:
        configurations = [
            (label, config)
            for label, config in configurations
            if label in args.workload
        ]
        if {label for label, _ in configurations} != set(args.workload):
            parser.error("Unknown workload selector")

    for label, config in configurations:
        problems, setup = {}, {}
        packages = {"reference": reference, "current": current_package}
        names = (
            ("reference", "current")
            if args.first_version == "reference"
            else ("current", "reference")
        )
        for name in names:
            package = packages[name]
            start = perf_counter()
            problems[name] = package.ThermodynamicCycleProblem(
                copy.deepcopy(config), out_dir=str(args.output / label / name)
            )
            setup[name] = perf_counter() - start
        baseline = snapshot(problems[names[0]], problems[names[0]].x0)
        assert snapshot(problems[names[1]], problems[names[1]].x0) == baseline, label
        shared = problems[names[0]].fluids

        def evaluate(name):
            problem = problems[name]
            if not args.shared_models:
                values = problem.fitness(problem.x0)
                return values, problem.cycle_data
            # This is the existing public functional API, with complete return
            # data and explicit native evaluators; no problem is patched.
            result = packages[name].CYCLE_TOPOLOGIES[problem.cycle_topology](
                dict(zip(problem.variable_names, problem.x0)),
                problem.fixed_parameters,
                problem.constraints,
                problem.objective_function,
                fluids=shared,
            )
            return None, result

        if args.shared_models:
            functional_signature = None
            for name in names:
                _, result = evaluate(name)
                jax.block_until_ready(result)
                signature = output_signature(result)
                if functional_signature is None:
                    functional_signature = signature
                assert signature == functional_signature, label
        bursts = []
        for cohort in range(args.cohorts):
            samples = {"reference": [], "current": []}
            for index in range(args.samples):
                order = names if (index + cohort) % 2 == 0 else names[::-1]
                for name in order:
                    start = perf_counter()
                    values, result = evaluate(name)
                    jax.block_until_ready(result)
                    samples[name].append(perf_counter() - start)
                    if args.shared_models:
                        values = np.hstack(
                            [
                                np.atleast_1d(result[key])
                                for key in (
                                    "objective_function",
                                    "equality_constraints",
                                    "inequality_constraints",
                                )
                            ]
                        )
                    assert np.all(np.isfinite(values))
            medians = {
                name: statistics.median(values) for name, values in samples.items()
            }
            bursts.append(
                dict(
                    raw_samples_seconds=samples,
                    medians_seconds=medians,
                    ratio=medians["current"] / medians["reference"],
                )
            )
        for problem in problems.values():
            assert snapshot(problem, problem.x0) == baseline, label
        if args.shared_models:
            for name in names:
                _, result = evaluate(name)
                assert output_signature(result) == functional_signature, label
        medians = {
            name: statistics.median(burst["medians_seconds"][name] for burst in bursts)
            for name in problems
        }
        ratio = medians["current"] / medians["reference"]
        report["results"][label] = dict(
            input=json.loads(
                json.dumps(canonical(config), default=plain, allow_nan=False)
            ),
            setup_seconds=setup,
            observations=bursts,
            aggregate_seconds=medians,
            ratio=ratio,
            timing_flag=ratio > 1.10,
            numerical_snapshot_sha256=hashlib.sha256(
                repr(baseline).encode()
            ).hexdigest(),
            physical_snapshots_identical=True,
            complete_shared_output_identical=True if args.shared_models else None,
        )
        save(destination, report)
        print(label, ratio, flush=True)
    assert identities == dict(
        reference=source_hash(reference_root), current=source_hash(current_source_root)
    )
    assert cpu_affinity()["actual_mask"] == affinity["actual_mask"]
    assert sha(__file__) == report["harness_sha256"]
    assert sha(ROOT / "tests/regression_support.py") == report["helper_sha256"]
    assert sha(ROOT / "poetry.lock") == report["lock_sha256"]
    report["source_and_affinity_unchanged"] = True
    report["timing_flags"] = [
        name for name, value in report["results"].items() if value["timing_flag"]
    ]
    save(destination, report)
    if report["timing_flags"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
