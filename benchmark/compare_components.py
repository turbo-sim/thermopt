"""Alternate complete eager component returns from two ThermOpt source trees.

Use installed JaxProp, independent native models, identical physical inputs and
the same CPU/x64 profile for both versions. All state properties, identities,
availability masks and statuses are compared outside timings. No projected
property kernels, patched backends or timing assertions in pytest are used.
"""

import argparse
from collections.abc import Mapping
import importlib.metadata as metadata
import importlib.util
import os
from pathlib import Path
import shutil
import statistics
import sys
from time import perf_counter

from run import ROOT, canonical, cpu_affinity, save, sha
from compare_versions import source_hash


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--samples", type=int, default=12)
    parser.add_argument("--cohorts", type=int, default=3)
    parser.add_argument("--cpu-affinity", type=lambda value: int(value, 0))
    parser.add_argument(
        "--first-version", choices=("reference", "current"), default="reference"
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

    if not Path(jp.__file__).resolve().is_relative_to(Path(sys.prefix).resolve()):
        raise RuntimeError("Use installed JaxProp without checkout shadowing")
    if not jax.config.x64_enabled or jax.default_backend() != "cpu":
        raise RuntimeError("This benchmark requires CPU/x64 execution")
    spec = importlib.util.spec_from_file_location(
        "thermopt_component_reference",
        reference_root / "thermopt/__init__.py",
        submodule_search_locations=[str(reference_root / "thermopt")],
    )
    reference = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = reference
    spec.loader.exec_module(reference)
    versions = dict(reference=reference, current=thermopt)
    order = (
        ("reference", "current")
        if args.first_version == "reference"
        else ("current", "reference")
    )
    sources = dict(reference=source_hash(reference_root), current=source_hash(ROOT))
    jp_root = Path(jp.__file__).resolve().parent
    jp_sources = {
        str(p.relative_to(jp_root)).replace("\\", "/"): sha(p)
        for p in sorted(jp_root.rglob("*.py"))
    }
    report = dict(
        schema=1,
        source_sha256=sources,
        jaxprop_source_sha256=jp_sources,
        harness_sha256=sha(__file__),
        lock_sha256=sha(ROOT / "poetry.lock"),
        cpu_affinity=affinity,
        samples=args.samples,
        cohorts=args.cohorts,
        first_version=args.first_version,
        warm_ratio_limit=1.10,
        environment=dict(
            python=sys.version,
            executable=sys.executable,
            x64=jax.config.x64_enabled,
            device=str(jax.devices()[0]),
            dependencies={
                name: metadata.version(name)
                for name in (
                    "jaxprop",
                    "jax",
                    "jaxlib",
                    "numpy",
                    "scipy",
                    "CoolProp",
                    "equinox",
                )
            },
        ),
        numerical_tolerance=dict(rtol=2e-10, atol=1e-7, equal_nan=True),
        interpretation="Complete native component returns, independent reusable models; cold/setup costs recorded separately. Shared-machine warm flags require repeats.",
        results={},
        model_setup={},
    )
    destination = output / "comparison.json"
    save(destination, report)

    def compare(actual, expected, path="result"):
        if isinstance(expected, jp.FluidState):
            assert isinstance(actual, jp.FluidState), path
            compare(actual.to_data(), expected.to_data(), path)
        elif isinstance(expected, Mapping):
            assert actual.keys() == expected.keys(), path
            for key in expected:
                compare(actual[key], expected[key], f"{path}/{key}")
        elif isinstance(expected, (tuple, list)):
            assert type(actual) is type(expected) and len(actual) == len(expected), path
            for index, (a, b) in enumerate(zip(actual, expected)):
                compare(a, b, f"{path}/{index}")
        elif expected is None or isinstance(expected, (str, set, frozenset)):
            assert actual == expected, path
        else:
            a, b = np.asarray(actual), np.asarray(expected)
            assert a.shape == b.shape, path
            if b.dtype.kind not in "fc":
                np.testing.assert_array_equal(a, b, err_msg=path)
            else:
                np.testing.assert_allclose(
                    a, b, rtol=2e-10, atol=1e-7, equal_nan=True, err_msg=path
                )

    tables = dict(
        h_min=200e3,
        h_max=850e3,
        p_min=40e3,
        p_max=350e3,
        N_h=32,
        N_p=32,
        N_p_sat=4,
        coarse_step=1,
        table_dir=str(ROOT / "benchmark/results/phase4/air_tables"),
    )
    for method in ("coolprop", "perfect_gas", "bicubic"):
        settings = dict(name="Air", method=method, identifier="benchmark_air")
        if method == "bicubic":
            settings["table_options"] = tables
        models = {}
        for version in order:
            start = perf_counter()
            models[version] = jp.Fluid(**settings)
            report["model_setup"][f"{method}/{version}"] = perf_counter() - start
        identities = {key: id(value) for key, value in models.items()}
        inlet = models[order[0]].get_state(jp.PT_INPUTS, 1e5, 300.0)
        turbine_inlet = models[order[0]].get_state(jp.PT_INPUTS, 2e5, 600.0)
        h_compressor, h_expander = float(inlet.enthalpy), float(turbine_inlet.enthalpy)
        jobs = [
            (
                f"{machine}/{kind}",
                machine,
                dict(
                    h_in=(
                        h_compressor if machine == "compression_process" else h_expander
                    ),
                    p_in=1e5 if machine == "compression_process" else 2e5,
                    p_out=2e5 if machine == "compression_process" else 1e5,
                    efficiency=0.85,
                    efficiency_type=kind,
                    mass_flow=2.0,
                ),
            )
            for machine in ("compression_process", "expansion_process")
            for kind in ("isentropic", "polytropic", "non-dimensional")
        ]
        jobs += [
            (
                "heat_exchanger",
                "heat_exchanger",
                dict(
                    h_in_hot=750e3,
                    h_out_hot=600e3,
                    p_in_hot=2e5,
                    p_out_hot=1.9e5,
                    h_in_cold=430e3,
                    h_out_cold=500e3,
                    p_in_cold=1e5,
                    p_out_cold=0.95e5,
                    num_steps=50,
                ),
            ),
            ("valve", "isenthalpic_valve", dict(p_out=1e5, N=50)),
        ]
        for label, function, inputs in jobs:
            key = f"{method}/{label}"

            def evaluate(version):
                module = sys.modules[
                    f"{versions[version].__name__}.components.basic_components"
                ]
                fluid = models[version]
                if function == "heat_exchanger":
                    return module.heat_exchanger(
                        fluid_hot=fluid, fluid_cold=fluid, **inputs
                    )
                if function == "isenthalpic_valve":
                    return module.isenthalpic_valve(
                        turbine_inlet, fluid=fluid, **inputs
                    )
                return getattr(module, function)(fluid=fluid, **inputs)

            cold, baseline = {}, None
            for version in order:
                start = perf_counter()
                result = evaluate(version)
                jax.block_until_ready(result)
                cold[version] = perf_counter() - start
                if baseline is None:
                    baseline = result
                else:
                    compare(result, baseline)
            observations = []
            for cohort in range(args.cohorts):
                samples = {version: [] for version in order}
                for index in range(args.samples):
                    for version in (
                        order if (cohort + index) % 2 == 0 else order[::-1]
                    ):
                        start = perf_counter()
                        result = evaluate(version)
                        jax.block_until_ready(result)
                        samples[version].append(perf_counter() - start)
                    compare(result, baseline)
                medians = {
                    name: statistics.median(values) for name, values in samples.items()
                }
                observations.append(
                    dict(
                        raw_samples_seconds=samples,
                        medians_seconds=medians,
                        ratio=medians["current"] / medians["reference"],
                    )
                )
            for version in order:
                compare(evaluate(version), baseline)
            medians = {
                name: statistics.median(
                    b["medians_seconds"][name] for b in observations
                )
                for name in order
            }
            ratio = medians["current"] / medians["reference"]
            report["results"][key] = dict(
                settings=canonical(settings),
                inputs=inputs,
                first_call_seconds=cold,
                observations=observations,
                aggregate_seconds=medians,
                ratio=ratio,
                timing_flag=ratio > 1.10,
                complete_outputs_agree=True,
            )
            save(destination, report)
            print(key, f"{ratio:.4f}", flush=True)
        assert identities == {key: id(value) for key, value in models.items()}
    assert sources == dict(
        reference=source_hash(reference_root), current=source_hash(ROOT)
    )
    assert cpu_affinity()["actual_mask"] == affinity["actual_mask"]
    assert sha(__file__) == report["harness_sha256"]
    assert sha(ROOT / "poetry.lock") == report["lock_sha256"]
    assert jp_sources == {
        str(p.relative_to(jp_root)).replace("\\", "/"): sha(p)
        for p in sorted(jp_root.rglob("*.py"))
    }
    report["source_and_protocol_unchanged"] = True
    report["numerical_checks_passed"] = True
    report["timing_flags"] = [
        key for key, value in report["results"].items() if value["timing_flag"]
    ]
    save(destination, report)
    if report["timing_flags"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
