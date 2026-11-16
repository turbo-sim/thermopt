"""Measure polytropic first use in independent CPU/x64 worker processes.

Separate import/model setup, inlet evaluation, explicit warmup, the first
complete component and sustained complete returns. Each worker owns one native
model and one source implementation. No compiler cache is shared between
versions. Optional compilation logs are diagnostic runs, not warm timing gates.
"""

import argparse
from contextlib import nullcontext
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
from time import perf_counter

from run import ROOT, cpu_affinity, save, sha
from compare_versions import source_hash


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--reference-root", type=Path)
    result.add_argument("--output", type=Path, required=True)
    result.add_argument("--samples", type=int, default=12)
    result.add_argument("--cohorts", type=int, default=3)
    result.add_argument("--cpu-affinity", type=lambda value: int(value, 0))
    result.add_argument("--model", choices=("perfect_gas", "bicubic"), action="append")
    result.add_argument(
        "--machine", choices=("compression", "expansion"), action="append"
    )
    result.add_argument("--warmup", choices=("none", "scalar", "path"), action="append")
    result.add_argument(
        "--table-directory",
        type=Path,
        default=ROOT / "benchmark/results/first_use_tables",
    )
    result.add_argument("--log-compiles", action="store_true")
    result.add_argument("--worker-root", type=Path, help=argparse.SUPPRESS)
    return result


def compare(actual, expected, path="result"):
    """Check exported full states, including metadata, NaNs and validity masks."""
    import numpy as np

    if isinstance(expected, dict):
        assert actual.keys() == expected.keys(), path
        for key in expected:
            compare(actual[key], expected[key], f"{path}/{key}")
    elif isinstance(expected, list):
        assert isinstance(actual, list) and len(actual) == len(expected), path
        for index, (a, b) in enumerate(zip(actual, expected)):
            compare(a, b, f"{path}/{index}")
    elif expected is None or isinstance(expected, (str, bool, int)):
        assert actual == expected, path
    else:
        np.testing.assert_allclose(
            actual, expected, rtol=2e-10, atol=1e-7, err_msg=path
        )


def plain(value):
    """Export full public state data without requesting unsupported properties."""
    from collections.abc import Mapping
    import jax
    import jaxprop as jp
    import numpy as np

    if isinstance(value, jp.FluidState):
        return plain(value.to_data())
    if isinstance(value, jp.Fluid):
        config = value.configuration
        return plain(
            dict(
                configuration={
                    name: getattr(config, name)
                    for name in (
                        "name",
                        "provider",
                        "model",
                        "backend",
                        "identifier",
                        "table_options",
                        "model_options",
                    )
                },
                reference_id=value.reference_id,
            )
        )
    if isinstance(value, Mapping):
        return {key: plain(item) for key, item in value.items()}
    if isinstance(value, (set, frozenset)):
        return [plain(item) for item in sorted(value)]
    if isinstance(value, (list, tuple)):
        return [plain(item) for item in value]
    if isinstance(value, np.ndarray) or isinstance(value, jax.Array):
        return plain(np.asarray(value).tolist())
    if isinstance(value, np.generic):
        return plain(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return repr(value)
    return value


def worker(args):
    affinity = cpu_affinity(args.cpu_affinity)
    os.environ.update(
        JAX_ENABLE_X64="true",
        JAX_PLATFORM_NAME="cpu",
        MPLBACKEND="Agg",
        OMP_NUM_THREADS="1",
    )
    root = args.worker_root.resolve()
    sys.path.insert(0, str(root))
    start = perf_counter()
    import jax
    import jaxprop as jp
    import numpy as np
    from thermopt.components import basic_components as components
    from thermopt.utilities import process_state_options

    imports = perf_counter() - start
    jax.config.update("jax_enable_compilation_cache", False)
    assert Path(components.__file__).resolve().is_relative_to(root / "thermopt")
    assert Path(jp.__file__).resolve().is_relative_to(Path(sys.prefix).resolve())
    assert jax.default_backend() == "cpu" and jax.config.x64_enabled
    identity = source_hash(root)
    settings = dict(name="Air", method=args.model[0], identifier="first_use_air")
    if args.model[0] == "bicubic":
        settings["table_options"] = dict(
            h_min=200e3,
            h_max=850e3,
            p_min=40e3,
            p_max=350e3,
            N_h=32,
            N_p=32,
            N_p_sat=4,
            coarse_step=1,
            table_dir=str(args.table_directory.resolve()),
        )
    start = perf_counter()
    fluid = jp.Fluid(**settings)
    setup = perf_counter() - start
    compression = args.machine[0] == "compression"
    p_in, temp, p_out = (1e5, 300.0, 2e5) if compression else (2e5, 600.0, 1e5)
    steps = 10 if compression else 50
    start = perf_counter()
    inlet = fluid.get_state(jp.PT_INPUTS, p_in, temp)
    jax.block_until_ready(inlet)
    inlet_seconds = perf_counter() - start
    h_in = float(inlet.enthalpy)
    start = perf_counter()
    if args.warmup[0] in ("scalar", "path"):
        # Warm the actual scalar representations used by both source versions.
        # Python/Python at entry, Python/NumPy in the legacy ODE, and consistent
        # NumPy/NumPy at the current host boundary are separate JIT signatures.
        for first, second in (
            (h_in, p_in),
            (h_in, np.float64(p_in)),
            (np.asarray(h_in), np.asarray(p_in)),
        ):
            state = fluid.get_state(jp.HmassP_INPUTS, first, second)
            jax.block_until_ready(state)
        for pressure in (p_out, np.asarray(p_out)):
            isentropic = fluid.get_state(jp.PSmass_INPUTS, pressure, state.entropy)
            jax.block_until_ready(isentropic)
    if args.warmup[0] == "path":
        # This is an explicit native full-state warmup at the reporting shape,
        # outside model construction and numerical iteration. No state adapter.
        path = fluid.get_state(
            jp.HmassP_INPUTS,
            np.full(steps, h_in),
            np.linspace(p_in, p_out, steps),
            **process_state_options(fluid, generalize_quality=not compression),
        )
        jax.block_until_ready(path)
    warmup = perf_counter() - start
    process = (
        components.compression_process if compression else components.expansion_process
    )

    def evaluate():
        return process(
            fluid,
            h_in,
            p_in,
            p_out,
            efficiency=0.85,
            efficiency_type="polytropic",
            num_steps=steps,
        )

    context = jax.log_compiles(True) if args.log_compiles else nullcontext()
    with context:
        start = perf_counter()
        first = evaluate()
        jax.block_until_ready(first)
        first_seconds = perf_counter() - start
    samples = []
    for _ in range(args.samples):
        start = perf_counter()
        last = evaluate()
        jax.block_until_ready(last)
        samples.append(perf_counter() - start)

    for result in (first, last):
        assert result["states"].shape == (steps,)
        np.testing.assert_array_equal(result["states"].status, 0)
        assert np.issubdtype(np.asarray(result["states"].status).dtype, np.integer)
        assert result["states"].identifier == fluid.identifier
        for name in ("temperature", "pressure", "density"):
            assert np.all(np.isfinite(result["states"][name])) and np.all(
                np.asarray(result["states"][name]) > 0
            )
        assert result["specific_work"] > 0
        np.testing.assert_allclose(result["state_out"].pressure, p_out)
    output = plain(first)
    compare(plain(last), output)
    assert (
        source_hash(root) == identity
        and cpu_affinity()["actual_mask"] == affinity["actual_mask"]
    )
    save(
        args.output,
        dict(
            protocol=2,
            source_root=str(root),
            source_sha256=identity,
            affinity=affinity,
            settings=settings,
            model=args.model[0],
            machine=args.machine[0],
            warmup=args.warmup[0],
            steps=steps,
            compilation_logging=args.log_compiles,
            imports_seconds=imports,
            setup_seconds=setup,
            inlet_seconds=inlet_seconds,
            warmup_seconds=warmup,
            first_component_seconds=first_seconds,
            warm_samples_seconds=samples,
            warm_median_seconds=statistics.median(samples),
            environment=dict(
                python=sys.version,
                executable=sys.executable,
                jaxprop_path=jp.__file__,
                persistent_compilation_cache=False,
                python_hash_seed=os.environ.get("PYTHONHASHSEED"),
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
            jaxprop_files={
                str(p.relative_to(Path(jp.__file__).parent)).replace("\\", "/"): sha(p)
                for p in sorted(Path(jp.__file__).parent.rglob("*.py"))
            },
            complete_output=output,
            numerical_checks_passed=True,
        ),
    )


def main():
    args = parser().parse_args()
    if args.samples < 4 or args.cohorts < 1:
        raise ValueError("Use at least four warm samples and one cohort")
    if args.worker_root:
        return worker(args)
    if args.reference_root is None or args.reference_root.resolve() == ROOT:
        raise ValueError("Specify a separate reference source tree")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(__file__, args.output / "harness.py")
    affinity = cpu_affinity(args.cpu_affinity)
    sources = dict(reference=args.reference_root.resolve(), current=ROOT)
    identities = {key: source_hash(root) for key, root in sources.items()}
    report = dict(
        schema=1,
        protocol=2,
        harness_sha256=sha(__file__),
        lock_sha256=sha(ROOT / "poetry.lock"),
        sources={key: str(path) for key, path in sources.items()},
        source_sha256=identities,
        cpu_affinity=affinity,
        samples=args.samples,
        cohorts=args.cohorts,
        compilation_logging=args.log_compiles,
        interpretation="Each invocation runs in its own fresh process. Compare first component calls separately from imports/setup/inlet/explicit warmup; no cold-performance budget is silently assumed.",
        results={},
    )
    destination = args.output / "comparison.json"
    save(destination, report)
    for method in args.model or ("perfect_gas", "bicubic"):
        for machine in args.machine or ("compression", "expansion"):
            for warmup in args.warmup or ("none", "scalar", "path"):
                label = f"{method}/{machine}/{warmup}"
                records = []
                for cohort in range(args.cohorts):
                    observation = {}
                    order = (
                        ("reference", "current")
                        if cohort % 2 == 0
                        else ("current", "reference")
                    )
                    for version in order:
                        stem = f"{method}_{machine}_{warmup}_{cohort}_{version}"
                        output = args.output / f"{stem}.json"
                        command = [
                            sys.executable,
                            str(Path(__file__).resolve()),
                            "--worker-root",
                            str(sources[version]),
                            "--output",
                            str(output),
                            "--samples",
                            str(args.samples),
                            "--model",
                            method,
                            "--machine",
                            machine,
                            "--warmup",
                            warmup,
                            "--table-directory",
                            str(args.table_directory.resolve()),
                        ]
                        if args.cpu_affinity is not None:
                            command += ["--cpu-affinity", hex(args.cpu_affinity)]
                        if args.log_compiles:
                            command.append("--log-compiles")
                        with (args.output / f"{stem}.stdout.log").open(
                            "w", encoding="utf-8"
                        ) as stdout, (args.output / f"{stem}.stderr.log").open(
                            "w", encoding="utf-8"
                        ) as stderr:
                            process = subprocess.run(
                                command,
                                cwd=ROOT,
                                stdout=stdout,
                                stderr=stderr,
                                env=dict(os.environ, PYTHONHASHSEED="0"),
                            )
                        if process.returncode:
                            raise RuntimeError(
                                f"Worker failed: inspect {stem}.stderr.log"
                            )
                        observation[version] = json.loads(
                            output.read_text(encoding="utf-8")
                        )
                        assert observation[version]["numerical_checks_passed"]
                        assert (
                            observation[version]["source_sha256"] == identities[version]
                        )
                    compare(
                        observation["current"]["complete_output"],
                        observation["reference"]["complete_output"],
                    )
                    assert (
                        observation["current"]["environment"]
                        == observation["reference"]["environment"]
                    )
                    assert (
                        observation["current"]["jaxprop_files"]
                        == observation["reference"]["jaxprop_files"]
                    )
                    records.append(observation)
                summary = {
                    version: {
                        name: statistics.median(r[version][name] for r in records)
                        for name in (
                            "setup_seconds",
                            "inlet_seconds",
                            "warmup_seconds",
                            "first_component_seconds",
                            "warm_median_seconds",
                        )
                    }
                    for version in sources
                }
                report["results"][label] = dict(
                    observations=records, medians=summary, complete_outputs_agree=True
                )
                save(destination, report)
                print(
                    label,
                    {
                        k: round(v["first_component_seconds"] * 1000, 2)
                        for k, v in summary.items()
                    },
                    flush=True,
                )
    assert identities == {key: source_hash(root) for key, root in sources.items()}
    assert (
        sha(__file__) == report["harness_sha256"]
        and sha(ROOT / "poetry.lock") == report["lock_sha256"]
    )
    assert cpu_affinity()["actual_mask"] == affinity["actual_mask"]
    report["source_and_protocol_unchanged"] = True
    report["numerical_checks_passed"] = True
    save(destination, report)


if __name__ == "__main__":
    main()
