"""Run complete eager cycle benchmarks and compare matching phase checkpoints.

Workload implementations currently live in development/bicubic/validate_060.py;
the migration plan consolidates them here without duplicating their equations.
No ordinary pytest assertion depends on these shared-machine wall times.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata as metadata
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
WORKLOAD = ROOT / "development/bicubic/validate_060.py"


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def save(path, data):
    Path(path).write_text(
        json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(value):
    """Ignore cache storage location, preserving all physical/solver settings."""
    if isinstance(value, dict):
        return {k: canonical(v) for k, v in value.items() if k != "table_dir"}
    if isinstance(value, list):
        return [canonical(v) for v in value]
    return value


def compare_checkpoint(current, previous):
    """Reject unmatched/invalid references before assessing known timing ratios."""
    if not previous.get("numerical_checks_passed") or (
        previous.get("comparison_signature") != current["comparison_signature"]
    ):
        raise ValueError(
            "incompatible baseline; compare matching protocols/environments/inputs"
        )
    timings = current["timings_seconds"]
    reference = previous["timings_seconds"]
    if not timings or timings.keys() != reference.keys():
        raise ValueError("incompatible baseline workload/model set")
    if not all(
        math.isfinite(value) and value > 0
        for value in [*timings.values(), *reference.values()]
    ):
        raise ValueError("timings must be finite and positive")
    ratios = {key: value / reference[key] for key, value in timings.items()}
    return ratios, [key for key, value in ratios.items() if value > 1.10]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--profile", choices=("smoke", "phase", "release"), default="smoke"
    )
    parser.add_argument("--grid-size", type=int)
    parser.add_argument("--samples", type=int, default=12)
    parser.add_argument(
        "--baseline", type=Path, help="Accepted preceding run's manifest.json"
    )
    parser.add_argument(
        "--output", type=Path, help="New directory; existing runs are never overwritten"
    )
    parser.add_argument(
        "--table-directory", type=Path, default=ROOT / "benchmark/results/tables"
    )
    args = parser.parse_args()
    if args.samples < 4 or args.samples % 2:
        parser.error(
            "Use an even sample count of at least four for balanced backend order"
        )
    grid = args.grid_size or (200 if args.profile == "release" else 100)
    if grid < 2:
        parser.error("Grid size must be at least two")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    output = (
        args.output or ROOT / f"benchmark/results/{stamp}_{args.profile}"
    ).resolve()
    output.mkdir(parents=True, exist_ok=False)
    args.table_directory = args.table_directory.resolve()

    # These benchmark processes have a declared CPU/x64 scientific profile.
    for name, value in {
        "JAX_ENABLE_X64": "true",
        "JAX_PLATFORM_NAME": "cpu",
        "MPLBACKEND": "Agg",
        "OMP_NUM_THREADS": "1",
    }.items():
        os.environ[name] = value
    import jax
    import jaxprop as jp
    import yaml

    if jp.__version__.split(".")[:2] != ["0", "6"]:
        raise RuntimeError("This checkpoint requires the JaxProp 0.6 interface")
    if jp.__version__ != metadata.version("jaxprop"):
        raise RuntimeError(
            "JaxProp import and installed distribution versions disagree"
        )
    if not Path(jp.__file__).resolve().is_relative_to(Path(sys.prefix).resolve()):
        raise RuntimeError(
            "Install JaxProp in this environment; remove checkout/PYTHONPATH shadowing"
        )
    if not jax.config.x64_enabled or jax.default_backend() != "cpu":
        raise RuntimeError("The declared benchmark profile is CPU/x64")

    jobs = [("recuperated_working", "sco2_recuperated", False)]
    if args.profile != "smoke":
        jobs += [
            ("recuperated_all", "sco2_recuperated", True),
            ("simple_polytropic", "sco2_simple_polytropic", False),
            ("recompression", "sco2_recompression", False),
            ("refrigeration", "refrigeration_air", True),
        ]
    environment = dict(
        python=platform.python_version(),
        platform=platform.platform(),
        machine=platform.machine(),
        processor=platform.processor(),
        logical_cpus=os.cpu_count(),
        device=str(jax.devices()[0]),
        x64=jax.config.x64_enabled,
        dependencies={
            n: metadata.version(n)
            for n in (
                "jaxprop",
                "jax",
                "jaxlib",
                "CoolProp",
                "numpy",
                "scipy",
                "pysolver_view",
                "equinox",
                "diffrax",
                "optimistix",
                "jaxtyping",
                "lineax",
                "matplotlib",
                "pandas",
            )
        },
    )
    manifest = dict(
        schema=1,
        started_at_utc=datetime.now(timezone.utc).isoformat(),
        profile=args.profile,
        grid_size=grid,
        samples=args.samples,
        environment=environment,
        jaxprop_path=str(Path(jp.__file__).resolve()),
        thermopt_root=str(ROOT),
        lock_sha256=sha(ROOT / "poetry.lock"),
        runner_sha256=sha(__file__),
        workload_sha256=sha(WORKLOAD),
        thermopt_source_sha256=hashlib.sha256(
            b"".join(p.read_bytes() for p in sorted((ROOT / "thermopt").rglob("*.py")))
        ).hexdigest(),
        warm_ratio_limit=1.10,
        results={},
        timings_seconds={},
        effective_inputs={},
        interpretation="Shared-machine timing flags require fresh repeat cohorts; never silently relax limits.",
    )
    destination = output / "manifest.json"
    save(destination, manifest)
    for label, case, all_streams in jobs:
        directory = output / label
        command = [
            sys.executable,
            str(WORKLOAD),
            "--case",
            case,
            "--grids",
            str(grid),
            "--samples",
            str(args.samples),
            "--output",
            str(directory),
            "--table-directory",
            str(args.table_directory),
        ]
        if all_streams:
            command.append("--all-fluids")
        if args.profile == "release" and case == "sco2_recuperated":
            command.append("--optimize")
        with (output / f"{label}.log").open("w", encoding="utf-8") as stream:
            process = subprocess.run(
                command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT
            )
        manifest["results"][label] = dict(exit_code=process.returncode, command=command)
        save(destination, manifest)
        if process.returncode:
            raise SystemExit(f"Numerical/workflow check failed: see {output/label}.log")
        result = read(directory / f"{case}_summary.json")
        for model, data in result["comparisons"].items():
            for backend, seconds in data["warm_median_seconds"].items():
                manifest["timings_seconds"][f"{label}/{model}/{backend}"] = seconds
        configs = [
            result["source_setup"]["configuration"],
            *[v["setup"]["configuration"] for v in result["comparisons"].values()],
        ]
        manifest["effective_inputs"][label] = [
            canonical(yaml.safe_load((ROOT / p).read_text(encoding="utf-8")))
            for p in configs
        ]
        save(destination, manifest)

    signature = {
        k: manifest[k]
        for k in (
            "schema",
            "profile",
            "grid_size",
            "samples",
            "environment",
            "runner_sha256",
            "workload_sha256",
            "effective_inputs",
        )
    }
    after = hashlib.sha256(
        b"".join(p.read_bytes() for p in sorted((ROOT / "thermopt").rglob("*.py")))
    ).hexdigest()
    if (
        after != manifest["thermopt_source_sha256"]
        or sha(__file__) != manifest["runner_sha256"]
        or sha(WORKLOAD) != manifest["workload_sha256"]
        or sha(ROOT / "poetry.lock") != manifest["lock_sha256"]
    ):
        manifest["comparison_status"] = (
            "source, harness or lock changed during this run"
        )
        save(destination, manifest)
        raise SystemExit(manifest["comparison_status"])
    manifest["source_and_protocol_unchanged_during_run"] = True
    manifest["comparison_signature"] = hashlib.sha256(
        json.dumps(signature, sort_keys=True).encode()
    ).hexdigest()
    manifest["numerical_checks_passed"] = True
    manifest["timing_flags"] = []
    if args.baseline:
        previous = read(args.baseline)
        try:
            ratios, flags = compare_checkpoint(manifest, previous)
        except ValueError as error:
            manifest["comparison_status"] = str(error)
            save(destination, manifest)
            raise SystemExit(manifest["comparison_status"])
        manifest["baseline"] = str(args.baseline.resolve())
        manifest["ratios_to_baseline"] = ratios
        manifest["timing_flags"] = flags
        manifest["comparison_status"] = (
            "repeat timing flags"
            if manifest["timing_flags"]
            else "within warm budgets in this cohort"
        )
    else:
        manifest["comparison_status"] = (
            "baseline candidate; no between-phase performance acceptance claimed"
        )
    manifest["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
    save(destination, manifest)
    print(destination)
    print(manifest["comparison_status"])
    if manifest["timing_flags"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
