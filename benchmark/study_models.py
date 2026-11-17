"""Compare complete repeated Rankine study evaluations against a source tree.

Workers use isolated package/module namespaces in separate interpreters. The
ABBA process order and repeated synchronized samples reduce shared-host drift;
all numerical residuals must agree before timing ratios are interpreted.
"""

import argparse
import hashlib
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]


def worker(root, output, samples):
    sys.path[:0] = [
        str(root),
        str(root / "case_studies/rankine_transcritical_vs_trilateral"),
    ]
    import numpy as np
    import thermopt as th
    import run_sensitivity as study
    from fast_evaluator import FastProblem

    assert Path(th.__file__).resolve().parent == root / "thermopt"
    start = perf_counter()
    config = study.configuration("butane", 50, 0.8, 0.8)
    problem = FastProblem(config["problem_formulation"], out_dir=str(output))
    setup = perf_counter() - start
    baseline = problem.fitness(problem.x0).copy()
    for _ in range(3):
        problem.fitness(problem.x0)
    times = []
    for _ in range(samples):
        start = perf_counter()
        actual = problem.fitness(problem.x0)
        times.append(perf_counter() - start)
        np.testing.assert_array_equal(actual, baseline)
    files = [
        *sorted((root / "thermopt").rglob("*.py")),
        root / "case_studies/rankine_transcritical_vs_trilateral/run_sensitivity.py",
        root / "case_studies/rankine_transcritical_vs_trilateral/fast_evaluator.py",
    ]
    record = dict(
        setup_seconds=setup,
        residuals=baseline.tolist(),
        warm_samples_seconds=times,
        warm_median_seconds=statistics.median(times),
        sources={
            str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in files
        },
    )
    (output / "result.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=60)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.samples < 4:
        parser.error("Use at least four samples")
    root, output = args.reference_root.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    os.environ.update(
        JAX_ENABLE_X64="true",
        MPLBACKEND="Agg",
        DISABLE_PLOTS="1",
        OMP_NUM_THREADS="1",
        PYTHONHASHSEED="0",
    )
    if args.worker:
        worker(root, output, args.samples)
        return
    if root == ROOT:
        parser.error("Reference must be an independent source tree")
    results = {"reference": [], "current": []}
    for index, label in enumerate(["reference", "current", "current", "reference"]):
        directory = output / f"{index}_{label}"
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--worker",
            "--reference-root",
            str(root if label == "reference" else ROOT),
            "--output",
            str(directory),
            "--samples",
            str(args.samples),
        ]
        with (output / f"{index}_{label}.log").open("w", encoding="utf-8") as log:
            subprocess.run(
                command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True
            )
        results[label].append(json.loads((directory / "result.json").read_text()))
    import numpy as np

    baseline = results["reference"][0]["residuals"]
    for entries in results.values():
        for entry in entries:
            np.testing.assert_allclose(
                entry["residuals"], baseline, rtol=1e-12, atol=1e-12
            )
    medians = {
        label: statistics.median(
            [value for run in entries for value in run["warm_samples_seconds"]]
        )
        for label, entries in results.items()
    }
    ratio = medians["current"] / medians["reference"]
    record = dict(
        order="reference/current/current/reference",
        samples_per_process=args.samples,
        same_residuals=True,
        medians_seconds=medians,
        ratio=ratio,
        warm_ratio_limit=1.10,
        timing_flag=ratio > 1.10,
        results=results,
        packages={
            name: metadata.version(name)
            for name in ("jaxprop", "CoolProp", "numpy", "scipy")
        },
        harness_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    (output / "comparison.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {"same_residuals": True, "medians_seconds": medians, "ratio": ratio},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
