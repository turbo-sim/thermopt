"""Explicitly record full solves; tests never refresh references themselves."""

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import sys
from time import perf_counter
from uuid import uuid4

os.environ["MPLBACKEND"] = "Agg"
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np

from regression_support import DATA, OPTIMIZATION_CASES, make_optimization, plain, snapshot


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    destination = DATA / "optimization_baselines.json"
    if destination.exists() and not args.overwrite:
        parser.error("Reference exists; review it before using --overwrite.")
    parent = ROOT / ".test_artifacts"
    parent.mkdir(exist_ok=True)
    directory = parent / ("solve-" + uuid4().hex)
    directory.mkdir()
    records = {
        "schema": 1,
        "versions": {name: importlib.metadata.version(name) for name in
                     ["thermopt", "jaxprop", "CoolProp", "scipy", "pysolver_view"]},
        "cases": {},
    }
    try:
        for name in OPTIMIZATION_CASES:
            start = perf_counter()
            cycle = make_optimization(name, directory / name)
            initial_objective = float(cycle.problem.f)
            cycle.run_optimization()
            solver = cycle.solver
            if not solver.success:
                raise RuntimeError(f"Cannot record an unsuccessful solve for {name}: {solver.message}")
            point = np.asarray(solver.x_final, dtype=float)
            final = snapshot(cycle.problem, point)
            item = {
                "success": bool(solver.success), "message": str(solver.message),
                "iterations": int(solver.iteration_count),
                "model_calls": int(solver.func_count_tot),
                "recording_seconds": perf_counter() - start,
                "initial_objective": initial_objective,
                "x_final": plain(point), "final": final,
            }
            records["cases"][name] = item
            # Preserve each completed solve if a later model fails.
            destination.write_text(json.dumps(records, indent=2, allow_nan=False) + "\n", encoding="utf-8")
            print(name, solver.success, solver.message, "iterations", solver.iteration_count,
                  "f", float(cycle.problem.f), "ineq", max([0., *cycle.problem.c_ineq]),
                  "seconds", item["recording_seconds"], flush=True)
    finally:
        resolved = directory.resolve()
        if resolved.parent != parent.resolve():
            raise RuntimeError(f"Unexpected solve artifact path: {resolved}")
        shutil.rmtree(resolved)


if __name__ == "__main__":
    main()
