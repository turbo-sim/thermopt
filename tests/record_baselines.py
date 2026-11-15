"""Record characterization data intentionally, never as part of pytest.

Existing files are protected unless --overwrite is supplied. Review numerical
changes before accepting a new reference. Run from any directory.
"""

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from uuid import uuid4

os.environ["MPLBACKEND"] = "Agg"
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from regression_support import (CASES, DATA, EXAMPLE_CASES, make_problem,
                                nearby_design, plain, snapshot)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    target = DATA / "cycle_baselines.json"
    if target.exists() and not args.overwrite:
        parser.error("Baseline exists. Inspect changes and use --overwrite deliberately.")
    CASES.mkdir(parents=True, exist_ok=True)
    for name, source in EXAMPLE_CASES.items():
        destination = CASES / f"{name}.yaml"
        if not destination.exists():
            if name == "orc_wet_expansion":
                raise FileNotFoundError("Restore the explicit wet-expansion test fixture; it is not the dry example.")
            shutil.copyfile(ROOT / "examples" / source, destination)
    record = {
        "schema": 1,
        "description": "Current HEOS characterization, not a proof of physical correctness or feasibility.",
        "versions": {name: importlib.metadata.version(name) for name in
                     ["thermopt", "jaxprop", "CoolProp", "numpy", "scipy", "pysolver_view"]},
        "thermopt_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "cases": {},
    }
    # Ordinary permissions avoid private-mode temp ACL issues on Windows.
    artifact_parent = ROOT / ".test_artifacts"
    artifact_parent.mkdir(exist_ok=True)
    out_dir = artifact_parent / ("record-" + uuid4().hex)
    out_dir.mkdir()
    try:
        for name in EXAMPLE_CASES:
            try:
                problem = make_problem(name, out_dir)
                initial = problem.x0.copy()
                perturbed = nearby_design(problem)
                item = {"variable_names": problem.variable_names,
                        "topology": problem.cycle_topology,
                        "x0": plain(initial),
                        "nearby_x": plain(perturbed),
                        "initial": snapshot(problem, initial),
                        "nearby": snapshot(problem, perturbed)}
                if name == "ptes_turbo":
                    item["model_corrections"] = [
                        "Compressor relative inlet velocity follows Figure 4: w/U = 0.38 + 3*phi - 7.8*phi**2; corrects the printed Equation 23 error."
                    ]
                record["cases"][name] = item
                print(name, "recorded", flush=True)
            except Exception as error:
                print(name, "FAILED", type(error).__name__, str(error), flush=True)
                record["cases"][name] = {"recording_error": type(error).__name__ + ": " + str(error)}
    finally:
        resolved = out_dir.resolve()
        if resolved.parent != artifact_parent.resolve():
            raise RuntimeError(f"Unexpected baseline artifact path: {resolved}")
        shutil.rmtree(resolved)
    failures = [name for name, item in record["cases"].items() if "recording_error" in item]
    if failures:
        raise RuntimeError(f"No baseline replaced; failed cases: {failures}")
    staging = target.with_suffix(".json.tmp")
    staging.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    staging.replace(target)
    print(target)


if __name__ == "__main__":
    main()
