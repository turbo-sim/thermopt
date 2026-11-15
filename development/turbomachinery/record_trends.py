"""Explicitly record current radial-machine regression behavior.

This is not a paper-validation dataset. Tests never invoke this recorder.
"""

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from development.turbomachinery import demo_compressor, demo_turbine


def main():
    """Record reviewed sweep references only after an explicit overwrite request."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    target = ROOT / "tests/data/turbomachinery_trends.json"
    if target.exists() and not args.overwrite:
        parser.error("Reference exists; review differences before using --overwrite.")
    module = ROOT / "thermopt/components/turbomachinery_nondimensional.py"
    record = {
        "purpose": "HEOS regression with the Figure 4 compressor velocity correction; not independent experimental validation.",
        "model_corrections": [
            "Compressor relative inlet velocity follows Figure 4: w/U = 0.38 + 3*phi - 7.8*phi**2; corrects the printed Equation 23 error."
        ],
        "module_sha256": hashlib.sha256(module.read_bytes()).hexdigest(),
        "versions": {name: importlib.metadata.version(name) for name in
                     ["thermopt", "jaxprop", "CoolProp", "numpy", "scipy"]},
        "compressor": demo_compressor.run_sweep(points=5),
        "turbine": demo_turbine.run_sweep(points=5),
    }
    staging = target.with_suffix(".json.tmp")
    staging.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    staging.replace(target)
    print(target)


if __name__ == "__main__":
    main()
