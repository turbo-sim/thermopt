"""Plot centrifugal-compressor design-point efficiency trends.

Run: .venv/Scripts/python development/turbomachinery/demo_compressor.py
"""

import argparse
import os
from pathlib import Path
import sys

os.environ.setdefault("MPLBACKEND", "Agg")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import jaxprop as props
import matplotlib.pyplot as plt
import numpy as np

from thermopt.components.turbomachinery_nondimensional import solve_compressor
from development.turbomachinery.demo_utils import save_results, solve_case


def run_sweep(points=31):
    """Return design-point compressor trends across geometry and scale choices.

    Each family samples flow coefficient; rows contain SI performance
    values and the solved efficiency residual. No files are written here.
    """
    fluid = props.Fluid("Air", backend="HEOS")
    families = {
        "pressure ratio": [(f"PR={pr:g}", 2., pr, {}) for pr in [2., 4., 6.]],
        "mass flow": [(f"m={mass:g} kg/s", mass, 2., {}) for mass in [.1, 2., 20.]],
        "design choices": [
            ("vaned, backsweep 40 deg", 2., 2., {}),
            ("vaneless, backsweep 40 deg", 2., 2., {"has_vaned": False}),
            ("vaned, backsweep 0 deg", 2., 2., {"backsweep": 0.}),
            ("vaned, clearance 4%", 2., 2., {"clearance_height_ratio": .04}),
        ],
    }
    rows = []
    for family, cases in families.items():
        for label, mass, pressure_ratio, design in cases:
            for phi in np.linspace(.02, .14, points):
                result = solve_case(solve_compressor, fluid, mass, 1e5, 300.,
                                    1e5 * pressure_ratio, flow_coefficient=float(phi), **design)
                rows.append({"family": family, "case": label, "flow_coefficient": float(phi),
                             "mass_flow_kg_s": mass, "pressure_ratio": pressure_ratio, **result})
    return rows


def main():
    """Parse sweep options and save CSV, metadata and efficiency figures."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points", type=int, default=31)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "results/compressor")
    args = parser.parse_args()
    if args.points < 3:
        parser.error("--points must be at least 3")
    rows = run_sweep(args.points)
    out = save_results(rows, args.output, "solve_compressor")
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=True)
    for ax, family in zip(axes, dict.fromkeys(row["family"] for row in rows)):
        selected = [row for row in rows if row["family"] == family]
        for label in dict.fromkeys(row["case"] for row in selected):
            curve = [row for row in selected if row["case"] == label]
            ax.plot([row["flow_coefficient"] for row in curve],
                    [100 * row["polytropic_efficiency"] for row in curve], label=label)
        ax.set(xlabel="Flow coefficient", title=family.capitalize())
        ax.grid(alpha=.25)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("Polytropic efficiency [%]")
    fig.suptitle("Centrifugal compressor | Air, inlet 1 bar / 300 K | current implementation")
    fig.tight_layout()
    for suffix in ["png", "svg"]:
        fig.savefig(out / f"efficiency_trends.{suffix}", dpi=160)
    plt.close(fig)
    print(f"Wrote {len(rows)} evaluated designs to {out}")


if __name__ == "__main__":
    main()
