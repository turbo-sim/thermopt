"""Plot radial-turbine design-point efficiency and speed trends.

Run: .venv/Scripts/python development/turbomachinery/demo_turbine.py
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

from thermopt.components.turbomachinery_nondimensional import solve_turbine
from development.turbomachinery.demo_utils import save_results, solve_case


def run_sweep(points=31):
    """Return design-point turbine trends across geometry and scale choices.

    Each family samples rotor-exit specific speed; rows contain SI performance
    values and the solved efficiency residual. No files are written here.
    """
    fluid = props.Fluid("Air", backend="HEOS")
    families = {
        "mass flow": [(f"m={mass:g} kg/s", mass, {}) for mass in [.1, 2., 20.]],
        "design choices": [
            ("scalloped, diffuser loss 0.35", 2., {}),
            ("unscalloped, diffuser loss 0.35", 2., {"is_scalloped": False}),
            ("scalloped, diffuser loss 0.50", 2., {"diffuser_loss_fraction": .5}),
            ("scalloped, clearance 4%", 2., {"clearance_height_ratio": .04}),
        ],
    }
    rows = []
    for family, cases in families.items():
        for label, mass, design in cases:
            for speed in np.linspace(.35, .85, points):
                result = solve_case(solve_turbine, fluid, mass, 2e5, 600., 1e5,
                                    specific_speed=float(speed), **design)
                rows.append({"family": family, "case": label, "specific_speed": float(speed),
                             "mass_flow_kg_s": mass, **result})
    return rows


def main():
    """Parse sweep options and save CSV, metadata and efficiency figures."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points", type=int, default=31)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "results/turbine")
    args = parser.parse_args()
    if args.points < 3:
        parser.error("--points must be at least 3")
    rows = run_sweep(args.points)
    out = save_results(rows, args.output, "solve_turbine")
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for ax, family in zip(axes[:2], ["mass flow", "design choices"]):
        selected = [row for row in rows if row["family"] == family]
        for label in dict.fromkeys(row["case"] for row in selected):
            curve = [row for row in selected if row["case"] == label]
            ax.plot([row["specific_speed"] for row in curve],
                    [100 * row["isentropic_efficiency"] for row in curve], label=label)
        ax.set(xlabel="Rotor-exit specific speed", ylabel="Overall isentropic efficiency [%]", title=family.capitalize())
        ax.grid(alpha=.25)
        ax.legend(fontsize=8)
    selected = [row for row in rows if row["family"] == "mass flow"]
    for label in dict.fromkeys(row["case"] for row in selected):
        curve = [row for row in selected if row["case"] == label]
        axes[2].plot([row["rpm"] / 1000 for row in curve],
                     [100 * row["isentropic_efficiency"] for row in curve], label=label)
    axes[2].set(xlabel="Rotational speed [krpm]", ylabel="Overall isentropic efficiency [%]", title="Speed and efficiency")
    axes[2].grid(alpha=.25)
    axes[2].legend(fontsize=8)
    fig.suptitle("Radial turbine | Air, inlet 2 bar / 600 K, outlet 1 bar | current implementation")
    fig.tight_layout()
    for suffix in ["png", "svg"]:
        fig.savefig(out / f"efficiency_trends.{suffix}", dpi=160)
    plt.close(fig)
    print(f"Wrote {len(rows)} evaluated designs to {out}")


if __name__ == "__main__":
    main()
