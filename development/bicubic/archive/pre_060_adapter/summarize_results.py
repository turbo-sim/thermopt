"""Collect the backend reruns and plot their recorded convergence histories."""

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt


RESULT_DIR = Path(__file__).resolve().parent / "results"


def summarize():
    summaries = {
        path.parent.name: json.loads(path.read_text(encoding="utf-8"))
        for path in sorted(RESULT_DIR.glob("*/summary.json"))
    }
    reference = summaries["coolprop"]["final"]["energy_analysis"]["cycle_efficiency"]
    rows = []
    figure, axes = plt.subplots(1, 2, figsize=(10, 4.2), layout="constrained")
    for name, summary in summaries.items():
        final = summary["final"]
        energy = final["energy_analysis"]
        reference_check = summary.get("cross_evaluations", {}).get(
            "bicubic_optimum_with_coolprop", final
        )
        row = {
            "run": name,
            "success": summary["success"],
            "efficiency_percent": 100 * energy["cycle_efficiency"],
            "difference_percentage_points": 100 * (energy["cycle_efficiency"] - reference),
            "iterations": summary["iterations"],
            "model_calls": summary["model_calls"],
            "compressor_inlet_MPa": final["states"]["compressor"]["state_in"]["p"] / 1e6,
            "compressor_inlet_K": final["states"]["compressor"]["state_in"]["T"],
            "expander_inlet_MPa": final["states"]["expander"]["state_in"]["p"] / 1e6,
            "expander_inlet_K": final["states"]["expander"]["state_in"]["T"],
            "recuperator_effectiveness": final["design_variables"]["recuperator_effectiveness"],
            "working_fluid_mass_flow_kg_s": energy["mass_flow_working_fluid"],
            "heater_heat_flow_MW": energy["heater_heat_flow"] / 1e6,
            "normalized_violation": final["maximum_normalized_constraint_violation"],
            "coolprop_recheck_efficiency_percent": 100 * reference_check["energy_analysis"]["cycle_efficiency"],
            "coolprop_recheck_cooler_pinch_K": reference_check["minimum_temperature_difference_K"]["cooler"],
        }
        rows.append(row)
        history_path = RESULT_DIR / name / "history.json"
        history = json.loads(history_path.read_text(encoding="utf-8"))
        label = "CoolProp HEOS" if name == "coolprop" else name.replace("_", " ")
        axes[0].plot(
            history["func_count_total"],
            [-100 * value for value in history["objective_value"]],
            label=label,
        )
        axes[1].plot(
            history["func_count_total"], history["constraint_violation"], label=label
        )
    axes[0].set_ylabel("Cycle efficiency [%]")
    axes[1].set_ylabel("Maximum normalized constraint violation")
    axes[1].set_yscale("symlog", linthresh=1e-9)
    axes[1].set_ylim(bottom=0)
    for axis in axes:
        axis.set_xlabel("Model evaluations")
        axis.grid(alpha=0.2)
        axis.legend(frameon=False)
    figure.savefig(RESULT_DIR / "comparison.png", dpi=200)
    plt.close(figure)
    with (RESULT_DIR / "comparison.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (RESULT_DIR / "comparison.json").write_text(
        json.dumps(rows, indent=2), encoding="utf-8"
    )
    for row in rows:
        print(
            f"{row['run']}: {row['efficiency_percent']:.9f}% "
            f"({row['difference_percentage_points']:+.9f} percentage points), "
            f"{row['iterations']} iterations"
        )


if __name__ == "__main__":
    summarize()
