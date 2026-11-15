"""Small helpers for reproducible design-point sweeps, using only jaxprop."""

import csv
import importlib.metadata
import json
from pathlib import Path

import numpy as np
import jaxprop as props


def solve_case(solve, fluid, mass_flow, p_in, temperature, p_out, **design):
    """Solve a design point through the functional API and flatten SI outputs.

    solve is solve_compressor or solve_turbine. The explicit jaxprop fluid
    provides the total inlet p,T flash; the solver uses outlet enthalpy.
    The returned row includes all eight performance values, outlet conditions,
    efficiency residual, and angular speed converted from rad/s to rpm.
    """
    inlet = fluid.get_state(props.PT_INPUTS, p_in, temperature)
    solution = solve(fluid, inlet, p_out, mass_flow, design)
    result = solution["performance"]
    return {**result,
            "outlet_temperature_K": float(solution["state_out"].T),
            "outlet_enthalpy_J_kg": float(solution["state_out"].h),
            "residual": solution["residual"],
            "rpm": float(result["angular_speed"] * 30 / np.pi)}


def save_results(rows, directory, model):
    """Write sweep rows to CSV with backend/version metadata; return its folder."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / "trends.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "model": model, "fluid": "Air", "backend": "HEOS",
        "iteration_variable": "outlet enthalpy", "rows": len(rows),
        "meaning": "Design-point parameter sweeps; not measured off-design maps or digitized paper data.",
        "versions": {name: importlib.metadata.version(name) for name in
                     ["thermopt", "jaxprop", "numpy", "scipy"]},
    }
    (directory / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return directory
