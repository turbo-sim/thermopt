"""Reproduce the seven component-model figures from Parisi et al. (2025).

Run from any directory: python development/turbomachinery/paper_figures.py
Published references are checked-in PDF-vector extractions. Curves calculated
here use the same correlation helpers as the production machine models.
"""

import argparse
import csv
import json
import os
from pathlib import Path
import sys

os.environ.setdefault("MPLBACKEND", "Agg")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import jaxprop as props
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.optimize import brentq

from thermopt.components.turbomachinery_nondimensional import (
    compressor_reference_coefficients, compressor_mach_correction,
    compressor_relative_velocity_ratio, compressor_diameter_to_chord_ratio,
    turbine_baseline_efficiency, turbine_kinetic_energy_fraction,
)

HERE = Path(__file__).resolve().parent
REFERENCE = HERE / "paper_reference/curves.json"
TITLES = {
    1: "Compressor baseline efficiency", 2: "Compressor work coefficient",
    3: "Compressor pressure-ratio effect", 4: "Compressor relative inlet velocity",
    5: "Compressor diameter/chord ratio", 6: "Turbine total-to-total efficiency",
    7: "Turbine total-to-static efficiency",
}
Y_LABELS = {
    1: "Polytropic efficiency", 2: "Work coefficient, $I$", 3: "Polytropic efficiency",
    4: "Relative velocity / tip speed, $w/U$", 5: "Diameter / chord, $D/c$",
    6: "Rotor total-to-total efficiency", 7: "Rotor total-to-static efficiency",
}
MODEL_LABELS = {1: "thermopt, Eq. 14", 2: "thermopt, Eqs. 12/17",
                4: "thermopt, corrected Eq. 23 / Fig. 4", 5: "thermopt, Eqs. 24/25",
                6: "thermopt, Eq. 29", 7: "thermopt, Eqs. 29/36/37"}
STYLE = {"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 12,
         "axes.labelsize": 11, "legend.fontsize": 8, "axes.grid": True,
         "grid.alpha": .18, "axes.axisbelow": True, "text.usetex": False,
         "mathtext.fontset": "dejavusans", "xtick.direction": "in", "ytick.direction": "in"}


def read_reference():
    """Load attributed publication vectors; no PDF parser or network is needed."""
    return json.loads(REFERENCE.read_text(encoding="utf-8"))


def air_heat_capacity_ratio():
    """Get the perfect-gas Air cp/cv ratio through jaxprop at 298 K and 1 bar.

    Figure 3 illustrates the Mach correction alone. Its perfect-gas form is
    independent of mass flow and absolute inlet temperature/pressure. The
    state conditions match the compressor study's Table 3 reference values;
    the caption does not fully specify a real-fluid design point.
    """
    fluid = props.FluidPerfectGas("Air", T_ref=298.0)
    state = fluid.get_state(props.PT_INPUTS, 1e5, 298.0)
    return float(state.cp / state.cv)


def mach_only_efficiency(flow_coefficient, pressure_ratio, gamma):
    """Solve Figure 3's baseline-plus-Mach efficiency illustration.

    Eqs. (9), (11), (12), (14), (18)-(20) give, for a perfect gas,
    T_out/T_in = PR**((gamma-1)/(gamma*eta_p)) and
    Ma_u**2 = (T_out/T_in-1)/((gamma-1)*I). Production helpers supply I,
    baseline eta_p and the Mach penalty. Reynolds, clearance and backsweep
    corrections are held at their reference values to isolate Mach effects.

    Return the highest-efficiency root, continuously connected to the
    baseline. The loss correlation can lose its physical root at large flow
    coefficient/high PR; return NaN when no root is found in the searched
    efficiency interval 0.2-1, rather than inventing or clipping a solution.
    Such points are outside the visible published curve in this comparison.
    """
    work, head = compressor_reference_coefficients(flow_coefficient)
    baseline = head / work

    def residual(efficiency):
        """Match trial polytropic efficiency to baseline plus tip-Mach penalty."""
        temperature_ratio = pressure_ratio**((gamma - 1) / (gamma * efficiency))
        tip_mach = np.sqrt((temperature_ratio - 1) / ((gamma - 1) * work))
        return baseline + compressor_mach_correction(flow_coefficient, tip_mach) - efficiency

    candidates = np.linspace(.999999, .2, 80)
    values = residual(candidates)
    changes = np.flatnonzero(values[:-1] * values[1:] <= 0)
    if not len(changes):
        return np.nan
    index = changes[0]
    return brentq(residual, candidates[index + 1], candidates[index], xtol=1e-12)


def evaluate_curve(figure, coordinates, key="proposed", *, gamma=None):
    """Evaluate one production correlation at dimensionless plot coordinates.

    Figures 1/2 use Table 4's reference coefficients; Figure 3 isolates Mach
    loss as documented above; Figures 4/5 use the production Reynolds-geometry
    submodels; Figures 6/7 use rotor baseline efficiency and exit kinetic
    fraction before scalloping, clearance or Reynolds penalties.
    """
    coordinates = np.asarray(coordinates, dtype=float)
    if figure in (1, 2):
        work, head = compressor_reference_coefficients(coordinates)
        return head / work if figure == 1 else work
    if figure == 3:
        pressure_ratio = float(key.removeprefix("pr_"))
        if gamma is None:
            gamma = air_heat_capacity_ratio()
        return np.array([mach_only_efficiency(value, pressure_ratio, gamma) for value in coordinates])
    if figure == 4:
        return compressor_relative_velocity_ratio(coordinates)
    if figure == 5:
        return compressor_diameter_to_chord_ratio(coordinates)
    if figure == 6:
        return turbine_baseline_efficiency(coordinates)
    if figure == 7:
        return turbine_baseline_efficiency(coordinates) / (1 + turbine_kinetic_energy_fraction(coordinates))
    raise ValueError(f"Unknown figure: {figure}")


def comparison_metrics(figure, source, specification, gamma):
    """Compare calculated values with independent PDF polyline vertices.

    Exclude points touching clipping boundaries: the PDF renderer can create
    linearly interpolated boundary intersections that are not evaluations of
    the underlying nonlinear curve. Record exclusions and failed evaluations.
    Tolerance 2e-4 is in the plotted dimensionless quantity, not percent.
    Figure 3 uses 5e-4 (0.05 efficiency percentage points), separately labeled
    as an approximate reconstruction because the caption omits full state
    conditions and this script uses a specified perfect-gas approximation.
    This is reproduction of publication graphics, not validation against raw
    experiments. Mean-line scatter is displayed but is not scored here.
    """
    x, expected = np.asarray(source["x"]), np.asarray(source["y"])
    actual = evaluate_curve(figure, x, source["key"], gamma=gamma)
    x_min, x_max, y_min, y_max = specification["axis_limits"]
    margin = (y_max - y_min) * 2e-4
    interior = (expected > y_min + margin) & (expected < y_max - margin)
    finite = np.isfinite(actual)
    selected = interior & finite
    difference = actual[selected] - expected[selected]
    maximum = float(np.max(np.abs(difference)))
    tolerance = 5e-4 if figure == 3 else 2e-4
    agrees = maximum <= tolerance and np.all(finite[interior])
    status = "close_under_stated_assumptions" if figure == 3 else "agrees_with_published_curve"
    return dict(figure=figure, series=source["key"], points_compared=int(selected.sum()),
                clipped_points_excluded=int((~interior).sum()),
                failed_model_points=int((interior & ~finite).sum()),
                rmse=float(np.sqrt(np.mean(difference**2))), max_abs_error=maximum,
                tolerance=tolerance, status=status if agrees else "discrepancy")


def generate_data(points=201):
    """Calculate all model curves and source comparisons without writing files.

    The requested grid controls plot smoothness only. Comparison statistics
    always evaluate the original PDF vertices, independent of that grid.
    """
    if points < 3:
        raise ValueError("points must be at least 3")
    reference = read_reference()
    gamma = air_heat_capacity_ratio()
    curves, comparisons = {}, []
    for number, specification in reference["figures"].items():
        figure = int(number)
        lower, upper = (.01, .15) if figure <= 5 else (.3, .95)
        if figure == 3:
            lower = .005  # Continue through the paper's clipped low-flow boundary.
        x = np.linspace(lower, upper, points)
        curves[number] = []
        for series in specification["series"]:
            if series["kind"] != "published_correlation" or not series["compared"]:
                continue
            y = evaluate_curve(figure, x, series["key"], gamma=gamma)
            curves[number].append(dict(key=series["key"], label=series["label"], x=x, y=y))
            comparisons.append(comparison_metrics(figure, series, specification, gamma))
    velocity = next(item for item in reference["figures"]["4"]["series"] if item.get("compared"))
    coefficients = np.polyfit(velocity["x"], velocity["y"], 2)
    return dict(reference=reference, curves=curves, comparisons=comparisons,
                gamma=gamma, figure4_published_quadratic=coefficients.tolist())


def plot_figure(figure, data, ax=None):
    """Draw a paper-numbered comparison, clearly separating computed and copied data.

    Reference polylines are gray and mean-line markers retain their published
    design groups. They are digitized graphics, not regenerated simulations.
    Turbine baseline curves are dashed outside the stated 0.35-0.85 range.
    Figure 4 uses the corrected correlation supported by the published line
    and mean-line markers; the error in printed Equation 23 is documented
    in the report rather than retained in the production calculation.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(7.6, 5.3))
    specification = data["reference"]["figures"][str(figure)]
    curves = data["curves"][str(figure)]
    colors = ["#4f127b", "#b5367a", "#d6456c", "#ee7733"]
    markers = {"theta_0": "#4f127b", "theta_20": "#b5367a", "theta_40": "#ee7733", "meanline": "#ee7733"}
    reference_handles = []
    for source in specification["series"]:
        if source["kind"] == "published_meanline":
            label = source["label"] + " (PDF)"
            if figure in (1, 2, 5):
                treatment = "baseline" if source["key"] == "theta_40" else "corrected"
                label = source["label"] + f", {treatment} (PDF)"
            ax.plot(source["x"], source["y"], linestyle="none", marker=source["marker"],
                    markersize=4, markerfacecolor="none", markeredgewidth=.7,
                    color=markers[source["key"]], label=label, zorder=4)
        elif figure == 3:
            index = next(i for i, curve in enumerate(curves) if curve["key"] == source["key"])
            ax.plot(source["x"], source["y"], color=colors[index], linestyle="none",
                    marker="o", markersize=2.5, markerfacecolor="none", markeredgewidth=.65, alpha=.65)
        else:
            style = "-" if source["compared"] or source["key"] == "kofskey" else ":"
            if source["key"] == "casey":
                style = "--"
            color = "#777777" if source["compared"] else "#333333"
            width = 3 if source["compared"] else 1.2
            ax.plot(source["x"], source["y"], color=color, linestyle=style, linewidth=width,
                    alpha=.45 if source["compared"] else .8, label=source["label"] + " (PDF)")
    for index, curve in enumerate(curves):
        x, y = curve["x"], curve["y"]
        if figure == 3:
            ax.plot(x, y, color=colors[index], linewidth=1.6, label=curve["label"] + " (model)")
        elif figure >= 6:
            valid = (x >= .35) & (x <= .85)
            ax.plot(x[valid], y[valid], color="#0072b2", linewidth=1.6, label=MODEL_LABELS[figure])
            for outside in (x < .35, x > .85):
                ax.plot(x[outside], y[outside], color="#0072b2", linewidth=1.3, linestyle="--")
        else:
            ax.plot(x, y, color="#0072b2", linewidth=1.6, label=MODEL_LABELS[figure])
    if figure == 3:
        reference_handles.append(Line2D([], [], color="#555555", marker="o", markersize=3,
                                        markerfacecolor="none", linestyle="none", label="PDF vector samples"))
    if figure >= 6:
        reference_handles.append(Line2D([], [], color="#0072b2", linestyle="--", label="Outside stated validity range"))
    x_min, x_max, y_min, y_max = specification["axis_limits"]
    ax.set(xlim=(x_min, x_max), ylim=(y_min, y_max),
           xlabel="Flow coefficient, $\\phi$" if figure <= 5 else "Rotor-exit specific speed, $\\omega_s$",
           ylabel=Y_LABELS[figure], title=f"Figure {figure} | {TITLES[figure]}")
    ax.set_xticks([0, .05, .10, .15] if figure <= 5 else np.arange(.3, .91, .1))
    y_ticks = {1: [.6, .7, .8, .9], 2: [.5, .6, .7, .8],
               3: [.7, .75, .8, .85, .9], 4: [.4, .5, .6, .7],
               5: [1, 2, 3, 4, 5], 6: [.8, .85, .9, .95], 7: [.8, .85, .9, .95]}
    ax.set_yticks(y_ticks[figure])
    ax.tick_params(top=True, right=True)
    handles, labels = ax.get_legend_handles_labels()
    handles += reference_handles
    labels += [handle.get_label() for handle in reference_handles]
    ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False)
    return ax.figure


def write_report(data, directory):
    """Save machine-readable comparison metrics and a human-readable findings report."""
    report = dict(source=data["reference"]["source"], gamma=data["gamma"],
                  figure3_method="Perfect-gas baseline plus Mach loss only; highest-efficiency root; no Reynolds/clearance/scallop corrections.",
                  figure4_published_quadratic=data["figure4_published_quadratic"],
                  comparisons=data["comparisons"])
    (directory / "comparison.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    lines = ["# Paper component-figure reproduction", "",
             "[Parisi et al. (2025)](https://doi.org/10.1016/j.enconman.2025.120085), Figures 1-7. Published graphics: CC BY 4.0. Reference curves/markers were extracted from PDF vectors, not regenerated from the original mean-line datasets.", "",
             "| Figure | Curve | Points | RMSE | Maximum absolute error | Outcome |",
             "|---|---|---:|---:|---:|---|"]
    for row in data["comparisons"]:
        lines.append(f"| {row['figure']} | {row['series']} | {row['points_compared']} | {row['rmse']:.6g} | {row['max_abs_error']:.6g} | {row['status']} |")
    lines += ["", "Errors use the dimensionless plotted quantity; efficiency errors are fractions. Mean-line scatter is shown for context and is not treated as an exact correlation target. Clipped polyline intersections are excluded from statistics; counts are in the JSON file.", "",
              "## Figure 4 and corrected Equation 23", "",
              "Thermopt uses $w/U = 0.38 + 3\\phi - 7.8\\phi^2$, matching the published Figure 4 curve and its supporting mean-line markers. Printed Equation 23 incorrectly gives $w/U = 0.42 + 5\\phi - 14\\phi^2$; the original implementation followed that error. The production docstring records this correction. Publication-vector reference coordinates remain unchanged.", "",
              "## Interpretation", "",
              f"Figure 3 isolates Mach loss, using jaxprop perfect-gas Air with $\\gamma={data['gamma']:.8f}$. Its caption does not specify a complete real-fluid machine design. This reconstruction uses the paper's equations with other corrections at reference values; it is not a complete HEOS design-point sweep. Its separate reproduction tolerance is 0.05 efficiency percentage points, while other comparisons use 0.0002 in their plotted quantity. All raw errors remain visible above.", "",
              "Figures 6 and 7 show rotor baseline efficiencies, before scalloping/clearance/Reynolds penalties. Figure 7 uses $\\eta_{{ts}}=\\eta_{{tt}}/(1+\\kappa)$, not the default diffuser loss fraction 0.35. The turbine coefficient remains 0.5, matching the published Figure 6. Dashed portions extend beyond the stated specific-speed range 0.35-0.85 solely to compare with the original figure.", "",
              "The Figure 1 legend lists a 1.4 percentage-point shift for the 20-degree mean-line designs, whereas Equation 16 gives 1.2 points. Published marker positions are retained as drawn; no second correction or silent adjustment is applied.", "",
              "All figures: calculated model curves plus digitized publication references. Matching a published curve is a reproduction check, not independent experimental validation.", ""]
    (directory / "comparison.md").write_text("\n".join(lines), encoding="utf-8")


def save_outputs(data, directory):
    """Write seven PNG/SVG pairs, an overview, evaluated CSV and comparison reports."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(STYLE):
        for number in range(1, 8):
            figure = plot_figure(number, data)
            figure.text(.1, .035, "Calculated: thermopt | References: published PDF vectors, Parisi et al. (2025), CC BY 4.0", fontsize=7)
            figure.subplots_adjust(left=.15, right=.64, bottom=.17, top=.9)
            stem = f"figure_{number:02d}"
            for extension in ("png", "svg"):
                figure.savefig(directory / f"{stem}.{extension}", dpi=180, bbox_inches="tight", pad_inches=.15)
            plt.close(figure)
        figure, axes = plt.subplots(4, 2, figsize=(15, 18))
        for number, ax in zip(range(1, 8), axes.flat):
            plot_figure(number, data, ax)
            ax.get_legend().remove()
        axes.flat[-1].axis("off")
        axes.flat[-1].text(.02, .85, "Blue: thermopt correlation\nGray/black: published reference curves\nMarkers: digitized mean-line designs\n\nFigure 4 uses the corrected correlation.\nFigure 3 colors distinguish pressure ratios.\n\nParisi et al. (2025), CC BY 4.0\ndoi:10.1016/j.enconman.2025.120085", va="top", fontsize=12, linespacing=1.8)
        figure.tight_layout(pad=2.5)
        figure.savefig(directory / "overview.png", dpi=150)
        plt.close(figure)
    with (directory / "calculated_curves.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["figure", "series", "coordinate", "value", "status"])
        writer.writeheader()
        for number, curves in data["curves"].items():
            for curve in curves:
                for x, y in zip(curve["x"], curve["y"]):
                    writer.writerow(dict(figure=number, series=curve["key"], coordinate=x,
                                         value=float(y) if np.isfinite(y) else "",
                                         status="solved" if np.isfinite(y) else "no_root_in_0.2_to_1"))
    write_report(data, directory)
    return directory


def main():
    """Parse plot resolution/output options and reproduce all seven component figures."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points", type=int, default=201)
    parser.add_argument("--output", type=Path, default=HERE / "results/paper_figures")
    args = parser.parse_args()
    if args.points < 3:
        parser.error("--points must be at least 3")
    data = generate_data(args.points)
    directory = save_outputs(data, args.output)
    print(f"Wrote seven figure reproductions, CSV and comparison reports to {directory}")
    for row in data["comparisons"]:
        print(f"Figure {row['figure']} {row['series']}: max error={row['max_abs_error']:.6g}, {row['status']}")


if __name__ == "__main__":
    main()
