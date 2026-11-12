"""Compare every rerun case and physical result with the saved v2 reference."""

import argparse
import csv
import json
import math
from pathlib import Path

from paths import RESULTS, STUDY


FAMILIES = ("butane", "pentane", "cyclopentane", "novec")
BASELINE_EFFICIENCY = 0.886
PARTIAL_EFFICIENCIES = (0.90, 0.80, 0.70, 0.60)

METRIC_TOLERANCES = {
    "system_efficiency_pct": (1e-4, 0.0),
    "cycle_efficiency_pct": (1e-4, 0.0),
    "fixed_reference_efficiency_pct": (1e-4, 0.0),
    "specific_system_power_kW_per_kg_s": (1e-4, 0.0),
    "net_cycle_power_MW": (1e-4, 0.0),
    "net_system_power_MW": (1e-4, 0.0),
    "power_at_194p44_kg_s_kW": (1e-4, 0.0),
    "source_mass_flow_kg_s": (1e-4, 0.0),
    "source_exit_C": (0.0, 0.02),
    "cooling_exit_C": (0.0, 0.02),
    "inlet_temperature_C": (0.0, 0.02),
    "cooling_rise_K": (0.0, 0.02),
    "inlet_superheat_K": (0.0, 0.02),
    "inlet_quality": (0.0, 1e-3),
    "outlet_quality": (0.0, 1e-3),
    "turbine_inlet_pressure_bar": (1e-3, 0.0),
    "condenser_inlet_pressure_bar": (1e-3, 0.0),
}


def point_name(mode, minimum, efficiency):
    efficiency_name = f"{100 * efficiency:g}".replace(".", "p")
    return f"{mode}_Tmin{minimum:03d}_eta{efficiency_name}"


def expected_cases(family, comparison=False):
    cases = {}
    temperatures = (70,) if comparison else range(140, 19, -5)
    for minimum in temperatures:
        designs = [("superheated", BASELINE_EFFICIENCY)]
        if comparison:
            designs.append(("partial_evaporation", BASELINE_EFFICIENCY))
        else:
            designs.extend(("partial_evaporation", efficiency)
                           for efficiency in PARTIAL_EFFICIENCIES)
            if minimum == 70:
                designs.append(("partial_evaporation", BASELINE_EFFICIENCY))
        for mode, efficiency in designs:
            name = f"{family}_{mode}" if comparison else point_name(mode, minimum, efficiency)
            cases[name] = dict(family=family, mode=mode, source_min_C=minimum,
                               turbine_efficiency=efficiency, success=minimum < 140)
    return cases


def load_table(path, failures):
    rows = {}
    try:
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            if not reader.fieldnames or "case" not in reader.fieldnames:
                failures.append(dict(file=str(path), error="Missing case column"))
                return rows
            for row in reader:
                case = row["case"]
                if not case or case in rows:
                    failures.append(dict(file=str(path), case=case,
                                         error="Missing or duplicate case name"))
                rows[case] = row
    except (OSError, csv.Error) as error:
        failures.append(dict(file=str(path), error=str(error)))
    return rows


def check_grid(dataset, role, rows, expected, failures):
    missing = sorted(set(expected) - set(rows))
    unexpected = sorted(set(rows) - set(expected))
    if missing or unexpected:
        failures.append(dict(dataset=dataset, role=role, error="Case grid differs",
                             missing_cases=missing, unexpected_cases=unexpected))
    for case in sorted(set(rows) & set(expected)):
        row, target = rows[case], expected[case]
        for field in ("family", "mode", "source_min_C", "turbine_efficiency", "success"):
            value = row.get(field, "")
            wanted = target[field]
            if field == "success":
                matches = value.lower() == str(wanted).lower()
            elif field in ("source_min_C", "turbine_efficiency"):
                try:
                    matches = math.isclose(float(value), wanted, rel_tol=0.0, abs_tol=1e-12)
                except (TypeError, ValueError):
                    matches = False
            else:
                matches = value == wanted
            if not matches:
                failures.append(dict(dataset=dataset, role=role, case=case,
                                     field=field, expected=wanted, actual=value,
                                     error="Case metadata differs"))
        if not target["success"]:
            if row.get("status") != "no_geothermal_temperature_drop" or row.get("attempts") != "0":
                failures.append(dict(dataset=dataset, role=role, case=case,
                                     error="Unavailable endpoint status or attempt count differs"))
            populated = [metric for metric in METRIC_TOLERANCES
                         if row.get(metric, "").strip().lower() not in ("", "nan")]
            if populated:
                failures.append(dict(dataset=dataset, role=role, case=case,
                                     metrics=populated, error="Unavailable endpoint has physical results"))


def compare_table(dataset, relative_path, expected, reference_dir, results_dir, failures):
    reference = load_table(reference_dir / relative_path, failures)
    current = load_table(results_dir / relative_path, failures)
    check_grid(dataset, "reference", reference, expected, failures)
    check_grid(dataset, "rerun", current, expected, failures)
    finite_cases = sorted(case for case, row in expected.items() if row["success"])
    summaries = []
    for metric, (rtol, atol) in METRIC_TOLERANCES.items():
        differences = []
        failed_cases = []
        for case in finite_cases:
            if case not in reference or case not in current:
                failed_cases.append(case)
                continue
            try:
                old = float(reference[case][metric])
                new = float(current[case][metric])
            except (KeyError, TypeError, ValueError):
                failures.append(dict(dataset=dataset, case=case, metric=metric,
                                     error="Missing or nonnumeric physical result"))
                failed_cases.append(case)
                continue
            if not math.isfinite(old) or not math.isfinite(new):
                failures.append(dict(dataset=dataset, case=case, metric=metric,
                                     error="Nonfinite physical result"))
                failed_cases.append(case)
                continue
            absolute = abs(new - old)
            relative = absolute / abs(old) if old else (0.0 if absolute == 0 else math.inf)
            allowed = atol + rtol * abs(old)
            differences.append(dict(case=case, absolute=absolute, relative=relative))
            if absolute > allowed:
                failed_cases.append(case)
                failures.append(dict(dataset=dataset, case=case, metric=metric,
                                     reference=old, rerun=new, absolute_difference=absolute,
                                     allowed_absolute_difference=allowed,
                                     error="Physical result exceeds tolerance"))
        worst_absolute = max(differences, key=lambda item: item["absolute"], default=None)
        worst_relative = max(differences, key=lambda item: item["relative"], default=None)
        max_relative = worst_relative["relative"] if worst_relative else None
        if max_relative == math.inf:
            max_relative = "infinity"
        summaries.append(dict(
            dataset=dataset, metric=metric, rtol=rtol, atol=atol,
            expected_cases=len(finite_cases), compared_cases=len(differences),
            passed=not failed_cases, failed_cases=len(failed_cases),
            maximum_absolute_difference=worst_absolute["absolute"] if worst_absolute else None,
            worst_absolute_case=worst_absolute["case"] if worst_absolute else "",
            maximum_relative_difference=max_relative,
            worst_relative_case=worst_relative["case"] if worst_relative else "",
        ))
    return summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, default=STUDY / "reference" / "v2",
                        help="Saved v2 tables, default %(default)s")
    parser.add_argument("--results-dir", type=Path, default=RESULTS,
                        help="Rerun tables and comparison report destination, default %(default)s")
    args = parser.parse_args()
    failures, summaries = [], []
    for family in FAMILIES:
        summaries.extend(compare_table(
            f"sensitivity_{family}", Path(f"sensitivity_{family}") / "sweep_results.csv",
            expected_cases(family), args.reference_dir, args.results_dir, failures))
    comparison_cases = {}
    for family in FAMILIES:
        comparison_cases.update(expected_cases(family, comparison=True))
    summaries.extend(compare_table(
        "fluid_comparison", Path("fluid_comparison") / "comparison.csv", comparison_cases,
        args.reference_dir, args.results_dir, failures))
    report = dict(
        passed=not failures,
        reference_dir=str(args.reference_dir.resolve()),
        results_dir=str(args.results_dir.resolve()),
        expected_sweep_cases=504, expected_finite_sweep_cases=484,
        expected_unavailable_sweep_cases=20, expected_comparison_cases=8,
        tolerance_rule="abs(rerun - reference) <= atol + rtol * abs(reference)",
        relative_difference_rule="abs(rerun - reference) / abs(reference); zero/zero is zero",
        tolerances={metric: dict(rtol=rtol, atol=atol)
                    for metric, (rtol, atol) in METRIC_TOLERANCES.items()},
        metrics=summaries, failures=failures,
    )
    args.results_dir.mkdir(parents=True, exist_ok=True)
    stem = args.results_dir / "reproducibility_comparison"
    with stem.with_suffix(".csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    stem.with_suffix(".json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n",
                                         encoding="utf-8")
    outcome = "PASSED" if report["passed"] else "FAILED"
    print(f"Reproducibility comparison {outcome}: 484 finite sweep cases, "
          f"20 unavailable endpoints, 8 comparison cases; {len(failures)} failures.")
    for summary in summaries:
        if not summary["passed"]:
            print(f"{summary['dataset']}/{summary['metric']}: "
                  f"{summary['failed_cases']} failed cases; "
                  f"max absolute difference {summary['maximum_absolute_difference']}; "
                  f"max relative difference {summary['maximum_relative_difference']}")
    print(f"Detailed comparison: {stem.with_suffix('.json')}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
