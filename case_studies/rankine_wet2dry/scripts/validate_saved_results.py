"""Reevaluate every accepted design and audit saved native records and sweep data."""
from copy import deepcopy
import hashlib
import json
import study as run
from paths import resolve_record_path

np, pd, th = run.np, run.pd, run.th
QUALITY_TOLERANCE = 1e-7
DENSE_TEMPERATURES_C = (135, 70, 20)
MONOTONIC_WARNING_PCT = 0.10
MONOTONIC_MATERIAL_PCT = 1.0


def fingerprint(paths):
    return {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def physics(config):
    """Compare physical inputs independently of initial and final design values."""
    config = deepcopy(config)
    for variable in config["design_variables"].values():
        variable.pop("value", None)
    return th.convert_numpy_to_python(config, precision=10)


def expansion_class(inlet_quality, outlet_quality):
    def phase(quality):
        if quality <= QUALITY_TOLERANCE:
            return "liquid"
        return "wet" if quality < 1 - QUALITY_TOLERANCE else "dry"

    return f"{phase(inlet_quality)} to {phase(outlet_quality)}"


def recorded_kkt(path):
    """Read native solver diagnostics without changing the acceptance criteria."""
    labels = ("First order optimality", "Equality feasibility", "Inequality feasibility",
              "Dual feasibility", "Complementary slackness")
    values, passes = {}, []
    for line in path.read_text(encoding="utf-8").splitlines():
        for label in labels:
            if line.lstrip().startswith(label):
                fields = line.split()
                values[label.lower().replace(" ", "_")] = float(fields[-4])
                passes.append(fields[-1] == "yes")
    return dict(values, kkt_report_available=len(passes) == len(labels),
                all_kkt_checks_pass=len(passes) == len(labels) and all(passes))


def evaluate(family, directory, row, saved):
    path = directory / row["case"] / "optimal_solution.yaml"
    config = run.read_formulation(path)
    expected = run.configuration(family, row["mode"], int(row["source_min_C"]),
                                 row["turbine_efficiency"], seed="generic")
    assert saved["signature"] == run.signature(expected), row["case"]
    assert physics(config) == physics(expected["problem_formulation"]), row["case"]
    backend = saved["backend"]
    assert backend["name"] == run.OPTIMIZATION_BACKEND
    assert backend["native_sources"] == run.backend_sources()
    if backend.get("initializer_sha256"):
        initializer = resolve_record_path(backend["initializer"])
        assert hashlib.sha256(initializer.read_bytes()).hexdigest() == backend["initializer_sha256"]
    case_file = path.parent / backend["input_file"]
    native_input = th.read_configuration_file(case_file)
    assert physics(native_input["problem_formulation"]) == physics(config)
    assert native_input["solver_options"] == expected["solver_options"]
    initial = run.read_formulation(case_file.parent / "initial_guess.yaml")
    assert physics(initial) == physics(config)
    for name, variable in initial["design_variables"].items():
        value = native_input["problem_formulation"]["design_variables"][name]["value"]
        assert np.isclose(variable["value"], value, rtol=1e-12, atol=1e-7), (row["case"], name)
    log = (case_file.parent / "convergence_history.txt").read_text(encoding="utf-8")
    assert "ThermodynamicCycleProblem" in log and "FastProblem" not in log
    assert row["optimization_backend"] == run.OPTIMIZATION_BACKEND
    assert row["problem_class"] == "thermopt.optimize_cycle.ThermodynamicCycleProblem"
    assert row["optimization_class"] == "thermopt.optimize_cycle.ThermodynamicCycleOptimization"
    assert row["optimization_called"] and row["function_evaluations"] > 0

    problem = th.ThermodynamicCycleProblem(config)
    problem.fitness(problem.x0)
    assert run.physical_cycle_valid(problem), row["case"]
    violation = run.violation(problem)
    assert violation <= run.FEASIBILITY_TOLERANCE, row["case"]
    energy, components = problem.cycle_data["energy_analysis"], problem.cycle_data["components"]
    inlet, outlet = components["expander"]["state_in"], components["expander"]["state_out"]
    source_exit = float(components["heater"]["hot_side"]["state_out"].T) - 273.15
    cooling_exit = float(components["cooler"]["cold_side"]["state_out"].T) - 273.15
    specific_power = float(energy["net_system_power"] / energy["mass_flow_heating_fluid"]) / 1000
    current = dict(
        system_efficiency_pct=100 * float(energy["system_efficiency"]),
        cycle_efficiency_pct=100 * float(energy["cycle_efficiency"]),
        fixed_reference_efficiency_pct=100000 * specific_power / run.fixed_reference_heat_J_per_kg(),
        inlet_quality=float(inlet.Q), outlet_quality=float(outlet.Q),
        inlet_superheat_K=float(inlet.superheating), source_exit_C=source_exit,
        cooling_exit_C=cooling_exit, source_mass_flow_kg_s=float(energy["mass_flow_heating_fluid"]),
        net_cycle_power_MW=float(energy["net_cycle_power"]) / 1e6,
        net_system_power_MW=float(energy["net_system_power"]) / 1e6,
        specific_system_power_kW_per_kg_s=specific_power,
        power_at_194p44_kg_s_kW=specific_power * 194.44,
        turbine_inlet_pressure_bar=float(inlet.p) / 1e5,
        condenser_inlet_pressure_bar=float(components["cooler"]["hot_side"]["state_in"].p) / 1e5,
        inlet_temperature_C=float(inlet.T) - 273.15,
    )
    for name, value in current.items():
        assert np.isclose(value, row[name], rtol=1e-10, atol=1e-7), (row["case"], name, value, row[name])
    assert abs(float(energy["energy_balance"])) < 0.1
    assert abs(float(energy["energy_balance"]) - row["energy_balance_W"]) < 1e-4
    assert np.isclose(current["net_cycle_power_MW"], 1.0, rtol=0, atol=1e-10)
    assert current["net_system_power_MW"] < current["net_cycle_power_MW"]
    assert source_exit >= row["source_min_C"] - 1e-6
    assert abs(cooling_exit - 25) < 2e-5
    if row["mode"] == "superheated":
        assert abs(float(inlet.superheating) - 2) < 2e-5
        assert float(outlet.superheating) >= -2e-5
        assert float(outlet.Q) >= 1 - QUALITY_TOLERANCE
    assert np.allclose([saved["variables"][name] for name in problem.variable_names],
                       problem.x0, rtol=0, atol=1e-7), row["case"]
    fraction = (problem.x0 - problem.lb) / (np.asarray(problem.ub) - problem.lb)
    distances = np.minimum(fraction, 1 - fraction)
    assert distances.min() >= -1e-8, row["case"]
    state_indices = [i for i, name in enumerate(problem.variable_names)
                     if name.startswith(("compressor_", "expander_"))]
    check = dict(
        case=row["case"], sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        mode=row["mode"], source_min_C=row["source_min_C"], turbine_efficiency=row["turbine_efficiency"],
        energy_balance_W=float(energy["energy_balance"]), max_constraint_violation=violation,
        efficiency_roundtrip_error_pct=current["system_efficiency_pct"] - row["system_efficiency_pct"],
        minimum_state_bound_distance=float(distances[state_indices].min()),
        minimum_all_bound_distance=float(distances.min()),
        state_bounds_near_limit=";".join(problem.variable_names[i] for i in state_indices if distances[i] < 1e-4),
        inlet_quality=float(inlet.Q), outlet_quality=float(outlet.Q),
        inlet_superheat_K=float(inlet.superheating), outlet_superheat_K=float(outlet.superheating),
        expansion_class=expansion_class(float(inlet.Q), float(outlet.Q)),
        source_exit_C=source_exit, cooling_exit_C=cooling_exit,
        optimization_backend=row["optimization_backend"], input_file=str(case_file.relative_to(run.HERE)),
    )
    check.update(recorded_kkt(path.parent / "optimization_report.txt"))
    pinch = None
    if row["source_min_C"] in DENSE_TEMPERATURES_C and (
        row["mode"] == "superheated" or row["turbine_efficiency"] in (.9, .6, run.BASELINE_EFFICIENCY)
    ):
        exchangers = ("heater", "cooler", "recuperator") if family == "novec" else ("heater", "cooler")
        dense = deepcopy(config)
        for name in exchangers:
            dense["fixed_parameters"][name]["num_elements"] = 501
        refined = th.ThermodynamicCycleProblem(dense)
        refined.fitness(refined.x0)
        pinch = dict(case=row["case"])
        for name in exchangers:
            pinch[name + "_model_K"] = float(np.min(components[name]["temperature_difference"]))
            pinch[name + "_501_K"] = float(np.min(refined.cycle_data["components"][name]["temperature_difference"]))
    return check, pinch


def audit_family(family):
    directory = run.RESULTS / f"sensitivity_{family}"
    table = pd.read_csv(directory / "sweep_results.csv", float_precision="round_trip")
    expected = {(mode, temperature, efficiency)
                for temperature in run.SOURCE_MINIMUMS_C
                for mode, efficiencies in (("superheated", (run.BASELINE_EFFICIENCY,)),
                                          ("partial_evaporation", run.TURBINE_EFFICIENCIES))
                for efficiency in efficiencies}
    expected.add(("partial_evaporation", 70, run.BASELINE_EFFICIENCY))
    assert len(table) == len(expected) == 126 and table["case"].is_unique
    assert set(zip(table["mode"], table["source_min_C"], table["turbine_efficiency"])) == expected
    assert table["family"].eq(family).all()
    assert table["case"].tolist() == [run.point_name(row["mode"], int(row["source_min_C"]), row["turbine_efficiency"])
                                      for row in table.to_dict("records")]
    endpoints = table.loc[table["source_min_C"].eq(140)]
    finite = table.loc[table["source_min_C"].lt(140)]
    accepted = finite.loc[finite["success"]]
    failed = finite.loc[~finite["success"]]
    assert len(endpoints) == 5 and endpoints["status"].eq("no_geothermal_temperature_drop").all()
    assert not endpoints["success"].any() and endpoints["attempts"].eq(0).all()
    assert endpoints[["system_efficiency_pct", "fixed_reference_efficiency_pct", "relative_gain_pct",
                      "source_exit_C", "actual_exploitation_pct", "inlet_quality", "outlet_quality"]].isna().all().all()
    assert len(accepted) + len(failed) == 121
    assert accepted["solver_success"].all() and accepted["physical_cycle_valid"].all()
    assert np.allclose(table["exploitation_pct"], 100 * (140 - table["source_min_C"]) / 125)
    assert np.allclose(accepted["actual_exploitation_pct"], 100 * (140 - accepted["source_exit_C"]) / 125)
    assert np.all(accepted["actual_exploitation_pct"] <= accepted["exploitation_pct"] + 1e-6)
    for side in ("inlet", "outlet"):
        assert np.allclose(accepted[f"{side}_quality_plot"], accepted[f"{side}_quality"].clip(0, 1))
    paths = [path for path in directory.rglob("*") if path.is_file()
             and (path.suffix in (".yaml", ".json", ".txt") or path.name in ("sweep_results.csv", "attempts.csv"))
             and path.name != "validation_summary.json"]
    before = fingerprint(paths)
    records = {}
    for row in table.to_dict("records"):
        folder = directory / row["case"]
        saved = json.loads((folder / "result.json").read_text(encoding="utf-8"))
        records[row["case"]] = saved
        for name, value in saved["row"].items():
            assert name in row, (row["case"], name)
            if isinstance(value, bool):
                assert bool(row[name]) == value, (row["case"], name)
            elif isinstance(value, (int, float)):
                assert np.isclose(row[name], value, rtol=1e-12, atol=1e-8, equal_nan=True), (row["case"], name)
            else:
                assert row[name] == value, (row["case"], name)
        if row["source_min_C"] == 140:
            assert saved["variables"] == {} and saved["signature"] == "zero_heat"
            assert not (folder / "optimal_solution.yaml").exists()
            assert not list(folder.glob("attempt_*"))
        else:
            attempts = pd.read_csv(folder / "attempts.csv")
            assert len(attempts) == row["attempts"]
    baseline = finite.loc[finite["mode"].eq("superheated")].set_index("source_min_C")
    checks, pinches = [], []
    for row in accepted.to_dict("records"):
        reference = baseline.loc[row["source_min_C"]]
        if reference["success"]:
            for metric in ("system_efficiency_pct", "fixed_reference_efficiency_pct", "specific_system_power_kW_per_kg_s"):
                gain = 100 * (row[metric] / reference[metric] - 1)
                assert abs(gain - row["relative_gain_pct"]) < 1e-8, (row["case"], metric)
        else:
            assert pd.isna(row["relative_gain_pct"])
        check, pinch = evaluate(family, directory, row, records[row["case"]])
        checks.append(check)
        if pinch is not None:
            pinches.append(pinch)
        print(f"Checked {family}/{row['case']}: {check['expansion_class']}", flush=True)
    assert before == fingerprint(paths), "A saved native record or sweep changed during validation."
    check_table, pinch_table = pd.DataFrame(checks), pd.DataFrame(pinches)
    check_table.to_csv(directory / "saved_solution_checks.csv", index=False)
    pinch_table.to_csv(directory / "pinch_grid_check.csv", index=False)
    anomalies = []
    for _, group in accepted.groupby(["mode", "turbine_efficiency"]):
        best_power, best_case = -np.inf, None
        for row in group.sort_values("source_min_C", ascending=False).to_dict("records"):
            power = row["specific_system_power_kW_per_kg_s"]
            loss = 100 * (1 - power / best_power) if np.isfinite(best_power) else 0
            if loss > 1e-6:
                severity = "numerical"
                if loss > MONOTONIC_WARNING_PCT:
                    severity = "warning"
                if loss > MONOTONIC_MATERIAL_PCT:
                    severity = "material"
                anomalies.append(dict(case=row["case"], earlier_case=best_case,
                                      relative_specific_power_loss_pct=loss, severity=severity))
            if power > best_power:
                best_power, best_case = power, row["case"]
    comparison_checks = []
    comparison_path = run.RESULTS / "fluid_comparison" / "comparison.csv"
    if comparison_path.exists():
        comparison = pd.read_csv(comparison_path, float_precision="round_trip")
        for mode in ("superheated", "partial_evaporation"):
            exported = comparison.loc[comparison["case"].eq(f"{family}_{mode}")]
            assert len(exported) == 1
            current = accepted.loc[accepted["mode"].eq(mode) & accepted["source_min_C"].eq(70)
                                   & np.isclose(accepted["turbine_efficiency"], run.BASELINE_EFFICIENCY)].iloc[0]
            assert exported.iloc[0]["source_case"] == current["case"]
            for metric in ("system_efficiency_pct", "specific_system_power_kW_per_kg_s", "relative_gain_pct"):
                assert abs(exported.iloc[0][metric] - current[metric]) < 1e-7, (family, mode, metric)
            comparison_checks.append(dict(mode=mode, source_case=current["case"],
                                          efficiency_difference_percentage_points=float(exported.iloc[0]["system_efficiency_pct"] - current["system_efficiency_pct"])))
    partial = accepted.loc[accepted["mode"].eq("partial_evaporation")
                           & accepted["turbine_efficiency"].isin(run.TURBINE_EFFICIENCIES)]
    counts = pd.Series([expansion_class(row.inlet_quality, row.outlet_quality)
                        for row in partial.itertuples()]).value_counts().to_dict()
    summary = dict(
        family=family, requested_points=len(table), accepted_points=len(accepted),
        unavailable_points=len(endpoints), failed_points=len(failed), failed_cases=failed["case"].tolist(),
        dense_points=len(pinch_table), expansion_classes={key: int(value) for key, value in counts.items()},
        optimization_backend=run.OPTIMIZATION_BACKEND, native_optimization_points=len(accepted),
        maximum_constraint_violation=float(check_table["max_constraint_violation"].max()),
        minimum_dense_approach_K={name.removesuffix("_501_K"): float(pinch_table[name].min())
                                  for name in pinch_table if name.endswith("_501_K")},
        saved_inputs_unchanged=True, fingerprinted_input_files=len(paths),
        independently_reevaluated_points=len(check_table),
        strict_kkt_passes=int(check_table["all_kkt_checks_pass"].sum()),
        cases_with_recorded_kkt=int(check_table["kkt_report_available"].sum()),
        max_energy_balance_error_W=float(check_table["energy_balance_W"].abs().max()),
        maximum_efficiency_roundtrip_error_pct=float(check_table["efficiency_roundtrip_error_pct"].abs().max()),
        minimum_state_bound_distance=float(check_table["minimum_state_bound_distance"].min()),
        cases_with_state_variable_near_bound=check_table.loc[check_table["state_bounds_near_limit"].ne(""), "case"].tolist(),
        monotonic_specific_power_anomalies=anomalies,
        monotonic_warning_threshold_pct=MONOTONIC_WARNING_PCT,
        monotonic_material_threshold_pct=MONOTONIC_MATERIAL_PCT,
        comparison_at_Tmin70_eta88p6=comparison_checks,
    )
    (directory / "validation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(summary, flush=True)
    return summary


def main():
    summaries = [audit_family(family) for family in run.FAMILIES]
    assert sum(item["requested_points"] for item in summaries) == 504
    assert sum(item["accepted_points"] + item["failed_points"] for item in summaries) == 484
    assert sum(item["unavailable_points"] for item in summaries) == 20
    print("Verified four complete studies and independently reevaluated every accepted design.", flush=True)


if __name__ == "__main__":
    main()
