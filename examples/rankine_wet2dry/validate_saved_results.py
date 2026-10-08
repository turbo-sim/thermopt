"""Reevaluate all three completed sweeps without optimizing or changing saved inputs.

Run after all three sweeps and the fluid-comparison export finish. The selected
YAML/JSON files and sweep CSVs are fingerprinted before and after evaluation.
Only audit tables and summaries are written. The 501-point heat-exchanger
calculations diagnose sampling error; they do not replace the 100-point results.
"""

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_fluid_comparison as run

np, pd, th = run.np, run.pd, run.th
RESULTS = Path(__file__).resolve().parent / "results"
FAMILIES = ("butane", "cyclopentane", "novec")
SOURCE_MINIMUMS_C = tuple(range(140, 19, -5))
TURBINE_EFFICIENCIES = (0.90, 0.80, 0.70, 0.60)
DENSE_TEMPERATURES_C = (135, 70, 20)
DENSE_EFFICIENCIES = (0.90, 0.60)
QUALITY_TOLERANCE = 1e-7
MONOTONIC_WARNING_PCT = 0.10  # Relative loss of specific system power.
MONOTONIC_MATERIAL_PCT = 1.0  # Diagnostic only: single-start local solutions are retained.


def expansion_class(inlet_quality, outlet_quality):
    """Classify from raw Q, whose native extension exceeds one in dry vapor."""
    def phase(quality):
        if quality <= QUALITY_TOLERANCE:
            return "liquid"
        return "wet" if quality < 1 - QUALITY_TOLERANCE else "dry"

    inlet, outlet = phase(inlet_quality), phase(outlet_quality)
    return f"{inlet}_to_{outlet}"


def fingerprint(paths):
    return {str(path.relative_to(RESULTS)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in paths}


def evaluate(family, output, row):
    directory = output / row["case"]
    config = th.read_configuration_file(directory / "optimal_solution.yaml")
    fixed, variables = config["fixed_parameters"], config["design_variables"]
    fluid_name, topology, _ = run.FAMILIES[family]
    assert config["cycle_topology"] == topology
    assert fixed["working_fluid"]["name"] == fluid_name
    assert fixed["working_fluid"]["backend"] == "HEOS"
    assert fixed["heating_fluid"]["name"] == fixed["cooling_fluid"]["name"] == "Water"
    assert fixed["net_power"] == 1e6
    assert fixed["heat_source"]["inlet_temperature"] == 413.15
    assert fixed["heat_source"]["inlet_pressure"] == fixed["heat_source"]["exit_pressure"] == 1e6
    assert np.isclose(fixed["heat_source"]["minimum_temperature"], row["source_min_C"] + 273.15)
    assert fixed["heat_sink"]["inlet_temperature"] == 293.15
    assert fixed["heat_sink"]["inlet_pressure"] == fixed["heat_sink"]["exit_pressure"] == 101325
    assert np.isclose(fixed["expander"]["efficiency"], row["turbine_efficiency"])
    assert fixed["expander"]["efficiency_type"] == "isentropic"
    assert np.isclose(variables["heat_source_exit_temperature"]["min"], row["source_min_C"] + 273.15)
    assert np.isclose(variables["heat_source_exit_temperature"]["max"], 413.14)
    assert np.isclose(variables["heat_sink_exit_temperature"]["min"], 298.14)
    assert np.isclose(variables["heat_sink_exit_temperature"]["max"], 298.15)
    if topology == "recuperated":
        assert variables["recuperator_effectiveness"]["min"] == 0
        assert variables["recuperator_effectiveness"]["max"] == 1
    exchanger_names = ("heater", "cooler", "recuperator") if topology == "recuperated" else ("heater", "cooler")

    def has_constraint(variable, direction, value):
        return any(c["variable"] == variable and c["type"] == direction and c["value"] == value
                   for c in config["constraints"])

    for name in exchanger_names:
        assert fixed[name]["num_elements"] == (25 if name == "recuperator" else 100)
        assert fixed[name]["pressure_drop_hot_side"] == fixed[name]["pressure_drop_cold_side"] == 0.01
        assert has_constraint(f"$components.{name}.temperature_difference", ">", 5)
    for name in ("compressor", "heat_source_pump", "heat_sink_pump"):
        assert fixed[name]["efficiency"] == 0.70
    assert has_constraint("$components.compressor.state_in.subcooling", ">", 1)
    assert has_constraint("$components.compressor.state_in.Q", "<", 0)
    if row["mode"] == "superheated":
        assert has_constraint("$components.expander.state_in.superheating", "=", 2)
        assert has_constraint("$components.expander.state_out.superheating", ">", 0)
    else:
        assert not any("$components.expander.state_in." in c["variable"]
                       for c in config["constraints"]), "Partial evaporation has a fixed inlet constraint."
    assert config["objective_function"] == dict(variable="$energy_analysis.system_efficiency", type="maximize")

    problem = th.ThermodynamicCycleProblem(config, out_dir=str(output))
    problem.fitness(problem.x0)
    assert run.physical_cycle_valid(problem)
    violation = run.violation(problem)
    assert violation <= 1e-5
    energy, components = problem.cycle_data["energy_analysis"], problem.cycle_data["components"]
    inlet, outlet = components["expander"]["state_in"], components["expander"]["state_out"]
    efficiency = 100 * float(energy["system_efficiency"])
    assert abs(efficiency - row["system_efficiency_pct"]) < 1e-7
    assert abs(float(inlet.Q) - row["inlet_quality"]) < 1e-7
    assert abs(float(outlet.Q) - row["outlet_quality"]) < 1e-7
    assert abs(float(energy["energy_balance"])) < 0.1
    assert np.isclose(float(energy["net_cycle_power"]), 1e6, rtol=0, atol=1e-5)
    assert np.isclose(float(energy["net_system_power"]) / 1e6, row["net_system_power_MW"], rtol=0, atol=1e-10)
    assert float(energy["net_system_power"]) < float(energy["net_cycle_power"])
    if row["mode"] == "superheated":
        assert abs(float(inlet.superheating) - 2) < 2e-5
        assert float(outlet.superheating) >= -2e-5
        assert float(outlet.Q) >= 1 - QUALITY_TOLERANCE

    source_exit_C = float(components["heater"]["hot_side"]["state_out"].T) - 273.15
    cooling_exit_C = float(components["cooler"]["cold_side"]["state_out"].T) - 273.15
    assert abs(source_exit_C - row["source_exit_C"]) < 1e-7
    assert source_exit_C >= row["source_min_C"] - 1e-6
    assert 24.99 - 1e-6 <= cooling_exit_C <= 25.00 + 1e-6
    assert abs(float(energy["mass_flow_heating_fluid"]) - row["source_mass_flow_kg_s"]) < 1e-6
    specific_power = float(energy["net_system_power"] / energy["mass_flow_heating_fluid"]) / 1000
    assert abs(specific_power - row["specific_system_power_kW_per_kg_s"]) < 1e-7
    assert abs(100000 * specific_power / run.fixed_reference_heat_J_per_kg() - row["fixed_reference_efficiency_pct"]) < 1e-8
    fraction = (problem.x0 - problem.lb) / (np.array(problem.ub) - problem.lb)
    distances = np.minimum(fraction, 1 - fraction)
    assert distances.min() >= -1e-8
    state_indices = [i for i, name in enumerate(problem.variable_names)
                     if name.startswith(("compressor_", "expander_"))]
    distance = float(distances[state_indices].min())
    saved = json.loads((directory / "result.json").read_text(encoding="utf-8"))
    assert saved["row"]["success"]
    assert np.allclose([saved["variables"][name] for name in problem.variable_names],
                       problem.x0, rtol=0, atol=1e-7)
    assert abs(saved["row"]["system_efficiency_pct"] - efficiency) < 1e-7
    checks = dict(
        case=row["case"], mode=row["mode"], source_min_C=row["source_min_C"],
        turbine_efficiency=row["turbine_efficiency"],
        efficiency_roundtrip_error_pct=efficiency - row["system_efficiency_pct"],
        max_constraint_violation=violation, minimum_state_bound_distance=distance,
        minimum_all_bound_distance=float(distances.min()),
        state_bounds_near_limit=";".join(problem.variable_names[i] for i in state_indices if distances[i] < 1e-4),
        inlet_quality=float(inlet.Q), outlet_quality=float(outlet.Q),
        expander_inlet_pressure_bar=float(inlet.p) / 1e5,
        expander_inlet_temperature_C=float(inlet.T) - 273.15,
        expander_outlet_temperature_C=float(outlet.T) - 273.15,
        inlet_superheat_K=float(inlet.superheating), outlet_superheat_K=float(outlet.superheating),
        expansion_class=expansion_class(float(inlet.Q), float(outlet.Q)),
        recuperator_effectiveness=float(variables["recuperator_effectiveness"]["value"])
            if topology == "recuperated" else np.nan,
        source_exit_C=source_exit_C, cooling_exit_C=cooling_exit_C,
        energy_balance_W=float(energy["energy_balance"]),
    )

    pinches = None
    dense = row["source_min_C"] in DENSE_TEMPERATURES_C and any(
        np.isclose(row["turbine_efficiency"], value) for value in DENSE_EFFICIENCIES)
    if dense:
        pinches = dict(case=row["case"])
        for name in exchanger_names:
            pinches[name + "_model_K"] = float(np.min(components[name]["temperature_difference"]))
        refined = deepcopy(config)
        for name in exchanger_names:
            refined["fixed_parameters"][name]["num_elements"] = 501
        dense_problem = th.ThermodynamicCycleProblem(refined, out_dir=str(output))
        dense_problem.fitness(dense_problem.x0)
        for name in exchanger_names:
            pinches[name + "_501_K"] = float(np.min(
                dense_problem.cycle_data["components"][name]["temperature_difference"]))
        coarse = th.read_configuration_file(run.HERE / "seed_results" / f"sensitivity_{family}"
                                             / row["case"] / "optimal_solution.yaml")
        for name in exchanger_names:
            coarse["fixed_parameters"][name]["num_elements"] = 501
        coarse_problem = th.ThermodynamicCycleProblem(coarse, out_dir=str(output))
        coarse_problem.fitness(coarse_problem.x0)
        for name in exchanger_names:
            pinches[name + "_coarse_501_K"] = float(np.min(
                coarse_problem.cycle_data["components"][name]["temperature_difference"]))
    return checks, pinches, exchanger_names


def audit_single_starts(family, output, table):
    count = 0
    for row in table.loc[table["source_min_C"].lt(140)].itertuples():
        folder = output / row.case
        initial = th.read_configuration_file(folder / "initial_guess.yaml")
        seed_folder = run.HERE / "seed_results" / f"sensitivity_{family}" / row.case
        seed = th.read_configuration_file(seed_folder / "optimal_solution.yaml")
        expected = deepcopy(seed)
        for component in ("heater", "cooler"):
            expected["fixed_parameters"][component]["num_elements"] = 100
        for name, variable in initial["design_variables"].items():
            assert np.isclose(variable["value"], expected["design_variables"][name]["value"], rtol=0, atol=1e-7)
            variable["value"] = expected["design_variables"][name]["value"]
        # The configuration reader converts numeric lists to NumPy arrays.
        canonical = lambda config: json.dumps(config, sort_keys=True, default=lambda value: value.tolist())
        assert canonical(initial) == canonical(expected), f"Initial state or physics changed beyond HX grids: {row.case}"
        marker = json.loads((folder / "attempt_started.json").read_text(encoding="utf-8"))
        assert marker["optimization_calls"] == 1
        assert row.attempts == 1
        attempts = pd.read_csv(folder / "attempts.csv")
        assert len(attempts) == 1 and attempts.iloc[0]["attempt"] == 1
        assert not list(folder.glob("attempt_*[0-9]"))
        count += 1
    assert count == 120
    return count


def audit_family(family):
    output = RESULTS / f"sensitivity_{family}"
    table = pd.read_csv(output / "sweep_results.csv", float_precision="round_trip")
    expected = {(mode, temperature, efficiency)
                for temperature in SOURCE_MINIMUMS_C
                for mode, efficiencies in (("superheated", (0.90,)),
                                          ("partial_evaporation", TURBINE_EFFICIENCIES))
                for efficiency in efficiencies}
    assert set(zip(table["mode"], table["source_min_C"], table["turbine_efficiency"])) == expected
    assert len(table) == 125 and table["case"].is_unique
    plan = pd.read_csv(output / "execution_plan.csv")
    assert len(plan) == 125 and set(plan["case"]) == set(table["case"])
    single_starts_verified = audit_single_starts(family, output, table)
    endpoint = table.loc[table["source_min_C"].eq(140)]
    finite = table.loc[table["source_min_C"].lt(140)]
    accepted = finite.loc[finite["success"]]
    failed = finite.loc[~finite["success"]]
    assert len(endpoint) == 5 and endpoint["status"].eq("no_geothermal_temperature_drop").all()
    assert not endpoint["success"].any() and endpoint["attempts"].eq(0).all()
    assert endpoint[["system_efficiency_pct", "fixed_reference_efficiency_pct", "relative_gain_pct", "source_exit_C",
                     "actual_exploitation_pct", "inlet_quality", "outlet_quality"]].isna().all().all()
    assert len(accepted) + len(failed) == 120 and accepted["success"].all()
    assert accepted["solver_success"].all() and accepted["physical_cycle_valid"].all()
    assert np.allclose(accepted["net_cycle_power_MW"], 1.0)
    assert np.allclose(table["exploitation_pct"], 100 * (140 - table["source_min_C"]) / 120)
    assert np.allclose(accepted["actual_exploitation_pct"], 100 * (140 - accepted["source_exit_C"]) / 120)
    assert np.all(accepted["actual_exploitation_pct"] <= accepted["exploitation_pct"] + 1e-6)
    for side in ("inlet", "outlet"):
        assert np.allclose(accepted[f"{side}_quality_plot"], accepted[f"{side}_quality"].clip(0, 1))

    paths = [output / "sweep_results.csv", output / "execution_plan.csv"]
    paths.extend(output / case / "result.json" for case in table["case"])
    paths.extend(output / case / "optimal_solution.yaml" for case in accepted["case"])
    before = fingerprint(paths)
    for case in endpoint["case"]:
        directory = output / case
        checkpoint = json.loads((directory / "result.json").read_text(encoding="utf-8"))
        assert checkpoint["variables"] == {} and checkpoint["row"]["attempts"] == 0
        assert not (directory / "optimal_solution.yaml").exists()
        assert not list(directory.glob("attempt_*"))

    baseline = finite.loc[finite["mode"].eq("superheated")].set_index("source_min_C")
    checks, pinches = [], []
    for _, row in accepted.iterrows():
        reference = baseline.loc[row["source_min_C"]]
        if reference["success"]:
            efficiency_gain = 100 * (row["system_efficiency_pct"] / reference["system_efficiency_pct"] - 1)
            power_gain = 100 * (row["specific_system_power_kW_per_kg_s"] / reference["specific_system_power_kW_per_kg_s"] - 1)
            fixed_gain = 100 * (row["fixed_reference_efficiency_pct"] / reference["fixed_reference_efficiency_pct"] - 1)
            assert abs(efficiency_gain - row["relative_gain_pct"]) < 1e-8
            assert abs(power_gain - row["relative_gain_pct"]) < 1e-8
            assert abs(fixed_gain - row["relative_gain_pct"]) < 1e-8
        else:
            assert pd.isna(row["relative_gain_pct"])
        check, pinch, exchanger_names = evaluate(family, output, row)
        checks.append(check)
        if pinch is not None:
            pinches.append(pinch)
        print(f"Checked {family}/{row['case']}: {check['expansion_class']}", flush=True)
    assert before == fingerprint(paths), "A saved input changed during validation."
    check_table, pinch_table = pd.DataFrame(checks), pd.DataFrame(pinches)
    assert len(pinch_table) == len(accepted.loc[accepted["source_min_C"].isin(DENSE_TEMPERATURES_C) & accepted["turbine_efficiency"].isin(DENSE_EFFICIENCIES)])
    check_table.to_csv(output / "saved_solution_checks.csv", index=False)
    pinch_table.to_csv(output / "pinch_grid_check.csv", index=False)

    anomalies = []
    for (mode, efficiency), group in accepted.groupby(["mode", "turbine_efficiency"]):
        best_power, best_case = -np.inf, None
        for _, row in group.sort_values("source_min_C", ascending=False).iterrows():
            power = float(row["specific_system_power_kW_per_kg_s"])
            if power < best_power:
                loss = 100 * (1 - power / best_power)
                if loss > 1e-6:
                    anomalies.append(dict(case=row["case"], earlier_case=best_case,
                                          relative_specific_power_loss_pct=loss,
                                          severity="material" if loss > MONOTONIC_MATERIAL_PCT else
                                          "warning" if loss > MONOTONIC_WARNING_PCT else "numerical"))
            if power > best_power:
                best_power, best_case = power, row["case"]

    comparison_checks = []
    comparison_path = RESULTS / "fluid_comparison" / "comparison.csv"
    if comparison_path.exists():
        comparison = pd.read_csv(comparison_path, float_precision="round_trip")
        for mode in ("superheated", "partial_evaporation"):
            selected = comparison.loc[comparison["case"].eq(f"{family}_{mode}")]
            assert len(selected) == 1
            if not bool(selected.iloc[0]["success"]):
                continue
            current = accepted.loc[accepted["mode"].eq(mode) & accepted["source_min_C"].eq(70)
                                   & np.isclose(accepted["turbine_efficiency"], 0.90)].iloc[0]
            difference = float(selected.iloc[0]["system_efficiency_pct"] - current["system_efficiency_pct"])
            assert abs(difference) < 1e-7
            comparison_checks.append(dict(mode=mode, efficiency_difference_percentage_points=difference))

    summary = dict(
        family=family, requested_slots=len(table), accepted_sweep_points=len(accepted),
        unavailable_zero_source_drop_points=len(endpoint), failed_sweep_points=len(failed),
        failed_cases=failed["case"].tolist(),
        baseline_dominated_cases=accepted.loc[accepted["baseline_dominated"], "case"].tolist(),
        saved_inputs_unchanged=True, fingerprinted_input_files=len(paths),
        single_optimization_per_finite_case_verified=True,
        exact_coarse_initial_states_verified=single_starts_verified,
        fixed_reference_heat_J_per_kg=run.fixed_reference_heat_J_per_kg(),
        template_optimizer_tolerance=1e-6,
        optimizer_tolerance_note="Verified in shared template; selected YAML stores only the cycle problem.",
        strict_kkt_passes=int(accepted["all_kkt_checks_pass"].sum()),
        max_constraint_violation=float(check_table["max_constraint_violation"].max()),
        max_energy_balance_error_W=float(check_table["energy_balance_W"].abs().max()),
        maximum_efficiency_roundtrip_error_pct=float(check_table["efficiency_roundtrip_error_pct"].abs().max()),
        minimum_state_bound_distance=float(check_table["minimum_state_bound_distance"].min()),
        cases_with_state_variable_near_bound=check_table.loc[
            check_table["state_bounds_near_limit"].ne(""), "case"].tolist(),
        expansion_class_counts={mode: {key: int(value) for key, value in subset["expansion_class"].value_counts().items()}
                                for mode, subset in check_table.groupby("mode")},
        quality_range={side: [float(check_table[f"{side}_quality"].min()),
                              float(check_table[f"{side}_quality"].max())] for side in ("inlet", "outlet")},
        dense_grid_points=len(pinches),
        minimum_sampled_approach_K={name: float(pinch_table[name + "_model_K"].min()) for name in exchanger_names},
        minimum_dense_approach_K={name: float(pinch_table[name + "_501_K"].min()) for name in exchanger_names},
        minimum_coarse_dense_approach_K={name: float(pinch_table[name + "_coarse_501_K"].min()) for name in exchanger_names},
        dense_grid_note="Diagnostic only: selected solutions use 100 heater/condenser samples and 25 recuperator samples.",
        monotonic_specific_power_anomalies=anomalies,
        monotonic_warning_threshold_pct=MONOTONIC_WARNING_PCT,
        monotonic_material_threshold_pct=MONOTONIC_MATERIAL_PCT,
        comparison_at_Tmin70_eta90=comparison_checks,
    )
    (output / "validation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)
    return summary


def main():
    template = th.read_configuration_file(run.HERE / "case_template.yaml")
    assert template["solver_options"]["tolerance"] == 1e-6
    summaries = [audit_family(family) for family in FAMILIES]
    assert sum(item["requested_slots"] for item in summaries) == 375
    assert sum(item["accepted_sweep_points"] + item["failed_sweep_points"] for item in summaries) == 360
    print("Three studies verified: 360 single attempts and 15 unavailable 140 C endpoints.", flush=True)


if __name__ == "__main__":
    main()
