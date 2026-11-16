"""Small acceptance check of grid bookkeeping and a few saved native solutions."""

from copy import deepcopy
import json

from fast_evaluator import FastProblem
import run_sensitivity as run

np, pd = run.np, run.pd


def main():
    table = pd.read_csv(run.RESULTS / "comparison.csv")
    expected = len(run.EXPLOITATION_PCT)*len(run.PUMP_EFFICIENCIES)*(1+len(run.BUTANE_TURBINE_EFFICIENCIES))
    assert len(table) == expected, (len(table), expected)
    assert not table.duplicated(["family", "exploitation_pct", "pump_efficiency", "turbine_efficiency"]).any()
    zero = table.loc[table.exploitation_pct.eq(0)]
    assert zero.attempts.eq(0).all() and ~zero.success.any()
    assert zero.relative_gain_pct.isna().all()
    accepted = table.loc[table.success]
    assert accepted.max_constraint_violation.le(run.FEASIBILITY_TOLERANCE).all()
    assert accepted.attempts.between(1, run.MAX_ATTEMPTS).all()
    assert accepted.source_exit_C.ge(accepted.source_min_C-1e-4).all()
    assert accepted.energy_balance_W.abs().lt(0.1).all()
    paired = accepted.loc[accepted.family.eq("butane") & accepted.propane_success]
    assert np.allclose(paired.relative_gain_pct,
                       100*(paired.fixed_reference_efficiency_pct/paired.propane_efficiency_pct-1))
    dry = accepted.loc[accepted.family.eq("propane")]
    assert dry.outlet_phase.eq("vapor").all()
    assert dry.dense_path_vapor_margin_J_kg.ge(-1).all()
    rows = []
    # Six points: both fluids, low/high pump efficiency, and a mid-grid condition.
    for family in run.FLUIDS:
        for exploitation, pump in ((100, 0.2), (100, 1.0), (50, 0.7)):
            selected = accepted.loc[accepted.family.eq(family) & accepted.exploitation_pct.eq(exploitation)
                                   & np.isclose(accepted.pump_efficiency, pump)
                                   & np.isclose(accepted.turbine_efficiency, 0.9)]
            if selected.empty:
                continue
            row = selected.iloc[0]
            directory = run.RESULTS / row["case"]
            checkpoint = json.loads((directory / "result.json").read_text(encoding="utf-8"))
            assert checkpoint["signature"] == run.signature(run.configuration(family, exploitation, pump, 0.9))
            config = run.th.read_configuration_file(directory / "optimal_solution.yaml")
            fixed = config["fixed_parameters"]
            assert config["cycle_topology"] == "simple"
            assert fixed["heater"]["num_elements"] == fixed["cooler"]["num_elements"] == 100
            assert np.isclose(fixed["compressor"]["efficiency"], pump)
            assert np.isclose(fixed["heat_source_pump"]["efficiency"], 0.7)
            native = run.ComparisonProblem(deepcopy(config), out_dir=str(directory))
            fast = FastProblem(deepcopy(config), out_dir=str(directory))
            native.fitness(native.x0)
            fast.fitness(fast.x0)
            objective_error = abs(float(native.f)-float(fast.f))
            constraint_error = float(np.max(np.abs(np.asarray(native.c_ineq)-np.asarray(fast.c_ineq))))
            assert objective_error < 1e-8, objective_error
            assert constraint_error < 1e-6, constraint_error
            assert run.violation(native) <= run.FEASIBILITY_TOLERANCE
            energy = native.cycle_data["energy_analysis"]
            efficiency = 100*float(energy["net_system_power"])/float(energy["mass_flow_heating_fluid"])/run.reference_heat()
            assert abs(efficiency-row.fixed_reference_efficiency_pct) < 1e-7
            assert abs(float(energy["recuperator_heat_flow"])) < 1e-9
            high_bounds = config["design_variables"]["expander_inlet_pressure"]
            assert high_bounds["min"] < row.critical_pressure_bar*1e5 < high_bounds["max"]
            assert np.allclose(native.x0, [checkpoint["variables"][key] for key in native.variable_names], rtol=0, atol=1e-6)
            rows.append(dict(case=row["case"], objective_error=objective_error,
                             max_constraint_difference=constraint_error, efficiency_roundtrip_error_pct=efficiency-row.fixed_reference_efficiency_pct))
    pd.DataFrame(rows).to_csv(run.RESULTS / "verification.csv", index=False)
    print(f"Checked {len(table)} grid records, {len(accepted)} accepted cases, and {len(rows)} native/fast reconstructions.")


if __name__ == "__main__":
    main()
