"""Optimize simple propane/dry and butane/flexible-inlet cycles at a 140 C source.

Run without arguments. Edit the globals below for a smaller grid or plot-only run.
Independent fluid/expander-efficiency sweeps may run in separate processes; each
sweep uses continuation and at most two restarts per point. Checkpoints are reused.
"""

from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
from copy import deepcopy
from dataclasses import replace
from functools import lru_cache
import hashlib
import io
import json
from pathlib import Path
import sys
import time
import traceback

import CoolProp.CoolProp as CP
import jaxprop as props
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysolver_view as psv
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
import thermopt as th

EXPLOITATION_PCT = tuple(range(0, 101, 5))
PUMP_EFFICIENCIES = tuple(round(x / 100, 2) for x in range(20, 101, 10))
BUTANE_TURBINE_EFFICIENCIES = (0.9, 0.8, 0.7, 0.6)
OPTIMIZE = True
WORKERS = 4
MAX_ATTEMPTS = 3  # Primary start, independent start, conditional recovery.
RESTART_GAIN_THRESHOLD = 0.002  # Third start if the first two differ by >0.2%.
PATH_POINTS = 17
FEASIBILITY_TOLERANCE = 1e-5
RESULTS = HERE / "results"
FLUIDS = {"propane": "Propane", "butane": "n-Butane"}


@lru_cache(maxsize=None)
def fluid(name):
    return props.Fluid(name=FLUIDS[name], backend="HEOS", exceptions=True)


@lru_cache(maxsize=1)
def reference_heat():
    return float(CP.PropsSI("H", "P", 1e6, "T", 413.15, "Water")
                 - CP.PropsSI("H", "P", 1e6, "T", 293.15, "Water"))


def expansion_path(problem, count=PATH_POINTS, subcritical_only=False):
    """Cumulative constant-isentropic-efficiency path, including both endpoints.

    h(p) = h_in - eta_t * [h_in - h(p, s_in)]. For propane constraints,
    sample only the saturation-capable part of the path. Above pcrit there
    is no two-phase region; a pseudocritical-superheat constraint would
    unnecessarily exclude physically admissible dense-fluid states.
    """
    expander = problem.cycle_data["components"]["expander"]
    inlet, outlet = expander["state_in"], expander["state_out"]
    name = problem.fixed_parameters["working_fluid"]["name"]
    upper = float(inlet.p)
    if subcritical_only:
        upper = min(upper, 0.999999 * CP.PropsSI("pcrit", name))
    pressure = np.geomspace(upper, float(outlet.p), count)
    ideal_h = CP.PropsSI("H", "P", pressure, "S", float(inlet.s), name)
    enthalpy = float(inlet.h) - expander["efficiency"] * (float(inlet.h) - ideal_h)
    return pressure, enthalpy


class ComparisonProblem(th.ThermodynamicCycleProblem):
    """Native cycle plus propane's sampled vapor-enthalpy constraints."""

    def fitness(self, x):
        result = super().fitness(x)
        # Compute saturation departures from the returned p/h/T. Some jaxprop
        # versions perturb their internal AbstractState while calculating
        # derivative properties of wet states, corrupting these auxiliary fields.
        name = self.fixed_parameters["working_fluid"]["name"]
        components = self.cycle_data["components"]
        for component, endpoint, departure, quality in (
            ("compressor", "state_in", "subcooling", 0),
            ("expander", "state_out", "superheating", 1),
        ):
            state = components[component][endpoint]
            sat_h = CP.PropsSI("H", "P", float(state.p), "Q", quality, name)
            sat_t = CP.PropsSI("T", "P", float(state.p), "Q", quality, name)
            sat_cp = CP.PropsSI("C", "P", float(state.p), "Q", quality, name)
            if quality == 0:
                value = (sat_h-float(state.h))/sat_cp if float(state.h) > sat_h else sat_t-float(state.T)
            else:
                value = (float(state.h)-sat_h)/sat_cp if float(state.h) < sat_h else float(state.T)-sat_t
            components[component][endpoint] = replace(state, **{departure: value})
        self.c_eq, self.c_ineq, self.constraint_data = th.evaluate_constraints(self.cycle_data, self.constraints)
        result = psv.combine_objective_and_constraints(self.f, self.c_eq, self.c_ineq)
        if self.fixed_parameters["working_fluid"]["name"] == FLUIDS["propane"]:
            pressure, enthalpy = expansion_path(self, subcritical_only=True)
            vapor_h = CP.PropsSI("H", "P", pressure, "Q", 1, FLUIDS["propane"])
            self.path_vapor_margin = enthalpy - vapor_h
            self.c_ineq = np.r_[self.c_ineq, -self.path_vapor_margin / 1e5]
            result = psv.combine_objective_and_constraints(self.f, self.c_eq, self.c_ineq)
        return result


class Cycle(th.ThermodynamicCycleOptimization):
    def read_config(self, config):
        return self.load_config(config)

    def setup_problem(self):
        from fast_evaluator import FastProblem
        problem = FastProblem(self.config["problem_formulation"], out_dir=self.out_dir)
        problem.fitness(problem.x0)
        return problem


def case_name(family, exploitation, pump, turbine):
    return f"{family}_E{exploitation:03d}_pump{round(100*pump):03d}_turb{round(100*turbine):02d}"


def configuration(family, exploitation, pump, turbine, seed_kind="high", variables=None):
    config = yaml.safe_load((HERE / "case_template.yaml").read_text(encoding="utf-8"))
    problem = config["problem_formulation"]
    fixed, design = problem["fixed_parameters"], problem["design_variables"]
    fixed["working_fluid"]["name"] = FLUIDS[family]
    fixed["compressor"]["efficiency"] = pump
    fixed["expander"]["efficiency"] = turbine
    minimum = 413.15 - 1.2 * exploitation
    fixed["heat_source"]["minimum_temperature"] = minimum
    design["heat_source_exit_temperature"].update(min=minimum, value=max(minimum, 318.15))
    wf = fluid(family)
    critical_p = float(wf.critical_point.p)
    cold = wf.get_state(props.QT_INPUTS, 0, 273.15)
    warm = wf.get_state(props.QT_INPUTS, 0, 363.15)
    pump_in = wf.get_state(props.PT_INPUTS,
                           float(wf.get_state(props.QT_INPUTS, 0, 303.15).p), 302.15)
    # Both high-pressure bounds straddle the critical pressure, for both fluids.
    # Low-pressure bounds also permit a broad condensing-temperature range.
    design["compressor_inlet_pressure"] = dict(min=float(cold.p), max=float(warm.p), value=float(pump_in.p))
    design["compressor_inlet_enthalpy"] = dict(min=float(cold.h)-20000, max=float(warm.h), value=float(pump_in.h))
    design["expander_inlet_pressure"] = dict(min=1.02*float(cold.p), max=3*critical_p, value=None)
    dilute_h = CP.PropsSI("H", "P", 101325, "T", 413.15, FLUIDS[family])
    design["expander_inlet_enthalpy"] = dict(min=float(cold.h), max=float(dilute_h)+50000, value=None)
    if family == "propane":
        p_high = (1.5 if seed_kind == "high" else 0.72) * critical_p
        inlet = wf.get_state(props.PT_INPUTS, p_high, 398.15 if seed_kind == "high" else 388.15)
        problem["constraints"].append(dict(variable="$components.expander.state_out.superheating",
                                          type=">", value=0, normalize=10))
    else:
        inlet = wf.get_state(props.QT_INPUTS, 0.12 if seed_kind == "high" else 1.0,
                             393.15 if seed_kind == "high" else 363.15)
    design["expander_inlet_pressure"]["value"] = float(inlet.p)
    design["expander_inlet_enthalpy"]["value"] = float(inlet.h)
    if variables:
        for key, value in variables.items():
            design[key]["value"] = float(np.clip(value, design[key]["min"], design[key]["max"]))
    return config


def violation(problem):
    return float(max(np.max(np.abs(problem.c_eq), initial=0), np.max(problem.c_ineq, initial=0)))


def phase(state, family):
    """Use actual saturation enthalpies; generalized Q is ambiguous above pcrit."""
    p, h = float(state.p), float(state.h)
    if p >= float(fluid(family).critical_point.p):
        return "supercritical_pressure", np.nan
    liquid, vapor = CP.PropsSI("H", "P", p, "Q", [0, 1], FLUIDS[family])
    q = (h-liquid)/(vapor-liquid)
    return ("liquid" if q < -1e-7 else "vapor" if q >= 1-1e-7 else "two_phase"), float(q)


def inlet_superheat(pressure, enthalpy, temperature_C, family):
    """True saturation departure; no saturation superheat above critical pressure."""
    if pressure >= float(fluid(family).critical_point.p):
        return np.nan
    name = FLUIDS[family]
    vapor_h = CP.PropsSI("H", "P", pressure, "Q", 1, name)
    if enthalpy < vapor_h:
        return float((enthalpy-vapor_h)/CP.PropsSI("C", "P", pressure, "Q", 1, name))
    return float(temperature_C+273.15-CP.PropsSI("T", "P", pressure, "Q", 1, name))


def measurements(cycle, family):
    solver = cycle.solver
    # Acceptance and all exported numbers use the native ThermoOpt evaluator.
    problem = ComparisonProblem(deepcopy(cycle.problem.configuration), out_dir=cycle.out_dir)
    cycle.problem = problem
    problem.fitness(solver.x_final)
    energy, components = problem.cycle_data["energy_analysis"], problem.cycle_data["components"]
    expander = components["expander"]
    inlet_phase, inlet_q = phase(expander["state_in"], family)
    outlet_phase, outlet_q = phase(expander["state_out"], family)
    valid = all(np.isfinite(energy[k]) and energy[k] > 0 for k in (
        "mass_flow_working_fluid", "mass_flow_heating_fluid", "mass_flow_cooling_fluid",
        "heater_heat_flow", "cooler_heat_flow", "net_system_power"))
    dense_margin = np.nan
    if family == "propane":
        p, h = expansion_path(problem, count=201, subcritical_only=True)
        dense_margin = float(np.min(h - CP.PropsSI("H", "P", p, "Q", 1, FLUIDS[family])))
        valid = valid and dense_margin >= -1.0
    constraint_violation = violation(problem)
    valid = valid and abs(float(energy["energy_balance"])) < 0.1
    high_p = float(expander["state_in"].p)
    low_p = float(expander["state_out"].p)
    pc = float(fluid(family).critical_point.p)
    regime = "transcritical" if low_p < pc < high_p else "subcritical" if high_p < pc else "other"
    return dict(
        success=bool(solver.success and constraint_violation <= FEASIBILITY_TOLERANCE and valid),
        solver_success=bool(solver.success), physical_cycle_valid=bool(valid),
        solver_message=str(solver.message), max_constraint_violation=constraint_violation,
        system_efficiency_pct=100*float(energy["system_efficiency"]),
        fixed_reference_efficiency_pct=100*float(energy["net_system_power"])/float(energy["mass_flow_heating_fluid"])/reference_heat(),
        cycle_efficiency_pct=100*float(energy["cycle_efficiency"]),
        specific_system_power_kW_per_kg_s=float(energy["net_system_power"])/float(energy["mass_flow_heating_fluid"])/1000,
        source_exit_C=float(components["heater"]["hot_side"]["state_out"].T)-273.15,
        actual_exploitation_pct=(413.15-float(components["heater"]["hot_side"]["state_out"].T))/1.2,
        source_mass_flow_kg_s=float(energy["mass_flow_heating_fluid"]),
        high_pressure_bar=high_p/1e5, low_pressure_bar=low_p/1e5,
        pressure_ratio=high_p/low_p, critical_pressure_bar=pc/1e5,
        critical_temperature_C=float(fluid(family).critical_point.T)-273.15,
        pressure_regime=regime, backwork_ratio_pct=100*float(energy["backwork_ratio"]),
        expander_power_MW=float(energy["expander_power"])/1e6,
        working_fluid_pump_power_MW=float(energy["compressor_power"])/1e6,
        auxiliary_pump_power_MW=float(energy["heat_source_pump_power"]+energy["heat_sink_pump_power"])/1e6,
        net_system_power_MW=float(energy["net_system_power"])/1e6,
        inlet_phase=inlet_phase, outlet_phase=outlet_phase,
        inlet_quality=inlet_q, outlet_quality=outlet_q,
        inlet_temperature_C=float(expander["state_in"].T)-273.15,
        inlet_superheat_K=inlet_superheat(high_p, float(expander["state_in"].h),
                                         float(expander["state_in"].T)-273.15, family),
        outlet_superheat_K=float(expander["state_out"].superheating),
        dense_path_vapor_margin_J_kg=dense_margin,
        energy_balance_W=float(energy["energy_balance"]),
        function_evaluations=int(solver.func_count_tot), gradient_evaluations=int(solver.grad_count),
        solver_seconds=float(solver.elapsed_time))


def signature(config):
    schema = deepcopy(config)
    for variable in schema["problem_formulation"]["design_variables"].values():
        variable.pop("value", None)
    return hashlib.sha256(json.dumps(dict(config=schema, path_points=PATH_POINTS,
        max_attempts=MAX_ATTEMPTS, restart_threshold=RESTART_GAIN_THRESHOLD,
        model_version=2), sort_keys=True).encode()).hexdigest()


def write_table(table, path):
    """Replace a completed table atomically so concurrent readers see whole rows."""
    temporary = path.with_suffix(path.suffix+".tmp")
    table.to_csv(temporary, index=False)
    temporary.replace(path)


def solve_point(family, exploitation, pump, turbine, warm=None):
    name = case_name(family, exploitation, pump, turbine)
    directory = RESULTS / name
    directory.mkdir(parents=True, exist_ok=True)
    cache = directory / "result.json"
    base = dict(case=name, family=family, exploitation_pct=exploitation,
                pump_efficiency=pump, turbine_efficiency=turbine, source_min_C=140-1.2*exploitation)
    if exploitation == 0:
        saved = dict(signature="zero_heat", row=dict(base, success=False, attempts=0,
                     status="no_available_heat"), variables={})
        cache.write_text(json.dumps(saved, indent=2), encoding="utf-8")
        return saved
    expected = signature(configuration(family, exploitation, pump, turbine))
    if cache.exists():
        saved = json.loads(cache.read_text(encoding="utf-8"))
        if saved["signature"] != expected:
            raise RuntimeError(f"Inputs changed for {name}; choose a new RESULTS directory.")
        return saved
    attempts, candidates = [], []
    starts = [("continuation" if warm else "high", warm), ("low", None)]
    for index in range(MAX_ATTEMPTS):
        if index == 2:
            accepted = [item for item in candidates if item["row"]["success"]]
            disagreement = len(accepted) == 2 and abs(accepted[0]["row"]["fixed_reference_efficiency_pct"]/
                accepted[1]["row"]["fixed_reference_efficiency_pct"]-1) > RESTART_GAIN_THRESHOLD
            if accepted and not disagreement:
                break
            # Retry the best result with a small, deterministic physical perturbation.
            if candidates:
                best = max(candidates, key=lambda c: (c["row"]["success"], c["row"]["fixed_reference_efficiency_pct"]))
                perturbed = dict(best["variables"])
                perturbed["expander_inlet_pressure"] *= 1.04
                perturbed["expander_inlet_enthalpy"] *= 1.015
                starts.append(("perturbed_best", perturbed))
            else:
                starts.append(("high", None))
        label, initial = starts[index]
        attempt_dir = directory / f"attempt_{index+1}"
        attempt_dir.mkdir(exist_ok=True)
        config = configuration(family, exploitation, pump, turbine,
                               seed_kind="low" if label == "low" else "high", variables=initial)
        log = io.StringIO()
        started = time.perf_counter()
        attempt = dict(base, attempt=index+1, seed=label, success=False)
        try:
            with redirect_stdout(log):
                cycle = Cycle(config, out_dir=str(attempt_dir))
                cycle.problem.save_current_configuration(str(attempt_dir / "initial_guess.yaml"))
                cycle.run_optimization()
                attempt.update(measurements(cycle, family))
                variables = dict(zip(cycle.problem.variable_names, map(float, cycle.solver.x_final)))
                cycle.problem.save_current_configuration(str(attempt_dir / "last_solution.yaml"))
                candidates.append(dict(row=attempt.copy(), variables=variables,
                                       config=deepcopy(cycle.problem.configuration)))
        except (ValueError, RuntimeError, FloatingPointError, OverflowError, np.linalg.LinAlgError) as error:
            attempt.update(solver_message=str(error), status="exception")
            (attempt_dir / "error.txt").write_text(traceback.format_exc(), encoding="utf-8")
        attempt["wall_seconds"] = time.perf_counter()-started
        attempts.append(attempt)
        (attempt_dir / "solver.log").write_text(log.getvalue(), encoding="utf-8")
        write_table(pd.DataFrame(attempts), directory / "attempts.csv")
    if candidates:
        selected = max(candidates, key=lambda c: (c["row"]["success"], c["row"]["fixed_reference_efficiency_pct"]))
        row, variables = selected["row"], selected["variables"]
        if row["success"]:
            (directory / "optimal_solution.yaml").write_text(
                yaml.safe_dump(th.convert_numpy_to_python(selected["config"], precision=12), sort_keys=False), encoding="utf-8")
    else:
        row, variables = attempts[-1], {}
    row.update(attempts=len(attempts), total_solver_seconds=sum(a["wall_seconds"] for a in attempts),
               status="accepted" if row["success"] else "unconverged")
    saved = dict(signature=expected, row=row, variables=variables)
    cache.write_text(json.dumps(saved, indent=2), encoding="utf-8")
    eta = row.get("fixed_reference_efficiency_pct", np.nan)
    print(f"{name}: {row['status']}, eta_fixed={eta:.4f}%, {len(attempts)} starts", flush=True)
    return saved


def run_sweep(task):
    family, turbine = task
    rows, previous_pump = [], {}
    for pump in sorted(PUMP_EFFICIENCIES, reverse=True):
        current_pump, previous_exploitation = {}, None
        for exploitation in sorted(EXPLOITATION_PCT, reverse=True):
            warm = previous_pump.get(exploitation, previous_exploitation)
            saved = solve_point(family, exploitation, pump, turbine, warm)
            rows.append(saved["row"])
            if saved["row"]["success"]:
                current_pump[exploitation] = saved["variables"]
                previous_exploitation = saved["variables"]
            write_table(pd.DataFrame(rows), RESULTS / f"sweep_{family}_turb{round(100*turbine):02d}.csv")
        previous_pump = current_pump
    return rows


def comparison(table):
    table = table.copy()
    # Normalize this auxiliary field when assembling older checkpoints too.
    # Above critical pressure it has no physical saturation-based definition.
    for index, row in table.loc[table.success].iterrows():
        path = RESULTS / row["case"] / "result.json"
        saved = json.loads(path.read_text(encoding="utf-8"))
        value = inlet_superheat(row.high_pressure_bar*1e5, saved["variables"]["expander_inlet_enthalpy"],
                               row.inlet_temperature_C, row.family)
        table.at[index, "inlet_superheat_K"] = value
        saved["row"]["inlet_superheat_K"] = value
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(saved, indent=2), encoding="utf-8")
        temporary.replace(path)
    propane = table.loc[table.family.eq("propane"),
                        ["exploitation_pct", "pump_efficiency", "success", "fixed_reference_efficiency_pct", "pressure_regime"]]
    propane = propane.rename(columns={"success": "propane_success", "fixed_reference_efficiency_pct": "propane_efficiency_pct",
                                      "pressure_regime": "propane_pressure_regime"})
    merged = table.merge(propane, on=["exploitation_pct", "pump_efficiency"], validate="many_to_one")
    merged["relative_gain_pct"] = np.where(
        merged.success & merged.propane_success & merged.propane_efficiency_pct.gt(0),
        100*(merged.fixed_reference_efficiency_pct/merged.propane_efficiency_pct-1), np.nan)
    merged["efficiency_difference_pp"] = np.where(merged.success & merged.propane_success,
        merged.fixed_reference_efficiency_pct-merged.propane_efficiency_pct, np.nan)
    write_table(merged, RESULTS / "comparison.csv")
    return merged


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    tasks = [("propane", 0.9)] + [("butane", eta) for eta in BUTANE_TURBINE_EFFICIENCIES]
    if OPTIMIZE:
        if WORKERS == 1:
            batches = [run_sweep(task) for task in tasks]
        else:
            with ProcessPoolExecutor(max_workers=WORKERS) as executor:
                batches = list(executor.map(run_sweep, tasks))
        table = pd.DataFrame([row for batch in batches for row in batch])
        write_table(table, RESULTS / "sweep_results.csv")
    else:
        table = pd.read_csv(RESULTS / "sweep_results.csv")
    table = comparison(table)
    write_table(table.drop(columns=["propane_success", "propane_efficiency_pct", "propane_pressure_regime", "relative_gain_pct", "efficiency_difference_pp"]),
                RESULTS / "sweep_results.csv")
    import plot_results
    plot_results.plot_all(table)
    import build_report
    build_report.make_report(table)
    print(table.groupby(["family", "turbine_efficiency"]).success.agg(["sum", "count"]))


if __name__ == "__main__":
    main()
