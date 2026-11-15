"""Shared configuration, optimization and native verification for the study."""

from thermopt.utilities.fluid_utils import get_property_fluid
import jaxprop as props
from contextlib import redirect_stdout
from copy import deepcopy
from functools import lru_cache
import hashlib
import html
import io
import json
from pathlib import Path
import sys
import time
import traceback
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore", message="FigureCanvasAgg is non-interactive")
import numpy as np
import pandas as pd
import yaml

from paths import STUDY as HERE, SCRIPTS, PROJECT, OUTPUT, RESULTS, LOGS, ARCHIVE, REFERENCE_INPUTS
sys.path.insert(0, str(PROJECT))
import thermopt as th

FAMILIES = {
    "butane": ("n-Butane", "simple", "n-Butane: simple cycle"),
    "pentane": ("n-Pentane", "simple", "n-Pentane: simple cycle"),
    "cyclopentane": ("Cyclopentane", "simple", "Cyclopentane: simple cycle"),
    "novec": ("Novec649", "recuperated", "Novec 649: recuperated cycle"),
}
BASELINE_EFFICIENCY = 0.886
TURBINE_EFFICIENCIES = (0.90, 0.80, 0.70, 0.60)
SOURCE_MINIMUMS_C = tuple(range(140, 19, -5))
INITIALIZERS = HERE / "inputs" / "initializers"
ARCHIVED_RESULTS = ARCHIVE / "coolprop_evaluator_2026-11-02" / "results"
OPTIMIZATION_BACKEND = "thermopt.native"
MAX_ATTEMPTS = 3
FEASIBILITY_TOLERANCE = 1e-5


def point_name(mode, minimum, efficiency):
    eta = f"{100*efficiency:g}".replace(".", "p")
    return f"{mode}_Tmin{minimum:03d}_eta{eta}"


@lru_cache(maxsize=1)
def fixed_reference_heat_J_per_kg():
    return float(get_property_fluid('Water').get_state(props.PT_INPUTS, 1000000.0, 413.15).h
                 - get_property_fluid('Water').get_state(props.PT_INPUTS, 1000000.0, 288.15).h)


def configuration(family, mode, minimum, efficiency, variables=None, seed="saved"):
    """Build physics from the template; saved designs supply initial values only."""
    config = th.convert_numpy_to_python(th.read_configuration_file(HERE / "case_template.yaml"), precision=14)
    problem = config["problem_formulation"]
    fixed, design = problem["fixed_parameters"], problem["design_variables"]
    name, topology, _ = FAMILIES[family]
    problem["cycle_topology"] = topology
    fixed["working_fluid"]["name"] = name
    fixed["heat_source"]["minimum_temperature"] = minimum+273.15
    fixed["expander"]["efficiency"] = efficiency
    # Retain the original broad subcritical state bounds, recomputed per fluid.
    pc = float(get_property_fluid(name).critical_point.p)
    pmin = get_property_fluid(name).get_state(props.QT_INPUTS, 0, 273.15).p
    hmin = get_property_fluid(name).get_state(props.QT_INPUTS, 0, 273.15).h
    hliqmax = get_property_fluid(name).get_state(props.QT_INPUTS, 0, 413.15).h
    hmax = get_property_fluid(name).get_state(props.PT_INPUTS, pmin, 413.15).h+50000
    pp = get_property_fluid(name).get_state(props.QT_INPUTS, 0, 302.65).p
    ph = get_property_fluid(name).get_state(props.PT_INPUTS, pp, 301.65).h
    ti = min(403.15, max(345., minimum+273.15-12))
    hp = get_property_fluid(name).get_state(props.QT_INPUTS, 0, ti).p
    hh = get_property_fluid(name).get_state(props.PT_INPUTS, hp, ti + 2).h if mode == "superheated" else get_property_fluid(name).get_state(props.QT_INPUTS, 0.3, ti).h
    design["heat_source_exit_temperature"].update(min=minimum+273.15, value=max(minimum+273.15, 343.15))
    design["compressor_inlet_pressure"] = dict(min=pmin, max=.95*pc, value=pp)
    design["compressor_inlet_enthalpy"] = dict(min=hmin, max=hliqmax, value=ph)
    design["expander_inlet_pressure"] = dict(min=pmin, max=.95*pc, value=hp)
    design["expander_inlet_enthalpy"] = dict(min=hmin, max=hmax, value=hh)
    c = problem["constraints"]
    c.extend([
        dict(variable="$components.cooler.cold_side.state_out.T", type="=", value=298.15, normalize=1.0),
        dict(variable="$components.expander.pressure_ratio", type=">", value=1., normalize=True),
        dict(variable="$components.expander.specific_work - $components.compressor.specific_work", type=">", value=0., normalize=10000.),
        dict(variable="$components.heater.cold_side.state_out.h - $components.heater.cold_side.state_in.h", type=">", value=0., normalize=100000.),
        dict(variable="$components.cooler.hot_side.state_in.h - $components.cooler.hot_side.state_out.h", type=">", value=0., normalize=100000.),
    ])
    if mode == "superheated":
        c.extend([
            dict(variable="$components.expander.state_in.superheating", type="=", value=2., normalize=True),
            dict(variable="$components.expander.state_out.superheating", type=">", value=0., normalize=2.),
        ])
    if topology == "recuperated":
        fixed["recuperator"] = dict(pressure_drop_hot_side=.01, pressure_drop_cold_side=.01, num_elements=25)
        design["recuperator_effectiveness"] = dict(min=0., max=1., value=.5)
        c.append(dict(variable="$components.recuperator.temperature_difference", type=">", value=5., normalize=True))
    if seed == "saved":
        path = INITIALIZERS / f"sensitivity_{family}" / point_name(mode, minimum, efficiency) / "optimal_solution.yaml"
        if path.exists():
            old = yaml.safe_load(path.read_text(encoding="utf-8"))["design_variables"]
            for variable_name, variable in design.items():
                if variable_name in old:
                    variable["value"] = float(np.clip(old[variable_name]["value"], variable["min"], variable["max"]))
    if variables:
        for k, value in variables.items():
            if k in design:
                design[k]["value"] = float(np.clip(value, design[k]["min"], design[k]["max"]))
    design["heat_sink_exit_temperature"]["value"] = 298.15
    return config


@lru_cache(maxsize=1)
def backend_sources():
    """Fingerprint the unchanged native model and component implementations."""
    package = Path(th.__file__).resolve().parent
    paths = [package / "optimize_cycle.py", *sorted((package / "cycles").glob("*.py")),
             *sorted((package / "components").glob("*.py"))]
    return {str(p.relative_to(package)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def write_provenance():
    import importlib.metadata as metadata
    import platform
    sources = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
               for p in REFERENCE_INPUTS.glob("*") if p.is_file()}
    record = dict(optimization_backend=OPTIMIZATION_BACKEND,
                  optimization_entrypoint="thermopt.ThermodynamicCycleOptimization",
                  problem_class="thermopt.optimize_cycle.ThermodynamicCycleProblem",
                  python=platform.python_version(), native_source_sha256=backend_sources(),
                  packages={p: metadata.version(p) for p in ("thermopt", "jaxprop", "pysolver_view", "CoolProp", "scipy", "numpy", "matplotlib")},
                  source_files_sha256=sources,
                  initialization="Versioned inputs/initializers states; every finite-heat case is optimized by native ThermoOpt.")
    (RESULTS / "environment_and_sources.json").write_text(json.dumps(record, indent=2), encoding="utf-8")


def create_native_cycle(config, directory):
    """Use the public ThermoOpt entry point directly, with a saved YAML input.

    No problem subclass, fitness override, surrogate or alternate cycle evaluator
    is installed. This exact input can be run independently with ThermoOpt.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    case_file = directory / "case.yaml"
    case_file.write_text(yaml.safe_dump(th.convert_numpy_to_python(config, precision=14), sort_keys=False), encoding="utf-8")
    cycle = th.ThermodynamicCycleOptimization(str(case_file), out_dir=str(directory))
    assert type(cycle.problem) is th.ThermodynamicCycleProblem
    return cycle


def violation(problem):
    return float(max(np.max(np.abs(problem.c_eq), initial=0), np.max(problem.c_ineq, initial=0)))


def physical_cycle_valid(problem):
    energy = problem.cycle_data["energy_analysis"]
    e = problem.cycle_data["components"]["expander"]
    return (all(np.isfinite(energy[k]) and energy[k] > 0 for k in (
        "mass_flow_working_fluid", "mass_flow_heating_fluid", "mass_flow_cooling_fluid",
        "heater_heat_flow", "cooler_heat_flow", "net_system_power"))
        and e["state_in"].p > e["state_out"].p and abs(energy["energy_balance"]) < .1)


def measurements(cycle):
    native, solver = cycle.problem, cycle.solver
    assert type(native) is th.ThermodynamicCycleProblem
    native.fitness(solver.x_final)
    e, c = native.cycle_data["energy_analysis"], native.cycle_data["components"]
    valid = physical_cycle_valid(native)
    vio = violation(native)
    source_flow, net = float(e["mass_flow_heating_fluid"]), float(e["net_system_power"])
    row = dict(success=bool(solver.success and vio <= FEASIBILITY_TOLERANCE and valid),
        solver_success=bool(solver.success), solver_message=str(solver.message),
        physical_cycle_valid=bool(valid), max_constraint_violation=vio,
        optimization_backend=OPTIMIZATION_BACKEND,
        problem_class=type(native).__module__+"."+type(native).__name__,
        optimization_class=type(cycle).__module__+"."+type(cycle).__name__,
        optimization_called=True,
        system_efficiency_pct=100*float(e["system_efficiency"]), cycle_efficiency_pct=100*float(e["cycle_efficiency"]),
        fixed_reference_efficiency_pct=100*net/source_flow/fixed_reference_heat_J_per_kg(),
        inlet_quality=float(c["expander"]["state_in"].Q), inlet_superheat_K=float(c["expander"]["state_in"].superheating),
        outlet_quality=float(c["expander"]["state_out"].Q),
        source_exit_C=float(c["heater"]["hot_side"]["state_out"].T)-273.15,
        cooling_exit_C=float(c["cooler"]["cold_side"]["state_out"].T)-273.15,
        cooling_rise_K=float(c["cooler"]["cold_side"]["state_out"].T-c["cooler"]["cold_side"]["state_in"].T),
        source_mass_flow_kg_s=source_flow, net_cycle_power_MW=float(e["net_cycle_power"])/1e6,
        net_system_power_MW=net/1e6, specific_system_power_kW_per_kg_s=net/source_flow/1000,
        power_at_194p44_kg_s_kW=net/source_flow/1000*194.44,
        achieved_heat_utilization_pct=100*float(e["heater_heat_flow"])/float(e["heater_heat_flow_max"]),
        energy_balance_W=float(e["energy_balance"]), function_evaluations=int(solver.func_count_tot),
        gradient_evaluations=int(solver.grad_count), solver_seconds=float(solver.elapsed_time),
        turbine_inlet_pressure_bar=float(c["expander"]["state_in"].p)/1e5,
        condenser_inlet_pressure_bar=float(c["cooler"]["hot_side"]["state_in"].p)/1e5,
        inlet_temperature_C=float(c["expander"]["state_in"].T)-273.15)
    for key in ("heater", "cooler") + (("recuperator",) if native.cycle_topology == "recuperated" else ()):
        row[key+"_model_K"] = float(np.min(c[key]["temperature_difference"]))
    return row


def signature(config):
    config = deepcopy(config)
    for d in config["problem_formulation"]["design_variables"].values():
        d.pop("value", None)
    return hashlib.sha256(json.dumps(dict(config=config, backend=OPTIMIZATION_BACKEND,
                                         native_sources=backend_sources()), sort_keys=True).encode()).hexdigest()


def solve_point(family, mode, minimum, efficiency, warm=None):
    name = point_name(mode, minimum, efficiency)
    directory = RESULTS / f"sensitivity_{family}" / name
    directory.mkdir(parents=True, exist_ok=True)
    cache = directory / "result.json"
    base = dict(case=name, family=family, mode=mode, source_min_C=minimum,
                turbine_efficiency=efficiency, exploitation_pct=100*(140-minimum)/125)
    if minimum == 140:
        saved = dict(signature="zero_heat", row=dict(base, success=False, status="no_geothermal_temperature_drop", attempts=0), variables={})
        cache.write_text(json.dumps(saved, indent=2), encoding="utf-8")
        return saved
    expected = signature(configuration(family, mode, minimum, efficiency))
    if cache.exists():
        saved = json.loads(cache.read_text(encoding="utf-8"))
        if saved["signature"] != expected:
            raise RuntimeError(f"Changed physical settings for {family}/{name}; use a new results directory.")
        return saved
    archived = INITIALIZERS / f"sensitivity_{family}" / name / "optimal_solution.yaml"
    archived_variables = None
    if archived.exists():
        previous = th.read_configuration_file(archived)
        archived_variables = {key: float(item["value"]) for key, item in previous["design_variables"].items()}
    starts = [("archived_initial_guess" if archived_variables else "saved", archived_variables),
              ("continuation", warm), ("generic", None)]
    attempts, candidates = [], []
    for index, (label, initial) in enumerate(starts[:MAX_ATTEMPTS]):
        if index and candidates and candidates[-1]["row"]["success"]:
            break
        if index == 1 and initial is None:
            continue
        attempt_dir = directory / f"attempt_{index+1}"
        attempt_dir.mkdir(exist_ok=True)
        config = configuration(family, mode, minimum, efficiency, variables=initial, seed="generic" if label == "generic" else "saved")
        attempt = dict(base, success=False, attempt=index+1, seed=label)
        output = io.StringIO()
        start = time.perf_counter()
        try:
            with redirect_stdout(output):
                cycle = create_native_cycle(config, attempt_dir)
                cycle.problem.save_current_configuration(str(attempt_dir / "initial_guess.yaml"))
                cycle.run_optimization()
                attempt.update(measurements(cycle))
                variables = dict(zip(cycle.problem.variable_names, map(float, cycle.solver.x_final)))
                cycle.problem.save_current_configuration(str(attempt_dir / "last_solution.yaml"))
                cycle.print_convergence_history(savefile=True)
                cycle.print_optimization_report(savefile=True)
                candidates.append(dict(row=attempt.copy(), variables=variables, config=deepcopy(cycle.problem.configuration), path=attempt_dir))
        except (ValueError, RuntimeError, FloatingPointError, OverflowError, np.linalg.LinAlgError) as error:
            attempt.update(solver_message=str(error), status="exception")
            (attempt_dir / "error.txt").write_text(traceback.format_exc(), encoding="utf-8")
        attempt["wall_seconds"] = time.perf_counter()-start
        attempts.append(attempt)
        (attempt_dir / "solver.log").write_text(output.getvalue(), encoding="utf-8")
        pd.DataFrame(attempts).to_csv(directory / "attempts.csv", index=False)
    if candidates:
        selected = max(candidates, key=lambda c: (c["row"]["success"], c["row"]["specific_system_power_kW_per_kg_s"]))
        row, variables = selected["row"], selected["variables"]
        if row["success"]:
            import shutil
            shutil.copy2(selected["path"] / "last_solution.yaml", directory / "optimal_solution.yaml")
            for file in ("convergence_history.txt", "optimization_report.txt"):
                shutil.copy2(selected["path"] / file, directory / file)
    else:
        row, variables = attempts[-1], {}
    row.update(attempts=len(attempts), status="accepted" if row["success"] else "unconverged",
               total_solver_seconds=sum(a["wall_seconds"] for a in attempts))
    saved = dict(signature=expected, row=row, variables=variables)
    saved["backend"] = dict(name=OPTIMIZATION_BACKEND, native_sources=backend_sources(),
                            input_file=f"attempt_{row.get('attempt', 1)}/case.yaml",
                            initializer=str(archived.relative_to(HERE)) if archived_variables else "original study",
                            initializer_sha256=hashlib.sha256(archived.read_bytes()).hexdigest() if archived_variables else None)
    cache.write_text(json.dumps(saved, indent=2), encoding="utf-8")
    print(f"{family}/{name}: {row['status']}, eta_fixed={row.get('fixed_reference_efficiency_pct', np.nan):.4f}%", flush=True)
    return saved


def add_plot_metrics(table):
    table = table.copy()
    table["fixed_reference_efficiency_pct"] = 100000*table["specific_system_power_kW_per_kg_s"]/fixed_reference_heat_J_per_kg()
    table["actual_exploitation_pct"] = 100*(140-table["source_exit_C"])/125
    table["inlet_quality_plot"] = table["inlet_quality"].clip(0, 1)
    table["outlet_quality_plot"] = table["outlet_quality"].clip(0, 1)
    baseline = table[table["mode"].eq("superheated")].set_index("source_min_C")
    table["relative_gain_pct"] = np.nan
    for i, row in table.iterrows():
        if row["source_min_C"] not in baseline.index:
            continue
        ref = baseline.loc[row["source_min_C"]]
        if row["success"] and ref["success"]:
            table.loc[i, "relative_gain_pct"] = 100*(row["system_efficiency_pct"]/ref["system_efficiency_pct"]-1)
    return table


def run_family(family):
    rows = []
    directory = RESULTS / f"sensitivity_{family}"
    directory.mkdir(parents=True, exist_ok=True)
    warm = {}
    for minimum in SOURCE_MINIMUMS_C:
        cases = [("superheated", BASELINE_EFFICIENCY)]+[("partial_evaporation", e) for e in TURBINE_EFFICIENCIES]
        if minimum == 70:
            cases.append(("partial_evaporation", BASELINE_EFFICIENCY))
        for mode, efficiency in cases:
            key = (mode, efficiency)
            saved = solve_point(family, mode, minimum, efficiency, warm.get(key))
            rows.append(saved["row"])
            if saved["row"]["success"]:
                warm[key] = saved["variables"]
        table = pd.DataFrame(rows)
        if "specific_system_power_kW_per_kg_s" in table:
            table = add_plot_metrics(table)
        table.to_csv(directory / "sweep_results.csv", index=False)
    return family, int(table["success"].sum()), len(table)


def save_figure(figure, filename, formats=(".png", ".pdf"), **kwargs):
    """Identical Matplotlib content in PNG, PDF, SVG and self-contained HTML."""
    filename = Path(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)
    th.savefig_in_formats(figure, filename, formats=list(formats), **kwargs)
    svg = io.StringIO()
    figure.savefig(svg, format="svg", bbox_inches="tight")
    svgtext = svg.getvalue()
    filename.with_suffix(".svg").write_text(svgtext, encoding="utf-8")
    title = html.escape(filename.stem.replace("_", " "))
    body = svgtext[svgtext.index("<svg"):]
    page = f'<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>{title}</title><style>html,body{{margin:0;background:white}}svg{{display:block;width:100%;height:auto}}</style><body>{body}</body></html>'
    filename.with_suffix(".html").write_text(page, encoding="utf-8")


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    if len(sys.argv) == 3 and sys.argv[1] == "--family":
        print("Completed:", run_family(sys.argv[2]), flush=True)
        return
    write_provenance()
    # File-backed subprocesses also work where Windows named pipes are denied.
    import subprocess
    processes = []
    for family in FAMILIES:
        log = (LOGS / f"optimization_{family}.log").open("w", encoding="utf-8")
        process = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), "--family", family],
                                   stdout=log, stderr=subprocess.STDOUT)
        processes.append((family, process, log))
    failures = []
    for family, process, log in processes:
        code = process.wait()
        log.close()
        print(f"Finished {family}: exit {code}", flush=True)
        if code:
            failures.append(family)
    if failures:
        raise RuntimeError(f"Failed fluid studies: {failures}; inspect output/logs/optimization_*.log")


if __name__ == "__main__":
    main()
