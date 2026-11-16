"""Recover failures and known feasible-reference losses with native ThermoOpt."""
from contextlib import redirect_stdout
from copy import deepcopy
import io
import json
import shutil
import time
import traceback
import study as run

np, pd = run.np, run.pd


def refine(family, row, variables, reason):
    directory = run.RESULTS / f"sensitivity_{family}" / row["case"]
    previous = json.loads((directory / "result.json").read_text(encoding="utf-8"))
    if reason in previous.get("recovery_reasons", []):
        return previous
    index = len(list(directory.glob("recovery_*")))+1
    attempt_dir = directory / f"recovery_{index}"
    attempt_dir.mkdir()
    (attempt_dir / "previous_result.json").write_text(json.dumps(previous, indent=2), encoding="utf-8")
    if (directory / "optimal_solution.yaml").exists():
        shutil.copy2(directory / "optimal_solution.yaml", attempt_dir / "previous_optimal_solution.yaml")
    config = run.configuration(family, row["mode"], int(row["source_min_C"]), row["turbine_efficiency"], variables=variables)
    attempt = dict(row, success=False, seed=reason)
    output, start = io.StringIO(), time.perf_counter()
    try:
        with redirect_stdout(output):
            cycle = run.create_native_cycle(config, attempt_dir)
            cycle.problem.save_current_configuration(str(attempt_dir / "initial_guess.yaml"))
            cycle.run_optimization()
            attempt.update(run.measurements(cycle))
            cycle.problem.save_current_configuration(str(attempt_dir / "last_solution.yaml"))
            cycle.print_convergence_history(savefile=True)
            cycle.print_optimization_report(savefile=True)
        if attempt["success"] and (not previous["row"]["success"] or attempt["specific_system_power_kW_per_kg_s"] > previous["row"]["specific_system_power_kW_per_kg_s"]):
            previous["row"].update(attempt)
            previous["row"]["status"] = "accepted"
            previous["variables"] = dict(zip(cycle.problem.variable_names, map(float, cycle.solver.x_final)))
            for src, dst in (("last_solution.yaml", "optimal_solution.yaml"), ("convergence_history.txt", "convergence_history.txt"), ("optimization_report.txt", "optimization_report.txt")):
                shutil.copy2(attempt_dir / src, directory / dst)
            previous["backend"]["input_file"] = f"recovery_{index}/case.yaml"
    except (ValueError, RuntimeError, FloatingPointError, OverflowError, np.linalg.LinAlgError) as error:
        attempt.update(solver_message=str(error))
        (attempt_dir / "error.txt").write_text(traceback.format_exc(), encoding="utf-8")
    attempt["wall_seconds"] = time.perf_counter()-start
    (attempt_dir / "solver.log").write_text(output.getvalue(), encoding="utf-8")
    (attempt_dir / "result.json").write_text(json.dumps(attempt, indent=2), encoding="utf-8")
    history = pd.read_csv(directory / "attempts.csv")
    history = pd.concat([history, pd.DataFrame([attempt])], ignore_index=True)
    history.to_csv(directory / "attempts.csv", index=False)
    previous["row"]["attempts"] = len(history)
    previous["row"]["total_solver_seconds"] = float(history["wall_seconds"].sum())
    previous.setdefault("recovery_reasons", []).append(reason)
    (directory / "result.json").write_text(json.dumps(previous, indent=2), encoding="utf-8")
    print(f"Refined {family}/{row['case']} ({reason}): {previous['row'].get('fixed_reference_efficiency_pct', np.nan):.4f}%", flush=True)
    return previous


def family(family):
    directory = run.RESULTS / f"sensitivity_{family}"
    table = pd.read_csv(directory / "sweep_results.csv")
    previous_best = {}
    for row in table.sort_values(["source_min_C", "mode"], ascending=[False, False]).to_dict("records"):
        if row["source_min_C"] == 140:
            continue
        saved = json.loads((directory / row["case"] / "result.json").read_text(encoding="utf-8"))
        baseline_case = run.point_name("superheated", int(row["source_min_C"]), run.BASELINE_EFFICIENCY)
        baseline = json.loads((directory / baseline_case / "result.json").read_text(encoding="utf-8"))
        if row["mode"] == "partial_evaporation" and baseline["row"]["success"]:
            dominated = row["turbine_efficiency"] >= .886 and (not saved["row"]["success"] or saved["row"]["specific_system_power_kW_per_kg_s"] < baseline["row"]["specific_system_power_kW_per_kg_s"]*.9999)
            if dominated or not saved["row"]["success"]:
                saved = refine(family, saved["row"], baseline["variables"], "standard_cycle_start")
        key = (row["mode"], row["turbine_efficiency"])
        best = previous_best.get(key)
        if best and (not saved["row"]["success"] or saved["row"]["specific_system_power_kW_per_kg_s"] < best["row"]["specific_system_power_kW_per_kg_s"]*.9999):
            saved = refine(family, saved["row"], best["variables"], "stricter_limit_start")
        if saved["row"]["success"] and (not best or saved["row"]["specific_system_power_kW_per_kg_s"] > best["row"]["specific_system_power_kW_per_kg_s"]):
            previous_best[key] = saved
    rows = [json.loads((directory / name / "result.json").read_text(encoding="utf-8"))["row"] for name in table["case"]]
    updated = run.add_plot_metrics(pd.DataFrame(rows))
    updated.to_csv(directory / "sweep_results.csv", index=False)
    print(f"Refinement completed: {family}", flush=True)


def main():
    if len(run.sys.argv) == 3 and run.sys.argv[1] == "--family":
        family(run.sys.argv[2])
        return
    import subprocess
    run.LOGS.mkdir(parents=True, exist_ok=True)
    processes = []
    for f in run.FAMILIES:
        log = (run.LOGS / f"refinement_{f}.log").open("w", encoding="utf-8")
        process = subprocess.Popen([run.sys.executable, str(run.SCRIPTS / "refine_results.py"), "--family", f], stdout=log, stderr=subprocess.STDOUT)
        processes.append((f, process, log))
    for f, process, log in processes:
        code = process.wait()
        log.close()
        if code:
            raise RuntimeError(f"Refinement failed: {f}")


if __name__ == "__main__":
    main()
