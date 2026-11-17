"""Export the eight 70 C designs at Cristina's 88.6% turbine efficiency."""
from copy import deepcopy
import json
import shutil
import study as run

np, pd, plt, th = run.np, run.pd, run.plt, run.th


def export_case(family, mode):
    case = run.point_name(mode, 70, run.BASELINE_EFFICIENCY)
    source = run.RESULTS / f"sensitivity_{family}" / case
    directory = run.RESULTS / "fluid_comparison" / f"{family}_{mode}"
    directory.mkdir(parents=True, exist_ok=True)
    sweep = pd.read_csv(run.RESULTS / f"sensitivity_{family}" / "sweep_results.csv")
    row = sweep.loc[sweep["case"].eq(case)].iloc[0].to_dict()
    if not row["success"]:
        raise RuntimeError(f"Comparison needs an accepted solution: {family}/{case}")
    for file in ("optimal_solution.yaml", "convergence_history.txt", "optimization_report.txt", "attempts.csv", "result.json"):
        shutil.copy2(source / file, directory / file)
    problem = th.ThermodynamicCycleProblem(run.read_formulation(source / "optimal_solution.yaml"), out_dir=str(directory))
    problem.fitness(problem.x0)
    assert run.physical_cycle_valid(problem) and run.violation(problem) < run.FEASIBILITY_TOLERANCE
    assert abs(100*problem.cycle_data["energy_analysis"]["system_efficiency"]-row["system_efficiency_pct"]) < 1e-7
    bounds = pd.DataFrame({"variable": problem.variable_names, "lower": problem.lb,
                           "value": problem.x0, "upper": problem.ub})
    fraction = (bounds["value"]-bounds["lower"])/(bounds["upper"]-bounds["lower"])
    bounds["distance_to_nearest_bound_fraction"] = np.minimum(fraction, 1-fraction)
    bounds["active_bound"] = np.where(fraction < 1e-4, "lower", np.where(fraction > 1-1e-4, "upper", ""))
    bounds.to_csv(directory / "design_variable_bounds.csv", index=False)
    problem.save_data_to_excel(filename=str(directory / "optimal_solution.xlsx"))
    problem.plot_cycle()
    problem.figure.suptitle(f"{run.FAMILIES[family][2]} - {mode.replace('_', ' ')}\nTurbine efficiency 88.6%; minimum source outlet 70 C")
    problem.figure.tight_layout(rect=(0, 0, 1, .93))
    run.save_figure(problem.figure, directory / "optimal_solution")
    plt.close(problem.figure)
    history = []
    for line in (source / "convergence_history.txt").read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) != 6:
            continue
        try:
            values = [float(x) for x in fields]
        except ValueError:
            continue
        history.append(values)
    if history:
        h = np.array(history)
        figure, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
        axes[0].plot(h[:, 1], h[:, 2], color=th.COLORS_MATLAB[1])
        axes[0].set_ylabel("Minimization objective")
        axes[1].semilogy(h[:, 1], np.maximum(h[:, 3], 1e-14), color=th.COLORS_MATLAB[0])
        axes[1].set_ylabel("Maximum constraint violation")
        for ax in axes:
            ax.set_xlabel("Function evaluations")
        figure.suptitle(f"{run.FAMILIES[family][2]}: {mode.replace('_', ' ')}")
        run.save_figure(figure, directory / "convergence_history")
        plt.close(figure)
    state_bounds = bounds[bounds["variable"].str.startswith(("compressor_", "expander_"))]
    row = dict(row, case=f"{family}_{mode}", source_case=case,
               minimum_state_bound_distance=float(state_bounds["distance_to_nearest_bound_fraction"].min()),
               active_bounds=", ".join(bounds.loc[bounds["active_bound"].ne(""), "variable"]))
    pd.DataFrame([row]).to_csv(directory / "summary.csv", index=False)
    return problem, row


def main():
    rows = []
    output = run.RESULTS / "fluid_comparison"
    output.mkdir(exist_ok=True)
    for family in run.FAMILIES:
        figure, axes = plt.subplots(2, 2, figsize=(10.5, 9.6), squeeze=False)
        for column, mode in enumerate(("superheated", "partial_evaporation")):
            problem, row = export_case(family, mode)
            rows.append(row)
            # Fresh instance ensures native plot artists belong to this grid.
            problem = th.ThermodynamicCycleProblem(deepcopy(problem.configuration))
            problem.fitness(problem.x0)
            problem.figure, problem.axes = figure, axes[:, column]
            problem.plot_cycle()
            label = "2 K inlet superheat" if mode == "superheated" else "Optimized inlet state"
            axes[0, column].set_title(label)
        figure.suptitle(f"{run.FAMILIES[family][2]}\n88.6% turbine efficiency; 70 C minimum source outlet")
        figure.tight_layout(rect=(0, .045, 1, .94))
        from matplotlib.lines import Line2D
        figure.legend(handles=[Line2D([], [], color=th.COLORS_MATLAB[i], label=label)
                      for i, label in ((1, "Working fluid"), (6, "Heat source"), (0, "Cooling water"))],
                      loc="lower center", ncol=3)
        run.save_figure(figure, output / f"comparison_grid_{family}")
        plt.close(figure)
    table = pd.DataFrame(rows)
    baseline = table.loc[table["case"].eq("butane_superheated"), "specific_system_power_kW_per_kg_s"].iloc[0]
    table["gain_vs_standard_butane_pct"] = 100*(table["specific_system_power_kW_per_kg_s"]/baseline-1)
    table.to_csv(output / "comparison.csv", index=False)
    print(table[["case", "specific_system_power_kW_per_kg_s", "relative_gain_pct", "gain_vs_standard_butane_pct"]].to_string(index=False))


if __name__ == "__main__":
    main()
