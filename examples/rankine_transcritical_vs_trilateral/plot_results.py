"""Reconstruct saved cycles and draw the pump/exploitation comparison maps."""

from pathlib import Path

from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator

import run_sensitivity as run

np, pd, plt = run.np, run.pd, run.plt
OUTPUT = run.RESULTS / "figures"


def save(figure, name):
    OUTPUT.mkdir(parents=True, exist_ok=True)
    run.th.savefig_in_formats(figure, OUTPUT / name, formats=[".png", ".pdf"])
    plt.close(figure)


def grid(table, family, turbine, column):
    selected = table.loc[table.family.eq(family) & np.isclose(table.turbine_efficiency, turbine)].copy()
    selected.loc[~selected.success, column] = np.nan
    return selected.pivot(index="pump_efficiency", columns="exploitation_pct", values=column).reindex(
        index=sorted(run.PUMP_EFFICIENCIES), columns=sorted(run.EXPLOITATION_PCT)).to_numpy(float)


def contour(figure, ax, values, levels, title, label, cmap="RdYlGn", norm=None, signed=False):
    x, y = sorted(run.EXPLOITATION_PCT), 100*np.array(sorted(run.PUMP_EFFICIENCIES))
    values = np.ma.masked_invalid(values)
    ax.set_facecolor("0.88")
    if values.count():
        filled = ax.contourf(x, y, values, levels=levels, cmap=cmap, norm=norm, corner_mask=False)
        valid_levels = [v for v in levels if values.min() < v < values.max()]
        if valid_levels:
            lines = ax.contour(x, y, values, levels=valid_levels, colors="0.22",
                               linewidths=[1.4 if signed and v == 0 else 0.5 for v in valid_levels],
                               linestyles="solid", corner_mask=False)
            ax.clabel(lines, fmt=lambda v: f"{v:+g}%" if signed and v else f"{v:g}%", fontsize=8)
        figure.colorbar(filled, ax=ax, label=label, shrink=0.9)
    else:
        ax.text(0.5, 0.5, "No accepted cases", transform=ax.transAxes, ha="center")
    ax.set(xlabel="Allowed heat-source exploitation [%]", ylabel="Working-fluid pump efficiency [%]",
           xlim=(0, 100), ylim=(20, 100), title=title)
    ax.set_xticks(np.arange(0, 101, 20))
    ax.set_yticks(np.arange(20, 101, 10))
    ax.grid(False)
    top = ax.secondary_xaxis("top", functions=(lambda e: 140-1.2*e, lambda t: (140-t)/1.2))
    top.set_xlabel("Minimum allowed source outlet temperature [°C]")
    top.set_xticks([140, 116, 92, 68, 44, 20])
    ax.title.set_position((0.5, 1.18))


def maps(table):
    efficiency = [grid(table, family, 0.9, "fixed_reference_efficiency_pct") for family in run.FLUIDS]
    finite = np.concatenate([g[np.isfinite(g)] for g in efficiency])
    levels = MaxNLocator(nbins=15).tick_values(0, max(finite.max(), 1)) if finite.size else np.linspace(0, 1, 11)
    for family, values in zip(run.FLUIDS, efficiency):
        figure, ax = plt.subplots(figsize=(9, 7), layout="constrained")
        contour(figure, ax, values, levels, f"{family.capitalize()}, expander efficiency 90%",
                "Net system efficiency on fixed 140 to 20°C heat basis [%]")
        figure.supxlabel("Gray: unavailable result; 0% exploitation has no available heat.", fontsize=9)
        save(figure, f"system_efficiency_{family}_turb90")
    gains = [grid(table, "butane", eta, "relative_gain_pct") for eta in run.BUTANE_TURBINE_EFFICIENCIES]
    finite = np.concatenate([g[np.isfinite(g)] for g in gains])
    low = min(-5, 5*np.floor(finite.min()/5)) if finite.size else -5
    high = max(5, 5*np.ceil(finite.max()/5)) if finite.size else 5
    levels = np.unique(np.r_[MaxNLocator(nbins=20).tick_values(low, high), 0])
    # One scale across all four panels, with the sign boundary at yellow.
    norm = TwoSlopeNorm(vmin=low, vcenter=0, vmax=high)
    for eta, values in zip(run.BUTANE_TURBINE_EFFICIENCIES, gains):
        figure, ax = plt.subplots(figsize=(9, 7), layout="constrained")
        contour(figure, ax, values, levels, f"Butane expander {100*eta:g}% vs propane expander 90%",
                "Relative net system efficiency change [%]", norm=norm, signed=True)
        boundaries = []
        x, y = sorted(run.EXPLOITATION_PCT), 100*np.array(sorted(run.PUMP_EFFICIENCIES))
        for field, family, efficiency, level, color, style, label in (
            ("high_pressure_bar", "propane", 0.9, float(run.fluid("propane").critical_point.p)/1e5,
             "#2166ac", "--", "Propane critical pressure"),
            ("inlet_quality", "butane", eta, 1.0, "0.15", ":", "Butane saturated vapor inlet"),
        ):
            data = np.ma.masked_invalid(grid(table, family, efficiency, field))
            if data.count() and data.min() < level < data.max():
                ax.contour(x, y, data, levels=[level], colors=[color], linestyles=[style],
                           linewidths=1.2, corner_mask=False)
                boundaries.append(Line2D([], [], color=color, linestyle=style, linewidth=1.2, label=label))
        if boundaries:
            ax.legend(handles=boundaries, loc="upper right", fontsize=8, framealpha=0.92)
        figure.supxlabel("100 × (butane efficiency / propane efficiency − 1); positive favors butane.\n"
                            "Gray: unavailable result; comparison uses a common heat reference.", fontsize=9)
        save(figure, f"relative_gain_butane_turb{round(100*eta):02d}")
    for column, label, name in (
        ("high_pressure_bar", "Optimized expander inlet pressure [bar]", "optimized_pressure"),
        ("backwork_ratio_pct", "Working-fluid pump power / expander power [%]", "backwork_ratio"),
        ("actual_exploitation_pct", "Actual heat-source exploitation [%]", "actual_exploitation"),
    ):
        values = [grid(table, family, 0.9, column) for family in run.FLUIDS]
        finite = np.concatenate([g[np.isfinite(g)] for g in values])
        levels = MaxNLocator(nbins=12).tick_values(0, max(finite.max(), 1)) if finite.size else np.linspace(0, 1, 11)
        figure, axes = plt.subplots(1, 2, figsize=(16, 6.8), layout="constrained")
        for ax, family, data in zip(axes, run.FLUIDS, values):
            contour(figure, ax, data, levels, f"{family.capitalize()}, expander efficiency 90%", label)
        save(figure, name)


def cycles(table):
    finite_exploitation = sorted(e for e in run.EXPLOITATION_PCT if e > 0)
    if not finite_exploitation:
        return
    cmap = LinearSegmentedColormap.from_list("pump_magma", plt.get_cmap("magma")(np.linspace(0.25, 0.75, 256)))
    norm = Normalize(vmin=min(run.PUMP_EFFICIENCIES)*100, vmax=max(run.PUMP_EFFICIENCIES)*100)
    styles = {"working_fluid": "-", "heating_fluid": "--", "cooling_fluid": ":"}
    figures, ranges = [], {family: {"s": [], "heat": []} for family in run.FLUIDS}
    for exploitation in (finite_exploitation[0], finite_exploitation[-1]):
        figure, axes = plt.subplots(2, 2, figsize=(13.5, 10.5))
        selected = table.loc[table.success & table.exploitation_pct.eq(exploitation)
                             & np.isclose(table.turbine_efficiency, 0.9)]
        for column, family in enumerate(run.FLUIDS):
            rows = selected.loc[selected.family.eq(family)].sort_values("pump_efficiency", ascending=False)
            for index, row in enumerate(rows.itertuples(index=False)):
                directory = run.RESULTS / row.case
                config = run.th.read_configuration_file(directory / "optimal_solution.yaml")
                problem = run.ComparisonProblem(config, out_dir=str(directory))
                problem.fitness(problem.x0)
                # Reconstruct the assumed expansion trajectory for the T-s curve.
                p, h = run.expansion_path(problem, count=101)
                problem.cycle_data["components"]["expander"]["states"] = problem.cycle_data["working_fluid"].get_state(run.props.HmassP_INPUTS, h, p)
                problem.plot_settings["fluid"]["plot_saturation_line"] = index == 0
                problem.plot_settings["fluid"]["plot_critical_point"] = index == 0
                problem.figure, problem.axes = figure, axes[:, column]
                problem.plot_cycle()
                color = cmap(norm(100*row.pump_efficiency))
                components = problem.cycle_data["components"]
                for name, line in problem.graphics["process_lines"][0].items():
                    if name.endswith("_hot_side"):
                        stream = components[name[:-9]]["hot_side"]
                    elif name.endswith("_cold_side"):
                        stream = components[name[:-10]]["cold_side"]
                    else:
                        stream = components[name]
                    line.set(color=color, linewidth=0.9, linestyle=styles[stream["states"]["identifier"]], alpha=0.95)
                    ranges[family]["s"].extend(np.asarray(line.get_xdata(), float))
                for points in problem.graphics["state_points"][0].values():
                    points.set_visible(False)
                for name, artists in problem.graphics["pinch_point_lines"][1].items():
                    for side, artist_name in (("hot_side", "hot_line"), ("cold_side", "cold_line")):
                        stream = components[name][side]
                        line = artists[artist_name]
                        line.set(color=color, linewidth=0.9, linestyle=styles[stream["states"]["identifier"]], alpha=0.95)
                        ranges[family]["heat"].extend(np.asarray(line.get_xdata(), float))
                    for artist_name, artist in artists.items():
                        if artist_name not in ("hot_line", "cold_line"):
                            artist.set_visible(False)
            axes[0, column].set_title(f"{family.capitalize()} T–s, expander 90%")
            axes[1, column].set_title(f"{family.capitalize()} T–Q, expander 90%")
            axes[0, column].set_xlabel("Working-fluid entropy [kJ/(kg K)]")
            axes[0, column].xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v/1000:g}"))
            axes[1, column].set_xlabel("Cumulative heat flow rate [MW]")
            axes[1, column].xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v/1e6:g}"))
            for ax in axes[:, column]:
                ax.set_ylabel("Temperature [°C]")
                ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v-273.15:g}"))
                ax.set_yticks(np.arange(20, 161, 20)+273.15)
                ax.set_ylim(288.15, 428.15)
                ax.grid(False)
        figure.suptitle(f"Allowed exploitation {exploitation}% · minimum source outlet {140-1.2*exploitation:g}°C", y=0.98)
        figure.legend(handles=[Line2D([], [], color="0.25", linestyle=styles[stream], label=label)
                       for stream, label in (("working_fluid", "Working fluid"), ("heating_fluid", "Heat source"),
                                             ("cooling_fluid", "Cooling water"))],
                       loc="lower center", bbox_to_anchor=(0.46, 0.02), ncol=3, frameon=False)
        figure.subplots_adjust(left=0.07, right=0.86, bottom=0.13, top=0.90, hspace=0.32, wspace=0.28)
        colorbar = figure.colorbar(ScalarMappable(norm=norm, cmap=cmap),
                                  cax=figure.add_axes((0.90, 0.15, 0.018, 0.73)),
                                  ticks=100*np.array(sorted(run.PUMP_EFFICIENCIES)))
        colorbar.set_label("Working-fluid pump efficiency [%]")
        figures.append((figure, axes, exploitation))
    for figure, axes, exploitation in figures:
        for column, family in enumerate(run.FLUIDS):
            entropy = ranges[family]["s"]
            heat = ranges[family]["heat"]
            if entropy:
                span = max(entropy)-min(entropy)
                axes[0, column].set_xlim(min(entropy)-0.03*span, max(entropy)+0.03*span)
            if heat:
                axes[1, column].set_xlim(-0.01*max(heat), 1.03*max(heat))
        save(figure, f"cycles_exploitation_{exploitation:03d}")


def plot_all(table):
    maps(table)
    cycles(table)


if __name__ == "__main__":
    plot_all(pd.read_csv(run.RESULTS / "comparison.csv"))
