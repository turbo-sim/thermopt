"""Write the study discussion and comparison tables from saved sweep results."""

import run_sensitivity as run

np, pd = run.np, run.pd


def markdown_table(frame, columns, formats=None):
    formats = formats or {}
    rows = ["| " + " | ".join(columns.values()) + " |",
            "| " + " | ".join("---" for _ in columns) + " |"]
    for row in frame.to_dict("records"):
        cells = []
        for key in columns:
            value = row[key]
            if pd.isna(value):
                cells.append("Unavailable")
            elif key in formats:
                cells.append(formats[key](value))
            elif isinstance(value, (float, np.floating)):
                cells.append(f"{value:.2f}")
            else:
                cells.append(str(value))
        rows.append("| " + " | ".join(cells) + " |")
    return "\n".join(rows)


def make_report(table):
    accepted = table.loc[table.success & table.exploitation_pct.gt(0)].copy()
    coverage = table.groupby(["family", "turbine_efficiency"], as_index=False).agg(
        accepted=("success", "sum"), total=("success", "count"),
        starts=("attempts", "sum"), seconds=("total_solver_seconds", "sum"))
    coverage.to_csv(run.RESULTS / "coverage.csv", index=False)
    deep = accepted.loc[accepted.exploitation_pct.eq(max(run.EXPLOITATION_PCT))]
    diagnostic = deep.loc[np.isclose(deep.turbine_efficiency, 0.9)
                          & deep.pump_efficiency.isin([0.2, 0.5, 0.7, 1.0])].sort_values(["pump_efficiency", "family"])
    diagnostic.to_csv(run.RESULTS / "deep_recovery_diagnostics.csv", index=False)
    diagnostic_display = diagnostic.copy()
    for location in ("inlet", "outlet"):
        diagnostic_display[f"{location}_state"] = diagnostic.apply(
            lambda row: f"{row[f'{location}_quality']:.2f} vapor fraction" if row[f"{location}_phase"] == "two_phase"
            else row[f"{location}_phase"].replace("_", " "), axis=1)
    gains = table.loc[table.family.eq("butane") & table.exploitation_pct.eq(max(run.EXPLOITATION_PCT))].pivot(
        index="pump_efficiency", columns="turbine_efficiency", values="relative_gain_pct").sort_index()
    gains.columns = [f"gain_turb{round(100*eta):02d}" for eta in gains.columns]
    gains = gains.reset_index()
    gains.to_csv(run.RESULTS / "deep_recovery_gains.csv", index=False)
    expander_break_even = []
    for pump in run.PUMP_EFFICIENCIES:
        selected = deep.loc[deep.family.eq("butane") & np.isclose(deep.pump_efficiency, pump)
                            & deep.relative_gain_pct.notna()].sort_values("turbine_efficiency")
        for lower, upper in zip(selected.to_dict("records"), selected.to_dict("records")[1:]):
            if lower["relative_gain_pct"] <= 0 <= upper["relative_gain_pct"]:
                fraction = -lower["relative_gain_pct"]/(upper["relative_gain_pct"]-lower["relative_gain_pct"])
                expander_break_even.append(dict(pump_efficiency=pump,
                    interpolated_break_even_expander_efficiency_pct=100*(lower["turbine_efficiency"]+
                        fraction*(upper["turbine_efficiency"]-lower["turbine_efficiency"])),
                    lower_sampled_expander_efficiency_pct=100*lower["turbine_efficiency"],
                    upper_sampled_expander_efficiency_pct=100*upper["turbine_efficiency"]))
    expander_break_even = pd.DataFrame(expander_break_even, columns=["pump_efficiency",
        "interpolated_break_even_expander_efficiency_pct", "lower_sampled_expander_efficiency_pct",
        "upper_sampled_expander_efficiency_pct"])
    expander_break_even.to_csv(run.RESULTS / "deep_recovery_expander_break_even.csv", index=False)
    threshold_rows = []
    for eta in run.BUTANE_TURBINE_EFFICIENCIES:
        for pump in run.PUMP_EFFICIENCIES:
            selected = table.loc[table.family.eq("butane") & np.isclose(table.turbine_efficiency, eta)
                                 & np.isclose(table.pump_efficiency, pump) & table.exploitation_pct.gt(0)].sort_values("exploitation_pct")
            if selected.empty:
                continue
            positive = selected.loc[selected.relative_gain_pct.gt(0)]
            threshold_rows.append(dict(turbine_efficiency=eta, pump_efficiency=pump,
                minimum_sampled_exploitation_for_gain_pct=positive.exploitation_pct.min() if len(positive) else np.nan,
                maximum_gain_pct=selected.relative_gain_pct.max(),
                gain_at_maximum_exploitation_pct=selected.iloc[-1].relative_gain_pct))
    thresholds = pd.DataFrame(threshold_rows)
    thresholds.to_csv(run.RESULTS / "sampled_break_even.csv", index=False)
    irregularities = []
    for keys, group in accepted.groupby(["family", "turbine_efficiency", "pump_efficiency"]):
        group = group.sort_values("exploitation_pct")
        previous = None
        for row in group.itertuples():
            if previous is not None:
                change = 100*(row.fixed_reference_efficiency_pct/previous.fixed_reference_efficiency_pct-1)
                if change < -0.01:
                    irregularities.append(dict(family=keys[0], turbine_efficiency=keys[1], pump_efficiency=keys[2],
                        previous_exploitation_pct=previous.exploitation_pct, exploitation_pct=row.exploitation_pct,
                        relative_efficiency_drop_pct=change))
            previous = row
    pd.DataFrame(irregularities, columns=["family", "turbine_efficiency", "pump_efficiency",
        "previous_exploitation_pct", "exploitation_pct", "relative_efficiency_drop_pct"]).to_csv(run.RESULTS / "local_solution_irregularities.csv", index=False)
    phase_counts = accepted.groupby(["family", "turbine_efficiency", "pressure_regime", "inlet_phase", "outlet_phase"]).size().reset_index(name="cases")
    phase_counts.to_csv(run.RESULTS / "phase_counts.csv", index=False)
    architecture_rows = []
    for eta in run.BUTANE_TURBINE_EFFICIENCIES:
        selected = accepted.loc[accepted.family.eq("butane") & np.isclose(accepted.turbine_efficiency, eta)
            & accepted.inlet_phase.isin(["two_phase", "liquid"]) & accepted.propane_pressure_regime.eq("transcritical")]
        architecture_rows.append(dict(turbine_efficiency=eta, paired_conditions=len(selected),
            positive_gains=int(selected.relative_gain_pct.gt(0).sum()),
            minimum_gain_pct=selected.relative_gain_pct.min(), maximum_gain_pct=selected.relative_gain_pct.max()))
    architectures = pd.DataFrame(architecture_rows)
    architectures.to_csv(run.RESULTS / "two_phase_vs_transcritical.csv", index=False)
    conclusions, recovery = [], ""
    middle = deep.loc[np.isclose(deep.pump_efficiency, 0.7) & np.isclose(deep.turbine_efficiency, 0.9)]
    if set(middle.family) == {"propane", "butane"}:
        propane = middle.loc[middle.family.eq("propane")].iloc[0]
        butane = middle.loc[middle.family.eq("butane")].iloc[0]
        conclusions.append(
            f"With maximum allowed exploitation and a 70%-efficient pump, butane at 90% expander efficiency reaches **{butane.fixed_reference_efficiency_pct:.2f}%** fixed-reference system efficiency versus **{propane.fixed_reference_efficiency_pct:.2f}%** for propane, a **{butane.relative_gain_pct:+.2f}%** relative gain. "
            f"Butane operates at **{butane.high_pressure_bar:.2f} bar** with a **{butane.backwork_ratio_pct:.2f}%** back-work ratio, versus **{propane.high_pressure_bar:.2f} bar** and **{propane.backwork_ratio_pct:.2f}%** for propane.")
        recovered_ratio = (butane.fixed_reference_efficiency_pct/butane.cycle_efficiency_pct)/(propane.fixed_reference_efficiency_pct/propane.cycle_efficiency_pct)
        conversion_ratio = butane.cycle_efficiency_pct/propane.cycle_efficiency_pct
        recovery = (f"At that 70%-pump, 90%-expander point, butane cools the source to **{butane.source_exit_C:.2f}°C** while propane leaves it at **{propane.source_exit_C:.2f}°C**. "
            f"Recovered heat per unit source flow changes by **{100*(recovered_ratio-1):+.2f}%**, and net power per unit absorbed heat changes by **{100*(conversion_ratio-1):+.2f}%**. "
            f"These factors multiply to give the **{butane.relative_gain_pct:+.2f}%** system-power gain. Thus the benefit combines source recovery and conversion, and can persist even when conversion efficiency by itself is slightly lower.")
    for family in run.FLUIDS:
        full = deep.loc[deep.family.eq(family) & np.isclose(deep.turbine_efficiency, 0.9)]
        low, high = full.loc[np.isclose(full.pump_efficiency, 0.2)], full.loc[np.isclose(full.pump_efficiency, 1)]
        if len(low) and len(high):
            decline = 100*(1-low.iloc[0].fixed_reference_efficiency_pct/high.iloc[0].fixed_reference_efficiency_pct)
            conclusions.append(f"For {family}, reducing pump efficiency from 100% to 20% at maximum allowed exploitation reduces optimized system efficiency by **{decline:.2f}%** at 90% expander efficiency.")
    propane_deep = deep.loc[deep.family.eq("propane")]
    subcritical_pumps = propane_deep.loc[propane_deep.pressure_regime.eq("subcritical"), "pump_efficiency"].sort_values()
    transcritical_pumps = propane_deep.loc[propane_deep.pressure_regime.eq("transcritical"), "pump_efficiency"].sort_values()
    pressure_finding = ("At maximum allowed exploitation, propane selects subcritical pressure at pump efficiencies "
        + ", ".join(f"{100*value:.0f}%" for value in subcritical_pumps)
        + " and transcritical pressure at " + ", ".join(f"{100*value:.0f}%" for value in transcritical_pumps) + ". "
        + "This transition is an optimization result within the permitted pressure bounds. Poor pump efficiency can change the preferred pressure regime as well as reducing efficiency.")
    figures = "results/figures"
    lines = ["# Propane and butane cycle comparison at a 140 degree source", "",
        "Simple propane cycles with dry expansion and simple n-butane cycles with a free expander inlet state are optimized to compare the benefit of source matching against working-fluid pump losses and prescribed expander efficiency. Positive relative gains favor butane at the same allowed source cooling and pump efficiency.", "",
        "\n\n".join(conclusions), "",
        "## Boundary conditions and comparison basis", "",
        "The source is water at 140°C and 10 bar. Cooling water enters at 20°C and exits at 24.99–25.00°C. Both cycles have no recuperator, 1% pressure loss on each exchanger side, at least 5 K approach at 100 uniform exchanger locations, and at least 1 K pump-inlet subcooling. Every design is sized to 1 MW after the working-fluid pump, before the external-water pumps. Those external pumps retain 70% efficiency.", "",
        "Allowed exploitation is `(140 − minimum source outlet temperature) / 120`, sampled in 5 percentage-point steps. The source outlet is optimized above this minimum, so 100% permission does not require cooling to 20°C. Working-fluid pump efficiency is sampled from 20% to 100% in 10 percentage-point steps. Propane's expander efficiency is 90%; butane's is 90%, 80%, 70%, or 60%.", "",
        "The free variables are source outlet temperature, cooling-water outlet temperature, working-fluid pump inlet pressure and enthalpy, and expander inlet pressure and enthalpy. High-pressure bounds extend to three times the critical pressure for both fluids and permit subcritical or supercritical heating. No wet inlet is imposed on butane. Propane has nonnegative exit superheat and 17 additional vapor-enthalpy constraints along the saturation-capable part of a cumulative constant-isentropic-efficiency expansion path. Accepted propane designs also pass a 201-point check of that path.", "",
        f"The figures use `η_fixed = net system power / [source mass flow × Δh_water(140→20°C at 10 bar)]`, with **Δh = {run.reference_heat()/1000:.3f} kJ/kg**. Net system power subtracts all three pumps. This measures electricity per unit source flow on a common heat reference. The native system efficiency, divided by heat available down to the allowed source outlet, remains the optimization objective and is exported separately. At a given grid point the denominators differ by a positive constant, so both objectives select the same optimum.", "",
        "Relative change is `100 × (η_fixed,butane / η_fixed,propane − 1)`. Since the denominator is common, this is also the net power change at fixed source mass flow. Percentage-point differences are exported separately. At 0% exploitation no heat can be released; these cases are recorded as unavailable and receive no optimization.", "",
        "## Optimization and numerical checks", "",
        "Each finite point receives a continuation or high-pressure start and an independent vapor-oriented start. A third, perturbed start is used when neither succeeds or accepted starts differ by more than 0.2%. The maximum is three starts per point. The best feasible converged result is retained; all attempts are logged. These limited restarts reduce sensitivity to initialization but do not establish a global optimum.", "",
        "The optimization evaluator computes the same simple-cycle heat balances and uniform exchanger constraints using low-level CoolProp. Native ThermoOpt reevaluates every completed candidate before acceptance and supplies all reported efficiencies and cycle diagrams. Saturation departures used in constraints are recomputed from returned pressure, enthalpy, and temperature because some jaxprop versions perturb their internal state during wet-state derivative calculations. No repository core code is changed.", "",
        markdown_table(coverage, {"family": "Fluid", "turbine_efficiency": "Expander efficiency", "accepted": "Accepted", "total": "Requested including zero heat", "starts": "Optimization calls"},
                       {"turbine_efficiency": lambda v: f"{100*v:.0f}%"}), "",
        f"The accepted solutions have maximum normalized constraint violation **{accepted.max_constraint_violation.max():.3g}** and maximum energy-balance residual **{accepted.energy_balance_W.abs().max():.3g} W**. Failed cases are masked in every plot and comparison.", "",
        "## Cycle pressure and phase selected by optimization", "",
        markdown_table(phase_counts, {"family": "Fluid", "turbine_efficiency": "Expander efficiency", "pressure_regime": "Pressure regime", "inlet_phase": "Inlet phase", "outlet_phase": "Outlet phase", "cases": "Cases"},
                       {"turbine_efficiency": lambda v: f"{100*v:.0f}%"}), "",
        pressure_finding, "",
        "Liquid-admission butane solutions approach the trilateral configuration. The assumed liquid or wet expander performance is still an input: the optimizer establishes the thermodynamic state and its required work, rather than demonstrating that a practical expander attains the prescribed efficiency.", "",
        "A subcritical propane result is part of the permitted comparison. Likewise a vapor-inlet butane result is a conventional Rankine-like cycle within the flexible formulation. The assumed butane expander efficiency applies to its entire sweep, including vapor-inlet points. These points should not be described as demonstrated advantages of two-phase admission.", "",
        "## Cycle diagrams at low and high allowed exploitation", "",
        f"![Cycles at minimum finite exploitation]({figures}/cycles_exploitation_{min(e for e in run.EXPLOITATION_PCT if e>0):03d}.png)", "",
        f"![Cycles at maximum allowed exploitation]({figures}/cycles_exploitation_{max(run.EXPLOITATION_PCT):03d}.png)", "",
        "Both figures show the 90%-efficient expander slices and all accepted pump efficiencies. Magma colors encode pump efficiency; solid, dashed, and dotted lines indicate working fluid, source, and cooling water. The expansion curves follow the assumed cumulative-efficiency path. External-stream traces in T–s use the associated working-fluid entropy coordinate, following the native plot convention. T–Q uses the native cumulative duty coordinate and the common 1 MW net-cycle sizing. Heat duties therefore change with the flow needed by each design.", "",
        "## Net system efficiency", "",
        f"![Propane system efficiency]({figures}/system_efficiency_propane_turb90.png)", "",
        f"![Butane system efficiency]({figures}/system_efficiency_butane_turb90.png)", "",
        "The two efficiency maps share their color levels. Actual exploitation can plateau below the allowed limit because additional cooling becomes less valuable than maintaining a favorable heat-addition profile and expansion work.", "",
        "## Relative change as butane expander efficiency decreases", ""]
    for eta in run.BUTANE_TURBINE_EFFICIENCIES:
        selected = accepted.loc[accepted.family.eq("butane") & np.isclose(accepted.turbine_efficiency, eta)]
        wet = selected.loc[selected.inlet_phase.eq("two_phase")]
        finite = selected.loc[selected.relative_gain_pct.notna()]
        positive = finite.loc[finite.relative_gain_pct.gt(0)]
        lines.extend([f"![Relative gain at butane expander {100*eta:.0f}%]({figures}/relative_gain_butane_turb{round(100*eta):02d}.png)", "",
            f"At **{100*eta:.0f}% butane expander efficiency**, butane wins at **{len(positive)} of {len(finite)}** accepted paired conditions. Relative changes range from **{finite.relative_gain_pct.min():+.2f}% to {finite.relative_gain_pct.max():+.2f}%**. Of {len(selected)} accepted butane designs, {len(wet)} have two-phase admission.", ""])
    gain_columns = {"pump_efficiency": "Pump efficiency"}
    gain_columns.update({f"gain_turb{round(100*eta):02d}": f"Butane expander {100*eta:.0f}%" for eta in run.BUTANE_TURBINE_EFFICIENCIES})
    lines.extend(["At maximum allowed exploitation, relative efficiency changes are:", "",
        markdown_table(gains, gain_columns, {"pump_efficiency": lambda v: f"{100*v:.0f}%",
            **{key: lambda v: f"{v:+.2f}%" for key in gain_columns if key != "pump_efficiency"}}), "",
        "The zero contour gives an interpolated break-even boundary. `sampled_break_even.csv` records the first sampled exploitation with positive gain for each pump and butane-expander efficiency. The 5% exploitation and 10% pump grids limit the precision of any boundary.", "",
        "The blue dashed line marks propane's critical-pressure boundary; the black dotted line marks the butane saturated-vapor inlet boundary. Both are interpolated from optimized grid values. Gains elsewhere can compare subcritical propane or vapor-admission butane. The subset that directly compares liquid or wet-admission butane with transcritical propane is:", "",
        markdown_table(architectures, {"turbine_efficiency": "Butane expander efficiency", "paired_conditions": "Paired conditions",
            "positive_gains": "Butane wins", "minimum_gain_pct": "Minimum gain [%]", "maximum_gain_pct": "Maximum gain [%]"},
            {"turbine_efficiency": lambda v: f"{100*v:.0f}%", "minimum_gain_pct": lambda v: f"{v:+.2f}",
             "maximum_gain_pct": lambda v: f"{v:+.2f}"}), "",
        "In this direct architecture comparison, the 70% and 60% butane expanders do not beat the 90% propane expander at any sampled condition. Positive gains in their full fluid maps must be read with the phase and pressure boundaries. In particular, the 70% butane expander is essentially tied at 20% pump efficiency and maximum allowed exploitation, where propane selects a subcritical cycle; that near-zero difference is smaller than the observed local-solution irregularities.", "",
        "Linear interpolation between neighboring prescribed butane expander efficiencies gives the following approximate break-even values at maximum allowed exploitation. Each fluid is reoptimized at the sampled efficiencies, so these estimates describe a boundary between optimized designs rather than a fixed-design expander calculation.", "",
        markdown_table(expander_break_even, {"pump_efficiency": "Pump efficiency",
            "interpolated_break_even_expander_efficiency_pct": "Approximate butane expander break even [%]"},
            {"pump_efficiency": lambda v: f"{100*v:.0f}%"}), "",
        "## Pressure levels and back work explain the tradeoff", "",
        markdown_table(diagnostic_display, {"family": "Fluid", "pump_efficiency": "Pump efficiency", "fixed_reference_efficiency_pct": "Efficiency [%]", "high_pressure_bar": "High pressure [bar]", "backwork_ratio_pct": "Back work [%]", "source_exit_C": "Source outlet [°C]", "inlet_state": "Inlet state", "outlet_state": "Outlet state"},
                       {"pump_efficiency": lambda v: f"{100*v:.0f}%"}), "",
        f"![Optimized pressures]({figures}/optimized_pressure.png)", "",
        f"![Backwork ratio]({figures}/backwork_ratio.png)", "",
        f"![Actual exploitation]({figures}/actual_exploitation.png)", "",
        "The comparison includes changes in both heat recovery and conversion efficiency. Back-work ratio is the working-fluid pump power divided by expander power. Absolute pressure, pressure ratio, inlet density, and specific expansion work all contribute; lower absolute pressure alone does not guarantee a smaller back-work ratio. The optimizer can respond to poor pump efficiency by changing pressure and inlet state, so a fixed-design pump-loss estimate does not predict the optimized transition by itself.", "",
        recovery, "",
        "## Local optimization and interpretation limits", ""])
    if irregularities:
        maximum = max(-row["relative_efficiency_drop_pct"] for row in irregularities)
        lines.extend([f"There are **{len(irregularities)}** neighboring exploitation pairs where fixed-reference efficiency falls by more than 0.01% despite allowing more source cooling; the largest decrease is **{maximum:.3f}%**. These are local-solution irregularities because relaxing the outlet-temperature limit retains every earlier design's feasibility. They are preserved in `local_solution_irregularities.csv` and should not be interpreted as a thermodynamic penalty for allowing deeper cooling.", ""])
    else:
        lines.extend(["No accepted neighboring exploitation points show an efficiency decrease exceeding 0.01% as more source cooling is allowed.", ""])
    lines.extend(["This is a thermodynamic comparison of two fluids and prescribed component efficiencies. It does not isolate fluid choice from cycle architecture or predict the efficiency, cost, leakage, volume ratio, or lifetime of a real wet expander. The triangular heat-addition profile can improve recovery, but the required expander performance must be achievable at the optimized inlet and outlet qualities. Propane's dry-path requirement is tied to the stated expansion-path model. Exchanger approaches are sampled constraints, so small gains near zero remain sensitive to discretization and local optimization.", ""])
    (run.HERE / "report.md").write_text("\n".join(lines), encoding="utf-8")
    return coverage, thresholds, phase_counts


if __name__ == "__main__":
    make_report(pd.read_csv(run.RESULTS / "comparison.csv"))
