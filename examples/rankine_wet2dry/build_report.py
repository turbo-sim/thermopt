"""Build the standalone report from completed saved studies; no optimizations."""

from pathlib import Path
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
FAMILIES = ("butane", "cyclopentane", "novec")
LABELS = {"butane": "Butane: simple cycle", "cyclopentane": "Cyclopentane: simple cycle",
          "novec": "Novec 649: recuperated cycle"}
TABLE_TEMPERATURES_C = (80, 70, 60, 50, 40, 20)


def row_at(table, mode, temperature, efficiency=0.9):
    subset = table.loc[table["mode"].eq(mode) & table["source_min_C"].eq(temperature)
                       & np.isclose(table["turbine_efficiency"], efficiency)]
    if len(subset) != 1 or not bool(subset.iloc[0]["success"]):
        raise ValueError(f"Missing accepted result: {mode}/{temperature}/{efficiency}")
    return subset.iloc[0]


def figure(family, stem, description):
    path = f"results/sensitivity_{family}/{stem}.png"
    if not (HERE / path).exists():
        raise FileNotFoundError(path)
    return f"![{description}]({path})\n"


def value_range(values):
    lower, upper = f"{values.min():.2f}", f"{values.max():.2f}"
    return lower if lower == upper else f"{lower}–{upper}"


def fluid_section(family, number, table, audit):
    title = LABELS[family]
    partial = table.loc[table["mode"].eq("partial_evaporation") & table["success"]]
    checks = pd.read_csv(HERE / "results" / f"sensitivity_{family}" / "saved_solution_checks.csv").set_index("case")
    b135, b70, b20 = [row_at(table, "superheated", t) for t in (135, 70, 20)]
    p135, p70, p20 = [row_at(table, "partial_evaporation", t) for t in (135, 70, 20)]
    cold = partial.loc[partial["source_min_C"].eq(20)]
    classes = audit["expansion_class_counts"]["partial_evaporation"]
    baseline_state, partial_state = checks.loc[b70.case], checks.loc[p70.case]
    recovery = p20["achieved_heat_utilization_pct"] / b20["achieved_heat_utilization_pct"]
    conversion = p20["system_efficiency_pct"] / b20["system_efficiency_pct"] / recovery
    blocks = [f"## {number}. {title}\n", f"### {number}.1 Superheated baseline\n",
        figure(family, "baseline_performance", f"{title}: baseline net efficiency on a fixed heat reference and actual exploitation"),
        f"At 90% expander efficiency, the baseline's fixed-reference net efficiency changes from "
        f"**{b135.fixed_reference_efficiency_pct:.2f}%** at the 135 °C reinjection limit to "
        f"**{b20.fixed_reference_efficiency_pct:.2f}%** at 20 °C. Its corresponding net power per "
        f"unit source flow changes from **{b135.specific_system_power_kW_per_kg_s:.2f} to "
        f"{b20.specific_system_power_kW_per_kg_s:.2f} kW per kg/s**. These two measures have "
        "exactly the same trend because their heat reference is fixed.\n",
        f"At the 20 °C limit, the baseline actually leaves the source at **{b20.source_exit_C:.2f} °C**, "
        f"or **{b20.actual_exploitation_pct:.2f}% actual exploitation**. Permission for further "
        "cooling does not require the optimizer to use it: the selected pressure levels, phase-change "
        "temperatures and exchanger approaches determine useful recovery.\n",
        figure(family, "baseline_cycle_evolution", f"{title}: superheated T–s and T–Q evolution at 90% expander efficiency"),
        "The overlays show how the pressure levels and heat-addition profile change as the "
        "reinjection limit is relaxed. Coincident low-limit curves indicate that a similar design "
        "is selected despite allowing more source cooling.\n",
        f"### {number}.2 Evolution with partial evaporation\n",
        figure(family, "partial_evaporation_cycle_evolution", f"{title}: optimized-inlet T–s and T–Q evolution at 90% efficiency"),
        f"The optimized inlet has **{100*p70.inlet_quality_plot:.2f}% vapor** at a 70 °C limit "
        f"and **{100*p20.inlet_quality_plot:.2f}%** at 20 °C. Values below 100% shorten the "
        "evaporation portion of external heating. The balance between sensible heating and "
        "evaporation can then match the cooling source differently from the superheated cycle. "
        "The expansion outlet must also be examined; wet admission does not imply dry discharge.\n",
        f"At the 70 °C limit, the partial-evaporation design admits fluid at "
        f"**{partial_state.expander_inlet_pressure_bar:.2f} bar and {partial_state.expander_inlet_temperature_C:.2f} °C**, "
        f"compared with **{baseline_state.expander_inlet_pressure_bar:.2f} bar and "
        f"{baseline_state.expander_inlet_temperature_C:.2f} °C** for the baseline. "
        "Inlet-state flexibility therefore changes pressure level and temperature of heat addition, "
        "not just the amount of liquid entering the expander.\n",
        f"### {number}.3 Net efficiency on the common heat basis\n",
        figure(family, "system_efficiency_contours", f"{title}: fixed-reference net-efficiency contours"),
        f"At 90% expander efficiency, the plotted efficiency is **{p135.fixed_reference_efficiency_pct:.2f}%** "
        f"at a 135 °C limit, **{p70.fixed_reference_efficiency_pct:.2f}%** at 70 °C and "
        f"**{p20.fixed_reference_efficiency_pct:.2f}%** at 20 °C. Unlike efficiency divided by "
        "heat available only down to each reinjection limit, this figure directly tracks net "
        "electricity per unit geothermal flow across the whole horizontal axis. Any plateau "
        "means that permitting further source cooling adds little useful power.\n",
        f"### {number}.4 Actual source exploitation\n",
        figure(family, "actual_exploitation_contours", f"{title}: actual heat-source exploitation"),
        f"At the 20 °C limit, the four partial-evaporation designs actually leave the source at "
        f"**{value_range(cold.source_exit_C)} °C**, giving "
        f"**{value_range(cold.actual_exploitation_pct)}%** "
        "actual exploitation. The map distinguishes allowed cooling from cooling selected by "
        "the optimizer. More recovered heat alone does not establish greater power: its "
        "temperature and the expansion efficiency also matter.\n",
        f"### {number}.5 Expander inlet and outlet quality\n",
        figure(family, "inlet_quality_contours", f"{title}: expander inlet vapor quality in percent"),
        f"Across the 96 partial-formulation designs, **{sum(value for key,value in classes.items() if key.startswith('wet_'))}** "
        f"have two-phase inlets, **{sum(value for key,value in classes.items() if key.startswith('dry_'))}** "
        f"have dry-vapor inlets and **{sum(value for key,value in classes.items() if key.startswith('liquid_'))}** "
        "have liquid inlets. The formulation optimizes inlet state rather than imposing a wet inlet. "
        "A displayed value of 100% includes superheated vapor.\n",
        figure(family, "outlet_quality_contours", f"{title}: expander outlet vapor quality in percent"),
        "| Expansion endpoint phases | Number of designs |\n|---|---:|\n" +
        "\n".join(f"| {key.replace('_', ' ')} | {value} |" for key,value in classes.items()) + "\n",
        f"At 90% efficiency, displayed outlet quality is **{100*p70.outlet_quality_plot:.2f}%** "
        f"at the 70 °C limit and **{100*p20.outlet_quality_plot:.2f}%** at 20 °C. " +
        ("The deep-recovery design remains two-phase at discharge; its advantage therefore requires an expander suited to substantial liquid at both ends.\n"
         if p20.outlet_quality_plot < 1 - 1e-7 else
         "The deep-recovery design reaches dry vapor at discharge. If its inlet is wet, this is a wet-to-dry expansion; the dry outlet does not remove the need to accommodate liquid admission.\n"),
        f"### {number}.6 Relative gain and its physical origin\n",
        figure(family, "relative_gain_contours", f"{title}: relative net-power gain against the same fluid's 90%-efficient superheated baseline"),
        f"At the **70 °C** limit and 90% expander efficiency, partial evaporation changes "
        f"fixed-reference efficiency from **{b70.fixed_reference_efficiency_pct:.2f}% to "
        f"{p70.fixed_reference_efficiency_pct:.2f}%**, a **{p70.relative_gain_pct:+.2f}% relative gain**. "
        f"The actual source outlets are **{b70.source_exit_C:.2f} °C** for the baseline and "
        f"**{p70.source_exit_C:.2f} °C** for partial evaporation. " +
        ("Because both reach essentially the same outlet, this gain comes from converting the recovered heat more effectively, not from a larger recovered heat quantity.\n"
         if abs(b70.source_exit_C-p70.source_exit_C)<0.02 else
         "The differing source outlets show that heat recovery also contributes; the gain cannot be interpreted as a change in conversion efficiency alone.\n"),
        f"At the **20 °C** limit, partial evaporation at 90% efficiency changes recovered "
        f"source heat per kilogram by **{100*(recovery-1):+.2f}%** and net system power per "
        f"unit absorbed heat by **{100*(conversion-1):+.2f}%**. These factors multiply, "
        f"giving **{p20.relative_gain_pct:+.2f}% net-power gain** at equal geothermal flow. "
        "This separates the value of additional source cooling from the value of converting "
        "the recovered heat.\n",
        "| Minimum source outlet [°C] | Nominal exploitation [%] | Gain at expander efficiency 90% | 80% | 70% | 60% |\n"
        "|---:|---:|---:|---:|---:|---:|\n" + "\n".join(
            f"| {temperature} | {100*(140-temperature)/120:.2f} | " + " | ".join(
                f"{row_at(table,'partial_evaporation',temperature,eta).relative_gain_pct:+.2f}%"
                for eta in (.9,.8,.7,.6)) + " |" for temperature in TABLE_TEMPERATURES_C) + "\n",
        f"### {number}.7 Efficiency penalty and break-even conditions\n"]
    for efficiency in (.8,.7,.6):
        subset = partial.loc[np.isclose(partial.turbine_efficiency,efficiency)]
        positive = subset.loc[subset.relative_gain_pct > 0]
        best = subset.loc[subset.relative_gain_pct.idxmax()]
        text = (f"At **{100*efficiency:.0f}% expander efficiency**, the largest modeled gain is "
                f"**{best.relative_gain_pct:+.2f}%**. ")
        if positive.empty:
            text += "No sampled reinjection limit beats the 90%-efficient superheated baseline."
        else:
            onset = positive.loc[positive.source_min_C.idxmax()]
            text += (f"The highest sampled reinjection limit with positive gain is **{onset.source_min_C:g} °C**, "
                     f"where the gain is **{onset.relative_gain_pct:+.2f}%**.")
        if abs(best.relative_gain_pct) < 1:
            text += " This is a near-tie and does not establish a robust design advantage."
        blocks.append(text + "\n")
    blocks.append("The zero contour is an approximate break-even boundary between the sampled efficiencies. "
                  "A positive gain means the inlet-state and heat-recovery benefits offset the prescribed "
                  "expander penalty. It does not demonstrate that a real expander can attain that efficiency "
                  "at the computed inlet and outlet qualities. Near-zero gains and abrupt phase transitions "
                  "are especially sensitive to local optimization and exchanger discretization.\n")
    anomalies = audit["monotonic_specific_power_anomalies"]
    if anomalies:
        largest = max(anomalies, key=lambda item: item["relative_specific_power_loss_pct"])
        if largest["relative_specific_power_loss_pct"] > 0.01:
            blocks.append(
                f"The largest decrease in saved power per source flow relative to an earlier, stricter "
                f"outlet limit is **{largest['relative_specific_power_loss_pct']:.2f}%**. "
                "Relaxing the limit preserves the earlier design's feasibility, so such a decrease "
                "is a local-solution irregularity rather than a thermodynamic penalty for allowing "
                "more cooling. It is retained under the single-start protocol.\n")
    dominated = partial.loc[np.isclose(partial.turbine_efficiency, 0.9) & partial.relative_gain_pct.lt(-0.01)]
    if not dominated.empty:
        blocks.append(
            f"At equal 90% expander efficiency, the partial formulation falls below the "
            f"superheated baseline by more than 0.01% at {len(dominated)} sampled "
            f"outlet limit{'s' if len(dominated) != 1 else ''} (largest loss "
            f"**{-dominated.relative_gain_pct.min():.2f}%**). The flexible formulation includes "
            "the baseline's inlet state, so these negative values identify local-solution "
            "limitations, not an inherent disadvantage of allowing partial evaporation.\n")
    if family == "novec":
        epsilon = [checks.loc[row.case,"recuperator_effectiveness"] for row in (b70,p70,b20,p20)]
        blocks.insert(blocks.index(f"### {number}.3 Net efficiency on the common heat basis\n"),
                      f"For Novec, recuperator effectiveness is **{epsilon[0]:.3f}/{epsilon[1]:.3f}** "
                      f"for baseline/partial evaporation at 70 °C and **{epsilon[2]:.3f}/{epsilon[3]:.3f}** "
                      "at 20 °C. Internal preheating reduces external heat demand, but can also prevent "
                      "deep cooling of the geothermal stream. The optimizer balances these effects; "
                      "a recuperated architecture need not use high effectiveness at every condition.\n")
    return "\n".join(blocks)


def main():
    tables = {family: pd.read_csv(HERE / "results" / f"sensitivity_{family}" / "sweep_results.csv",
                                  float_precision="round_trip") for family in FAMILIES}
    audits = {family: json.loads((HERE / "results" / f"sensitivity_{family}" / "validation_summary.json")
                                 .read_text(encoding="utf-8")) for family in FAMILIES}
    comparison = pd.read_csv(HERE / "results/fluid_comparison/comparison.csv", float_precision="round_trip")
    assert len(comparison) == 6 and all(a["accepted_sweep_points"] == 120 for a in audits.values())
    reference_heat_kJ_per_kg = audits[FAMILIES[0]]["fixed_reference_heat_J_per_kg"] / 1000
    best_comparison = comparison.loc[comparison.fixed_reference_efficiency_pct.idxmax()]
    same_source_outlet = all(abs(row_at(tables[f], mode, 70).source_exit_C - 70) < .02
                            for f in FAMILIES for mode in ("superheated", "partial_evaporation"))
    sections = ["# Partial evaporation in geothermal Rankine cycles: butane, cyclopentane and Novec 649\n",
        "This study refines the evaporator and condenser to 100 elements to investigate when allowing liquid at an expander inlet increases geothermal "
        "electricity production enough to compensate for reduced expander efficiency. It compares "
        "simple n-butane and cyclopentane cycles with recuperated Novec 649, using each fluid's "
        "superheated cycle as the reference for its sensitivity study. A fixed heat reference lets "
        "efficiency contours be compared across every permitted source outlet temperature. Both "
        "inlet and outlet phase states are examined to distinguish wet-to-dry from wet-to-wet expansion.\n",
        "## 1. Boundary conditions and optimization formulation\n",
        "### 1.1 Common conditions\n",
        "| Parameter | Value |\n|---|---|\n"
        "| Source water | 140 °C, 10 bar |\n"
        "| Source minimum outlet | 70 °C for the six-case comparison; 140–20 °C in 5 K steps for sensitivities |\n"
        "| Cooling water | 20 °C inlet; 24.99–25.00 °C outlet |\n"
        "| Net cycle power used to size every design | 1 MW after the working-fluid pump, before external-water pumps |\n"
        "| Working fluids and architectures | Simple n-butane; simple cyclopentane; recuperated Novec 649 |\n"
        "| Expander efficiency | 90% for superheated baselines; 90%, 80%, 70%, 60% for partial evaporation |\n"
        "| Pump efficiencies | 70% |\n"
        "| Exchanger pressure losses | 1% per active side |\n"
        "| Minimum exchanger approach | 5 K at 100 evaporator/condenser locations; 25 recuperator locations |\n"
        "| Working-fluid pump inlet | At least 1 K subcooling |\n"
        "| Superheated formulation | 2 K inlet superheat and dry expander outlet |\n"
        "| Partial-evaporation formulation | Optimized inlet state; no fixed quality or dry-outlet constraint |\n",
        "### 1.2 Optimization variables\n",
        "| Variable | Role |\n|---|---|\n"
        "| Source outlet temperature | Amount of source cooling, subject to the reinjection limit |\n"
        "| Cooling-water outlet temperature | Approximately 5 K total warming |\n"
        "| Working-fluid pump inlet pressure and enthalpy | Low-pressure state and subcooling |\n"
        "| Expander inlet pressure and enthalpy | High-pressure state and inlet phase |\n"
        "| Recuperator effectiveness | 0–1, Novec only |\n",
        "Mass flows follow from the states, heat balances and 1 MW sizing target. Pressure and enthalpy "
        "bounds are broad and property-based. Each initial state comes from a saved coarse-grid optimum; "
        "the partial-evaporation formulation does not constrain the inlet quality. The cooling-water pump heats the water slightly "
        "before the condenser, so the condenser's own rise is slightly less than 5 K.\n",
        "### 1.3 Fixed-reference efficiency and comparison basis\n",
        "The efficiency contours and baseline performance curves use net system power after all pumps "
        "divided by the source heat available between **140 and 20 °C**, at the source pressure of 10 bar:\n\n"
        "```text\nηfixed = net system power / [source mass flow × (hwater(140 °C, 10 bar) − hwater(20 °C, 10 bar))]\n```\n",
        f"The reference enthalpy drop is **{reference_heat_kJ_per_kg:.3f} kJ/kg** for every case. This efficiency is directly "
        "proportional to net system power per source flow, so deeper permissible cooling is assessed "
        "on the same energy basis. The 1 MW cycle sizing does not imply equal geothermal flow.\n",
        "ThermoOpt's original system efficiency instead uses heat available down to the case's "
        "reinjection limit. It remains the optimization objective and is retained in the data. At "
        "a fixed limit the two denominators differ only by a constant, so maximizing either selects "
        "the same design. Relative efficiency gain at a matching limit equals net system power "
        "gain at equal geothermal flow, regardless of which denominator is used.\n",
        "## 2. Comparison of fluid and cycle configurations\n",
        "All six designs use a 70 °C reinjection limit and 90% expander efficiency. This table uses "
        "**superheated cyclopentane** as the common baseline; the sensitivities below use each fluid's own baseline.\n",
        "| Case | Inlet vapor [%] | Outlet vapor [%] | Source outlet [°C] | Gain vs cyclopentane baseline [%] | ThermoOpt system efficiency [%] | Fixed-reference net efficiency [%] |\n"
        "|---|---:|---:|---:|---:|---:|---:|\n" + "\n".join(
            f"| {r.case.replace('_',' ')} | {100*np.clip(r.inlet_quality,0,1):.2f} | "
            f"{100*np.clip(r.outlet_quality,0,1):.2f} | {r.source_exit_C:.2f} | {r.relative_gain_pct:+.2f} | "
            f"{r.system_efficiency_pct:.4f} | {r.fixed_reference_efficiency_pct:.4f} |"
            for r in comparison.itertuples()) + "\n",
        "Each fluid has a separate 2 × 2 comparison: superheated and partial-evaporation "
        "designs occupy the left and right columns, with T–s diagrams above T–Q diagrams.\n",
        "\n".join(
            f"![{LABELS[f]}: superheated and partial evaporation, T–s above T–Q]"
            f"(results/fluid_comparison/comparison_grid_{f}.png)\n"
            for f in ("cyclopentane", "butane", "novec")),
        "Relative to each fluid's own superheated design, the partial-evaporation gains at this limit are " +
        ", ".join(f"**{row_at(tables[f], 'partial_evaporation', 70).relative_gain_pct:+.2f}% for {LABELS[f].split(':')[0]}**"
                  for f in FAMILIES) + ". " +
        f"The highest fixed-reference efficiency is **{best_comparison.fixed_reference_efficiency_pct:.2f}%**, "
        f"for **{best_comparison['case'].replace('_', ' ')}**. " +
        ("All six designs cool the source to 70 °C. Their power-per-flow differences therefore arise from "
         "conversion of the same recovered heat per kilogram, including pump consumption. "
         "Changing inlet state alters the pressure level and heat-addition profile, while recuperation "
         "in Novec also redistributes heat internally.\n" if same_source_outlet else
         "The optimized source outlets differ, so the comparison includes both heat-recovery and "
         "heat-conversion effects.\n"),
        "The phase paths and heat-input profiles show how flexible inlet conditions modify each "
        "cycle's use of the source. Differences between butane and cyclopentane compare two simple "
        "architectures; differences involving Novec combine fluid choice and recuperation. They "
        "therefore do not isolate fluid properties alone. Native T–s plots place external streams "
        "on working-fluid entropy coordinates; expander lines connect endpoint states.\n",
        "## 3. Sensitivity rationale and initialization\n",
        "### 3.1 Physical question and grid\n",
        "Wet admission can reduce the isothermal evaporation duty and alter source matching, while "
        "an efficiency penalty reduces work recovered in expansion. The sensitivity tests where "
        "the former benefit outweighs the latter. There are **24 finite-heat temperatures × five "
        "designs = 120 optimized points per fluid**, or **360 total**. The 140 °C limit permits no "
        "geothermal temperature drop and makes ThermoOpt's available-heat denominator zero. Its "
        "15 requested combinations are recorded as unavailable, not as zero-efficiency cycles.\n",
        "### 3.2 Single-start grid refinement\n",
        "Each combination of fluid, inlet formulation, reinjection limit and expander efficiency "
        "starts from the highest-efficiency accepted saved solution with matching inputs and bounds "
        "in the previous v3/v2 studies. Its coarse-grid state is validated and copied before the "
        "evaporator and condenser are changed from 25 to 100 elements. Novec's recuperator retains "
        "25 elements. The initial design variables are transferred exactly, without perturbation.\n",
        "There is **one optimization per finite-heat case**, with no multistart, neighboring-case "
        "recovery or retry. A failed run remains flagged and is omitted from performance contours. "
        "The six-case comparison reuses matching sweep solutions. Local-optimum irregularities "
        "are reported without starting another optimization. All 360 finite-heat cases converged "
        "on their single attempt and passed the saved-state physical and sampled-constraint checks.\n",
        "### 3.3 Plot conventions\n",
        "Nominal exploitation is `100 × (140 − minimum outlet) / 120`; actual exploitation replaces "
        "the minimum with the optimized outlet. Both use 20 °C as the reference. Quality is shown "
        "in percent: 0% includes liquid and 100% includes superheated vapor. Dark blue indicates "
        "more liquid. Contour lines label the interpolated bands. The T–s/T–Q overlays use every "
        "other temperature at 90% expander efficiency, thin lines, and magma restricted to "
        "colormap coordinates 0.25–0.75. Their native cumulative heat-duty coordinate uses the "
        "1 MW sizing basis.\n"]
    for number,family in enumerate(FAMILIES,4):
        sections.append(fluid_section(family,number,tables[family],audits[family]))
    sections.extend(["## 7. Cross-fluid interpretation and limits\n",
        "At a 20 °C reinjection limit and 90% expander efficiency:\n",
        "| Configuration | Baseline power per source flow [kW per kg/s] | Partial power per source flow [kW per kg/s] | Gain [%] | Partial fixed-reference efficiency [%] |\n"
        "|---|---:|---:|---:|---:|\n" + "\n".join(
            f"| {LABELS[f]} | {row_at(tables[f],'superheated',20).specific_system_power_kW_per_kg_s:.2f} | "
            f"{row_at(tables[f],'partial_evaporation',20).specific_system_power_kW_per_kg_s:.2f} | "
            f"{row_at(tables[f],'partial_evaporation',20).relative_gain_pct:+.2f} | "
            f"{row_at(tables[f],'partial_evaporation',20).fixed_reference_efficiency_pct:.2f} |" for f in FAMILIES) + "\n",
        "Relative gain and absolute electricity production answer different questions. A larger "
        "gain can reflect a weaker baseline rather than a better final cycle. The fixed-reference "
        "efficiency columns permit direct comparison of power per geothermal flow across both "
        "fluids and reinjection limits. The phase maps are equally important: a wet-to-wet design "
        "and a wet-to-dry design place different requirements on an expander. Prescribed efficiencies "
        "explore worthwhile component performance, not what any particular machine can achieve.\n",
        "All accepted designs were reevaluated from their saved configurations to check physical "
        "flows, energy balances, state bounds, efficiency definitions and sampled constraints. "
        "They remain selected local solutions; contour interpolation does not establish globally "
        "optimal transition boundaries.\n",
        "Every case receives one optimization from its saved coarse-grid seed. These results "
        "therefore include both the effect of finer exchanger sampling and movement of the local "
        "solution. No additional starting points are used to improve an unfavorable result.\n",
        "### 7.1 Effect of evaporator and condenser refinement\n",
        "The table compares fixed-reference net efficiency before and after refinement at the "
        "70 and 20 °C source outlet limits, with 90% expander efficiency. The coarse value "
        "belongs to the actual saved seed selected for that case. A negative change means "
        "lower efficiency in the 100-element result. The two uniform grids are not nested, "
        "and movement between local solutions can also affect this comparison; the efficiency "
        "difference alone does not isolate discretization error.\n",
        "| Fluid | Formulation | Minimum source outlet [°C] | Coarse efficiency [%] | Fine efficiency [%] | Change [percentage points] |\n"
        "|---|---|---:|---:|---:|---:|\n" + "\n".join(
            f"| {LABELS[f].split(':')[0]} | {mode.replace('_',' ')} | {temperature} | "
            f"{row_at(tables[f],mode,temperature).coarse_fixed_reference_efficiency_pct:.4f} | "
            f"{row_at(tables[f],mode,temperature).fixed_reference_efficiency_pct:.4f} | "
            f"{row_at(tables[f],mode,temperature).fixed_reference_efficiency_pct-row_at(tables[f],mode,temperature).coarse_fixed_reference_efficiency_pct:+.4f} |"
            for f in FAMILIES for temperature in (70,20) for mode in ("superheated","partial_evaporation")) + "\n",
        "Relative gains can increase even when both absolute efficiencies decrease, if the "
        "superheated baseline is affected more strongly. All gain contours therefore compare "
        "the refined partial-evaporation result with its refined superheated baseline, using "
        "the same 100-element evaporator and condenser in both.\n",
        "### 7.2 Remaining sampling and optimization limits\n",
        "The **5 K approach is imposed at 100 evaporator/condenser locations and 25 recuperator locations, not continuously**. "
        "Nine selected designs per fluid and their coarse-grid seeds were reevaluated at 501 locations without reoptimization:\n",
        "| Fluid | Coarse seed: lowest checked heater approach [K] | Fine result: lowest checked heater approach [K] | Fine result: lowest checked condenser approach [K] | Fine result: lowest checked recuperator approach [K] |\n"
        "|---|---:|---:|---:|---:|\n" + "\n".join(
            f"| {LABELS[f].split(':')[0]} | {audits[f]['minimum_coarse_dense_approach_K']['heater']:.3f} | "
            f"{audits[f]['minimum_dense_approach_K']['heater']:.3f} | "
            f"{audits[f]['minimum_dense_approach_K']['cooler']:.3f} | " +
            (f"{audits[f]['minimum_dense_approach_K']['recuperator']:.3f}" if f == "novec" else "—") + " |"
            for f in FAMILIES) + "\n",
        "These subset diagnostics do not locate the worst approach over all designs. The reported "
        "results retain the requested 100-element evaporator and condenser model. Approaches below 5 K mean that a finer "
        "discretization and reoptimization are needed before treating precise break-even efficiencies "
        "or small gains as design conclusions. Close fluid rankings should be read with the same caution.\n",
        "Numerical data: [six-case comparison](results/fluid_comparison/comparison.csv), " + ", ".join(
            f"[{f} sensitivity](results/sensitivity_{f}/sweep_results.csv)" for f in FAMILIES) + ". "
        "Audit summaries and exchanger checks accompany each sensitivity. Figures also have PDF "
        "versions. See [README.md](README.md) for scripts, configuration and reproduction instructions.\n"])
    text = "\n".join(sections)
    headings = [line[3:] for line in text.splitlines() if line.startswith("## ")]
    import re
    contents = "**Contents**\n\n" + "\n".join(
        f"- [{heading}](#{re.sub(r'[^\w -]', '', heading.lower()).replace(' ', '-')})" for heading in headings)
    text = text.replace("## 1. Boundary", contents + "\n\n## 1. Boundary",1)
    (HERE / "report.md").write_text(text,encoding="utf-8")
    print("Wrote report.md from final saved data.")


if __name__ == "__main__":
    main()
