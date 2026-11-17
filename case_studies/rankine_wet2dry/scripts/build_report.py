"""Build matching Markdown and HTML reports from the saved native results."""
import html
import json
import os
from pathlib import Path
import study as run

pd, np = run.pd, run.np
NODES = []


def heading(text, level=2):
    NODES.append(dict(kind="heading", text=text, level=level))


def paragraph(text):
    NODES.append(dict(kind="paragraph", text=text))


def table(headers, rows):
    NODES.append(dict(kind="table", headers=headers, rows=rows))


def figure(path, caption):
    assert (run.HERE / (path+".png")).exists(), path
    assert (run.HERE / (path+".html")).exists(), path
    NODES.append(dict(kind="figure", path=path, caption=caption))


def at(t, mode, temperature, efficiency):
    rows = t.loc[t["mode"].eq(mode) & t["source_min_C"].eq(temperature) & np.isclose(t["turbine_efficiency"], efficiency)]
    return rows.iloc[0] if len(rows) == 1 and rows.iloc[0]["success"] else None


def value(row, key, digits=2):
    return "unavailable" if row is None else f"{row[key]:.{digits}f}"


def main():
    NODES.clear()
    tables = {f: pd.read_csv(run.RESULTS / f"sensitivity_{f}" / "sweep_results.csv") for f in run.FAMILIES}
    audits = {f: json.loads((run.RESULTS / f"sensitivity_{f}" / "validation_summary.json").read_text(encoding="utf-8")) for f in run.FAMILIES}
    comparison = pd.read_csv(run.RESULTS / "fluid_comparison/comparison.csv")
    heading("Partial evaporation in geothermal Rankine cycles", 1)
    paragraph("Updated EMPOWER study: n-butane, n-pentane, cyclopentane and Novec 649. The boundary conditions follow Nishith's 20 October 2026 request and Cristina's 19 October clarification. The standard-cycle reference and the 70 °C comparison use Cristina's 88.6% turbine isentropic efficiency. Partial-evaporation sensitivities retain 90%, 80%, 70% and 60% efficiency.")
    gains = []
    for f in run.FAMILIES:
        row = comparison.loc[comparison["case"].eq(f+"_partial_evaporation")].iloc[0]
        gains.append(f"{run.FAMILIES[f][0]} {row.relative_gain_pct:+.2f}%")
    best = comparison.loc[comparison.specific_system_power_kW_per_kg_s.idxmax()]
    paragraph("At equal 88.6% turbine efficiency and a 70 °C source-outlet limit, partial evaporation changes net power per geothermal flow by "+", ".join(gains)+f" relative to each fluid's own standard cycle. The highest output among these eight saved designs is {best.specific_system_power_kW_per_kg_s:.3f} kW per kg/s, for {best['case'].replace('_',' ')}. Deeper permitted source cooling can increase the advantage by improving heat recovery as well as conversion. These boundary-condition changes include the cooling-water inlet and turbine reference, so differences from the previous report cannot be attributed to pinch alone.")
    heading("1. Boundary conditions and interpretation")
    table(["Parameter", "Setting"], [
        ["Source water", "140 °C; 10 bar"],
        ["Cooling water", "15 °C inlet; exactly 25 °C return"],
        ["External heater / condenser minimum approach", "4.5 K; 100 locations each"],
        ["Novec internal recuperator", "5 K minimum approach; 25 locations; optimized effectiveness"],
        ["Source minimum outlet", "70 °C comparison; 140 to 20 °C in 5 K sensitivity steps"],
        ["Standard-cycle turbine efficiency", "88.6% throughout"],
        ["Partial-evaporation turbine efficiency", "90%, 80%, 70%, 60%; additional 88.6% cases at 70 °C"],
        ["Standard inlet / outlet", "2 K inlet superheat; dry expander outlet"],
        ["Partial evaporation", "Optimized inlet pressure and enthalpy; wet inlet and dry outlet not imposed"],
        ["Architectures", "Simple n-butane, n-pentane and cyclopentane; recuperated Novec 649"],
        ["Pumps / pressure drops", "70% isentropic pump efficiencies; 1% per active exchanger side"],
        ["Working-fluid pump inlet", "At least 1 K subcooling"],
        ["Sizing", "1 MW expander output minus working-fluid pump; external pumps deducted subsequently"],
        ["Optimization", "SLSQP, tolerance 1e-6, maximum 250 iterations per attempt"],
    ])
    paragraph("The 4.5 K instruction applies to the external hot exchangers and condensation side for both cycle types. The existing 5 K internal-recuperator assumption is retained because the exchange does not explicitly revise it. The original 2 K superheat, pump efficiencies, pressure losses and architecture assumptions are also retained. The spreadsheets do not fully specify these quantities, generator losses or auxiliary-power accounting, so this study aligns the stated boundary conditions rather than claiming an exact reproduction of EXERGY electrical outputs.")
    paragraph("A 70 °C source limit is a lower bound, not a fixed outlet. The optimizer can leave the source warmer. Permission to cool the source to 20 °C likewise does not mean the optimum reaches 20 °C.")
    heading("Optimization variables", 3)
    table(["Variable", "Role"], [
        ["Source outlet temperature", "Optimized above the case's minimum outlet"],
        ["Cooling-water outlet temperature", "Retained as a variable and constrained to exactly 25 °C"],
        ["Working-fluid pump inlet pressure and enthalpy", "Low-pressure state and at least 1 K subcooling"],
        ["Expander inlet pressure and enthalpy", "High-pressure state and admission phase"],
        ["Recuperator effectiveness", "Optimized between 0 and 1 for Novec only"],
    ])
    paragraph("Mass flows follow from the optimized states, heat balances and 1 MW cycle-sizing target. State bounds are broad and property-based. The standard-cycle superheat and dry-outlet constraints restrict the same inlet-state variables that remain flexible in the partial-evaporation formulation.")
    heading("2. Normalization and computation")
    paragraph(f"Fixed-reference efficiency is net system power divided by source mass flow and the water enthalpy drop from 140 to 15 °C at 10 bar. This reference heat is {run.fixed_reference_heat_J_per_kg()/1000:.3f} kJ/kg. It is common to all designs, and efficiency therefore tracks net power per unit geothermal flow directly. Nominal exploitation is 100 × (140 − minimum outlet) / 125; actual exploitation replaces the minimum by the achieved outlet. A 20 °C limit corresponds to 96% nominal exploitation on this reference.")
    paragraph("ThermoOpt's native system-efficiency objective uses available source heat down to each case's minimum outlet. At a fixed limit, maximizing that objective is equivalent to maximizing fixed-reference efficiency or net power per unit source flow. Reported powers at 194.44 kg/s scale each design to the source flow in Cristina's workbook. They remain model net powers before any separately specified generator losses.")
    paragraph("Every finite-heat case is optimized through the public thermopt.ThermodynamicCycleOptimization entry point using its native ThermodynamicCycleProblem, cycle equations and component models. No custom fitness function or alternate cycle evaluator is used. The previously converged designs supply initial guesses only. Each attempt preserves its complete case.yaml, initial and final configurations, convergence history, solver log and optimization report. Positive physical flows and work, constraint satisfaction and energy balance are checked before acceptance. Native source-code hashes are recorded so the implementation can be verified. A numerical solution is not proof of global optimality or real two-phase expander performance.")
    table(["Fluid", "Accepted finite-heat cases", "Failed finite-heat cases", "Zero-heat unavailable"],
          [[run.FAMILIES[f][0], a["accepted_points"], a["failed_points"], a["unavailable_points"]] for f, a in audits.items()])
    paragraph("There are 484 finite-heat designs: 24 temperatures × five designs × four fluids, plus four partial-evaporation designs at 88.6% for the 70 °C comparison. The 20 requested combinations at a 140 °C minimum outlet cannot release geothermal heat and are unavailable. Each fluid's sensitivity gain uses its standard cycle at 88.6% and the matching source-outlet limit.")
    heading("3. Eight-case comparison at 70 °C and 88.6% turbine efficiency")
    table(["Case", "Inlet vapor %", "Outlet vapor %", "Source out °C", "Net kW/(kg/s)", "Gain vs standard n-butane %", "Gain vs same-fluid standard %"],
          [[r["case"].replace("_", " "), f"{100*np.clip(r['inlet_quality'],0,1):.2f}", f"{100*np.clip(r['outlet_quality'],0,1):.2f}", f"{r['source_exit_C']:.2f}", f"{r['specific_system_power_kW_per_kg_s']:.3f}", f"{r['gain_vs_standard_butane_pct']:+.2f}", f"{r['relative_gain_pct']:+.2f}"] for r in comparison.to_dict("records")])
    paragraph("The common cross-fluid reference is standard n-butane, matching Cristina's spreadsheet. The same-fluid column isolates the change in inlet-state formulation at equal turbine efficiency. Superheated vapor is displayed as 100% quality; raw quality extensions remain in the data. Partial evaporation may produce wet-to-wet, wet-to-dry or dry-to-dry endpoints.")
    paragraph("The comparison exporter writes optimal_solution.xlsx with native cycle states and energy analysis, and design_variable_bounds.csv with the saved optimum, lower and upper bounds, fractional distance to the nearest bound and active-bound flags. Each case has its own directory under output/results/fluid_comparison/. The HTML comparison page offers the generated downloads when available. A bound is flagged when the optimum lies within 0.01% of its variable range.")
    for f in run.FAMILIES:
        figure(f"output/results/fluid_comparison/comparison_grid_{f}", f"{run.FAMILIES[f][2]}: standard and partial evaporation; T-s above T-Q; turbine efficiency 88.6%.")
    heading("4. Comparison with Cristina's revised workbook")
    paragraph("The original and revised workbooks report the same performance values. The revision adds condenser-inlet pressures and clarifies turbine-inlet pressure. The overlapping cases below are compared at 194.44 kg/s source flow. Pressures and inlet qualities in that workbook are optimized outputs, so they are validation targets rather than imposed inputs. The standard cyclopentane condenser-inlet pressure is blank in the workbook.")
    reference = {"butane_superheated": (6795.782384193805,13.55,2.77,1.), "pentane_superheated": (6396.,4.34,.78,1.), "cyclopentane_superheated": (5917.,2.81,None,1.), "pentane_partial_evaporation": (7061.,8.33,.81,.31), "cyclopentane_partial_evaporation": (6988.,6.09,.51,.19)}
    rows = []
    for r in comparison.to_dict("records"):
        if r["case"] not in reference:
            continue
        p, hp, lp, q = reference[r["case"]]
        rows.append([r["case"].replace("_", " "), f"{p:.1f}", f"{r['power_at_194p44_kg_s_kW']:.1f}", f"{hp:.2f} / {r['turbine_inlet_pressure_bar']:.2f}", f"{'unavailable' if lp is None else f'{lp:.2f}'} / {r['condenser_inlet_pressure_bar']:.2f}", f"{100*q:.1f} / {100*np.clip(r['inlet_quality'],0,1):.1f}"])
    table(["Case", "Workbook net kWe", "Model net kW", "Turbine p: workbook / model bar", "Condenser p: workbook / model bar", "Inlet vapor: workbook / model %"], rows)
    paragraph("Source: summary results_rev1.xlsx, Summary_results!C5:I17; assumptions clarified in RE: EMPOWER Project - ORC Conference, 19–20 October 2026. Net electrical and model mechanical powers may have different generator and auxiliary-loss definitions. The workbook's comparison row has a stray bar unit label; its values are dimensionless relative gains.")
    for number, (f, t) in enumerate(tables.items(), start=5):
        heading(f"{number}. {run.FAMILIES[f][2]}")
        b70 = at(t,"superheated",70,.886)
        b20 = at(t,"superheated",20,.886)
        p70 = at(t,"partial_evaporation",70,.9)
        p20 = at(t,"partial_evaporation",20,.9)
        equal_partial = at(t,"partial_evaporation",70,run.BASELINE_EFFICIENCY)
        paragraph(f"At a 70 °C source limit, standard-cycle fixed-reference efficiency is {value(b70,'fixed_reference_efficiency_pct')}%. Partial evaporation at 90% turbine efficiency reaches {value(p70,'fixed_reference_efficiency_pct')}%, a {value(p70,'relative_gain_pct')}% gain over the 88.6%-efficient standard reference. This sensitivity comparison includes the difference in turbine efficiency; the preceding eight-case comparison holds turbine efficiency equal.")
        paragraph(f"At equal 88.6% turbine efficiency and a 70 °C source limit, standard admission is {value(b70,'turbine_inlet_pressure_bar')} bar and {value(b70,'inlet_temperature_C')} °C. Partial-evaporation admission is {value(equal_partial,'turbine_inlet_pressure_bar')} bar and {value(equal_partial,'inlet_temperature_C')} °C. Flexible admission changes pressure level and the temperature of heat addition along with the vapor fraction.")
        paragraph(f"At the 20 °C limit, standard and 90%-efficient partial cycles produce {value(b20,'specific_system_power_kW_per_kg_s')} and {value(p20,'specific_system_power_kW_per_kg_s')} kW per kg/s of source flow. Their achieved source outlets are {value(b20,'source_exit_C')} and {value(p20,'source_exit_C')} °C. The partial-expansion inlet/outlet vapor fractions are {value(p20,'inlet_quality',4)} / {value(p20,'outlet_quality',4)} before display clipping.")
        recovery = p20.achieved_heat_utilization_pct/b20.achieved_heat_utilization_pct
        conversion = p20.cycle_efficiency_pct/b20.cycle_efficiency_pct
        power_ratio = p20.specific_system_power_kW_per_kg_s/b20.specific_system_power_kW_per_kg_s
        assert np.isclose(recovery*conversion, power_ratio, rtol=1e-8)
        paragraph(f"At the 20 °C limit, recovered source heat per kilogram changes by {100*(recovery-1):+.2f}% and net power per unit absorbed heat changes by {100*(conversion-1):+.2f}%. Their ratios multiply to give the {100*(power_ratio-1):+.2f}% net-power gain at equal geothermal flow. This comparison uses a 90%-efficient partial turbine and the 88.6%-efficient standard turbine.")
        if f == "novec":
            effectiveness = [run.read_formulation(
                run.RESULTS / "sensitivity_novec" / row["case"] / "optimal_solution.yaml"
            )["design_variables"]["recuperator_effectiveness"]["value"]
                for row in (b70, equal_partial, b20, p20)]
            table(["Minimum source outlet °C", "Standard effectiveness, 88.6% turbine", "Partial turbine efficiency %", "Partial effectiveness"], [
                [70, f"{effectiveness[0]:.3f}", 88.6, f"{effectiveness[1]:.3f}"],
                [20, f"{effectiveness[2]:.3f}", 90, f"{effectiveness[3]:.3f}"],
            ])
            paragraph("Internal recuperation reduces external heat demand and can limit how deeply the geothermal stream is cooled. Optimized effectiveness balances these effects, so its selected value changes with the permitted source outlet and inlet-state formulation.")
        table(["Minimum source outlet °C", "90% turbine: gain %", "80% turbine: gain %", "70% turbine: gain %", "60% turbine: gain %"],
              [[str(temp)]+[value(at(t,"partial_evaporation",temp,eff),"relative_gain_pct") for eff in run.TURBINE_EFFICIENCIES] for temp in (80,70,60,50,40,20)])
        for efficiency in (.8,.7,.6):
            subset = t[t["mode"].eq("partial_evaporation") & t["success"] & np.isclose(t["turbine_efficiency"],efficiency)]
            positive = subset[subset["relative_gain_pct"].gt(0)]
            if len(positive):
                paragraph(f"At {100*efficiency:.0f}% turbine efficiency, the highest sampled source-outlet limit with a positive gain is {positive.source_min_C.max():.0f} °C. The largest saved gain is {subset.relative_gain_pct.max():.2f}%.")
            else:
                paragraph(f"At {100*efficiency:.0f}% turbine efficiency, no accepted sampled case beats the 88.6%-efficient standard reference.")
        paragraph("Expansion endpoint counts across the 96 partial-evaporation sensitivity cases: "+", ".join(f"{key}: {n}" for key,n in audits[f]["expansion_classes"].items())+". These endpoint classes do not certify the phase path inside a real expander.")
        captions = {
            "baseline_performance":"Standard-cycle efficiency and achieved heat-source exploitation at 88.6% turbine efficiency.",
            "baseline_cycle_evolution":"Standard-cycle T-s and T-Q evolution at 88.6% turbine efficiency.",
            "partial_evaporation_cycle_evolution":"Partial-evaporation T-s and T-Q evolution at 90% turbine efficiency.",
            "system_efficiency_contours":"Fixed-reference net system efficiency for the partial-evaporation sensitivity grid.",
            "actual_exploitation_contours":"Achieved heat-source exploitation on the 140-to-15 °C reference.",
            "inlet_quality_contours":"Expander inlet vapor fraction. Superheated states are displayed as 100%.",
            "outlet_quality_contours":"Expander outlet vapor fraction. Wet admission does not imply dry discharge.",
            "relative_gain_contours":"Net-power gain against the same fluid's standard cycle at 88.6% turbine efficiency.",
        }
        for stem, caption in captions.items():
            figure(f"output/results/sensitivity_{f}/{stem}", caption)
    heading("9. Validation and limits")
    reproduction_path = run.RESULTS / "reproducibility_comparison.json"
    if reproduction_path.exists():
        reproduction = json.loads(reproduction_path.read_text(encoding="utf-8"))
        status = "passes" if reproduction["passed"] else "fails"
        power_checks = [item for item in reproduction["metrics"]
                        if item["metric"] == "specific_system_power_kW_per_kg_s"]
        relative_difference = max(item["maximum_relative_difference"] for item in power_checks)
        absolute_difference = max(item["maximum_absolute_difference"] for item in power_checks)
        paragraph(f"The consolidated run {status} comparison with the saved v2 reference across {reproduction['expected_finite_sweep_cases']} finite-heat sweep cases, {reproduction['expected_unavailable_sweep_cases']} unavailable endpoints and {reproduction['expected_comparison_cases']} comparison cases. The largest specific-power difference is {100*relative_difference:.3e}% relative and {absolute_difference:.3e} kW per kg/s absolute. Efficiencies, powers and source mass flow use 0.01% relative tolerance, pressures use 0.1%, temperatures use 0.02 K absolute tolerance and raw qualities use 0.001. Full metric comparisons and tolerances are recorded in output/results/reproducibility_comparison.json and its CSV companion.")
    kkt_fields = {"cases_with_recorded_kkt", "strict_kkt_passes"}
    if all(kkt_fields <= audit.keys() for audit in audits.values()):
        table(["Fluid", "Cases with native KKT records", "All strict KKT checks passed"],
              [[run.FAMILIES[f][0], a["cases_with_recorded_kkt"], a["strict_kkt_passes"]]
               for f, a in audits.items()])
        paragraph("The native first-order KKT diagnostics use a threshold of five times the solver tolerance. Case acceptance follows the native solver-success, physical-flow, energy-balance and constraint checks used in the v2 reference. The additional first-order diagnostics qualify the reported local solutions, and global optimality remains unproven.")
    table(["Fluid", "Native ThermoOpt optimizations", "Dense heater min K", "Dense condenser min K", "Dense recuperator min K"],
          [[run.FAMILIES[f][0],a['native_optimization_points'],f"{a['minimum_dense_approach_K']['heater']:.3f}", f"{a['minimum_dense_approach_K']['cooler']:.3f}", f"{a['minimum_dense_approach_K']['recuperator']:.3f}" if f == "novec" else "n.a."] for f,a in audits.items()])
    paragraph("The minimum approaches are enforced on the optimization grids. Selected accepted designs at 135, 70 and 20 °C are reevaluated at 501 locations without reoptimization. Approaches below the prescribed limit reveal sampling error between optimization locations. These are subset diagnostics, not a global worst-case bound. Small gains, precise break-even contours and close fluid rankings need grid refinement before design conclusions. Contour lines interpolate between four sampled turbine efficiencies; the plots are not additional optimizations.")
    paragraph("The study compares three simple-fluid architectures and one recuperated architecture. Differences involving Novec combine working-fluid choice and recuperation. Real wet-expander efficiencies, internal phase trajectories, equipment sizes and economics are not predicted by a prescribed constant isentropic efficiency.")
    heading("10. Embedding figures")
    paragraph("Every report figure is saved as PNG, PDF, SVG and a self-contained HTML file using the same Matplotlib drawing. The HTML uses inline SVG, requires no network or JavaScript, and scales with its container. Embed the HTML in an iframe or the SVG in an image/object element. These are static vector figures; they do not provide Plotly-style hover or data filtering. The standalone HTML report embeds its figures and data. output/figures.html is a local gallery.")
    md, body = [], []
    for n in NODES:
        kind = n["kind"]
        if kind == "heading":
            md.append("#"*n["level"]+" "+n["text"])
            body.append(f'<h{n["level"]}>{html.escape(n["text"])}</h{n["level"]}>')
        elif kind == "paragraph":
            md.append(n["text"])
            body.append("<p>"+html.escape(n["text"])+"</p>")
        elif kind == "table":
            md.append("| "+" | ".join(n["headers"])+" |\n|"+"---|"*len(n["headers"])+"\n"+"\n".join("| "+" | ".join(map(str,row))+" |" for row in n["rows"]))
            body.append('<div class="table"><table><thead><tr>'+"".join("<th>"+html.escape(str(h))+"</th>" for h in n["headers"])+"</tr></thead><tbody>"+"".join("<tr>"+"".join("<td>"+html.escape(str(v))+"</td>" for v in row)+"</tr>" for row in n["rows"])+"</tbody></table></div>")
        elif kind == "figure":
            path, caption = n["path"], n["caption"]
            md.append(f"![{caption}]({path}.png)\n\n[Embeddable HTML]({path}.html) · [SVG]({path}.svg) · [PDF]({path}.pdf)")
            body.append(f'<figure><img loading="lazy" src="{path}.svg" alt="{html.escape(caption)}"><figcaption>{html.escape(caption)} <a href="{path}.html">HTML</a> · <a href="{path}.svg">SVG</a> · <a href="{path}.pdf">PDF</a></figcaption></figure>')
    css = 'body{font:17px/1.6 system-ui,sans-serif;color:#243446;max-width:1120px;margin:45px auto;padding:0 24px}h1,h2{line-height:1.2;color:#153a55}h2{margin-top:3em}table{border-collapse:collapse;font-size:13px;width:100%}th{background:#163e58;color:white}td,th{padding:9px 12px;border-bottom:1px solid #d8e2e7;text-align:left}tr:nth-child(even){background:#f3f7fa}.table{overflow:auto}figure{margin:30px 0}img{width:100%;height:auto}figcaption{color:#526777;font-size:14px}a{color:#076994}@media print{figure,table{break-inside:avoid}body{margin:0}}'
    (run.HERE / "report_partial_evaporation.md").write_text("\n\n".join(md)+"\n", encoding="utf-8")
    (run.HERE / "report_partial_evaporation.html").write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Partial evaporation: updated EMPOWER study</title><style>'+css+'</style><body>'+"\n".join(body)+"</body></html>", encoding="utf-8")
    (run.RESULTS / "report_content.json").write_text(json.dumps(NODES, ensure_ascii=False, indent=2), encoding="utf-8")
    import build_html_report
    build_html_report.build()
    all_figures = sorted(run.RESULTS.rglob("*.svg"))
    import xml.etree.ElementTree as ET
    manifest = []
    for path in all_figures:
        root = ET.fromstring(path.read_text(encoding="utf-8"))
        viewbox = list(map(float, root.attrib["viewBox"].split()))
        manifest.append(dict(stem=path.with_suffix("").relative_to(run.HERE).as_posix(),
                             aspect_ratio=viewbox[2]/viewbox[3], formats=["png", "pdf", "svg", "html"]))
    (run.RESULTS / "figure_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    gallery = [f'<figure><img loading="lazy" src="{Path(os.path.relpath(p, run.OUTPUT)).as_posix()}" alt="{html.escape(p.stem)}"><figcaption><a href="{Path(os.path.relpath(p.with_suffix(".html"), run.OUTPUT)).as_posix()}">{html.escape(p.parent.name+" / "+p.stem)}</a></figcaption></figure>' for p in all_figures]
    (run.OUTPUT / "figures.html").write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>EMPOWER figure gallery</title><style>'+css+'</style><body><h1>EMPOWER figure gallery</h1>'+"".join(gallery)+"</body></html>", encoding="utf-8")
    print(f"Built reports and gallery: {len(NODES)} content blocks, {len(all_figures)} figures")


if __name__ == "__main__":
    main()
