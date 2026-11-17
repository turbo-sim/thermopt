"""Build an offline, multipage report from saved native ThermoOpt results.

This module formats saved results and draws scalar summary charts. It does not
implement or evaluate a thermodynamic cycle. Final builds reject alternate-model
results; --preview-archive is a clearly marked layout-only preview.
"""
from pathlib import Path
import html
import json
import os
import shutil
import sys
import study as run
from paths import WEB_ASSETS, HTML_PAGES

np, pd, plt = run.np, run.pd, run.plt
ESC = html.escape
LABELS = {"butane": "n-Butane", "pentane": "n-Pentane", "cyclopentane": "Cyclopentane", "novec": "Novec 649"}
SOURCE = run.RESULTS
OUTPUT = HTML_PAGES
PREVIEW = False
TABLES = {}
AUDITS = {}
COMPARISON = None
PAGES = []


def url(path):
    return Path(os.path.relpath(Path(path), OUTPUT)).as_posix()


def number(value, digits=2, signed=False):
    if value is None or not np.isfinite(value):
        return "n.a."
    return f"{value:+.{digits}f}" if signed else f"{value:.{digits}f}"


def selected(family, mode, minimum=70, efficiency=.886):
    t = TABLES[family]
    rows = t[t["mode"].eq(mode) & t["source_min_C"].eq(minimum) & np.isclose(t["turbine_efficiency"], efficiency)]
    return rows.iloc[0] if len(rows) == 1 and rows.iloc[0]["success"] else None


def val(row, key, digits=2, signed=False):
    return "n.a." if row is None else number(row[key], digits, signed)


def link(path, text, download=False):
    attr = " download" if download else ""
    return f'<a href="{ESC(url(path), quote=True)}"{attr}>{ESC(text)}</a>'


def section(title, subtitle=""):
    return f'<div class="section-head"><h2>{ESC(title)}</h2>'+ (f'<p>{ESC(subtitle)}</p>' if subtitle else "")+"</div>"


def note(text, warning=False):
    return f'<div class="callout{" warning" if warning else ""}">{text}</div>'


def metric(label, value, unit="", detail="", accent=False):
    return f'<div class="metric"><div class="metric-label">{ESC(label)}</div><div class="metric-value{" teal" if accent else ""}">{ESC(str(value))}<em>{ESC(unit)}</em></div><div class="metric-detail">{ESC(detail)}</div></div>'


def plot(stem, title, caption, featured=False):
    base = (OUTPUT / "assets" / stem) if PREVIEW and stem.startswith("report_summary/") else (SOURCE / stem)
    if not base.with_suffix(".svg").exists():
        raise FileNotFoundError(base.with_suffix(".svg"))
    expand = '<svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.4" aria-hidden="true"><path d="M6 2H2v4M10 2h4v4M2 10v4h4M10 14h4v-4M2 2l4 4m8-4-4 4M2 14l4-4m8 4-4-4"/></svg>'
    links = "".join(link(base.with_suffix("."+ext), ext.upper(), True) for ext in ("svg", "png", "pdf", "html"))
    return f'<figure class="plot{" featured" if featured else ""}"><div class="plot-top"><h3>{ESC(title)}</h3><button class="expand" data-zoom aria-label="Expand {ESC(title,quote=True)}">{expand} Expand</button></div><img src="{ESC(url(base.with_suffix(".svg")),quote=True)}" alt="{ESC(title,quote=True)}" loading="lazy"><figcaption>{ESC(caption)}<br><span class="plot-links">{links}</span></figcaption></figure>'


def data_table(headers, rows, classes="", table_id="", attributes=None):
    ident = f' id="{table_id}"' if table_id else ""
    body = []
    for i, row in enumerate(rows):
        attrs = ""
        if attributes:
            attrs = "".join(f' data-{k}="{ESC(str(v),quote=True)}"' for k,v in attributes[i].items())
        cells = []
        for j,v in enumerate(row):
            numeric = j > 0 and "methods" not in classes and 'class="file-links"' not in str(v)
            cells.append(f'<td class="num">{v}</td>' if numeric else f'<td>{v}</td>')
        body.append(f"<tr{attrs}>"+"".join(cells)+"</tr>")
    return f'<div class="table-wrap"><table class="data-table {classes}"{ident}><thead><tr>'+"".join(f"<th>{h}</th>" for h in headers)+"</tr></thead><tbody>"+"".join(body)+"</tbody></table></div>"


def tabset(items):
    buttons, panels = [], []
    for i, (key, label, content) in enumerate(items):
        buttons.append(f'<button id="tab-{key}" class="tab" role="tab" aria-controls="{key}" aria-selected="{str(i==0).lower()}" tabindex="{0 if i==0 else -1}">{ESC(label)}</button>')
        panels.append(f'<section class="tab-panel" id="{key}" role="tabpanel" aria-labelledby="tab-{key}"{(" hidden" if i else "")}>{content}</section>')
    return '<div class="tabs" role="tablist" data-tabs>'+"".join(buttons)+"</div>"+"".join(panels)


def navigation(current):
    groups = [("Study", [("index", "Overview", "01"), ("comparison", "Compare the designs", "02")]),
              ("Working fluids", [(f, LABELS[f], f"0{i+3}") for i,f in enumerate(run.FAMILIES)]),
              ("Evidence", [("methods", "Assumptions & method", "07"), ("validation", "Validation & data", "08")])]
    text = []
    for label, entries in groups:
        text.append(f'<div class="nav-label">{ESC(label)}</div>')
        for key, title, index in entries:
            active = ' aria-current="page"' if key == current else ""
            text.append(f'<a class="nav-link" href="{key}.html"{active}><span class="nav-index">{index}</span>{ESC(title)}</a>')
    return "".join(text)


def page(key, title, eyebrow, intro, content, conditions=True):
    count = sum(int(t["success"].sum()) for t in TABLES.values())
    preview = '<div class="preview-note">Layout preview using the previous dataset. The native ThermoOpt rerun is in progress.</div>' if PREVIEW else ""
    documents = SOURCE.parent if PREVIEW else run.HERE
    pdf = link(documents / "report_partial_evaporation.pdf", "Download PDF", True)
    condition = '<div class="condition-line"><span>140 °C source</span><span>15 → 25 °C cooling water</span><span>4.5 K external approach</span><span>88.6% standard reference</span></div>' if conditions else ""
    html_text = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>{ESC(title)} · EMPOWER</title><meta name="description" content="{ESC(intro,quote=True)}"><link rel="stylesheet" href="assets/report.css"><script defer src="assets/report.js"></script></head>
<body><a class="skip" href="#main">Skip to content</a><button class="mobile-scrim" aria-label="Close navigation"></button>
<aside class="sidebar" aria-label="Report navigation"><a class="brand" href="index.html"><i class="brand-mark" aria-hidden="true"></i><div><strong>EMPOWER</strong><span>Thermal cycle study</span></div></a><nav>{navigation(key)}</nav><div class="sidebar-footer"><strong>Partial evaporation</strong><br>Four fluids · {count} optimized designs<br>{'Layout preview' if PREVIEW else 'Native ThermoOpt'} · November 2026</div></aside>
<div class="app"><header class="topbar"><button class="menu-toggle" aria-expanded="false" aria-label="Open navigation">☰</button><div class="breadcrumb">Research report<span>/</span><b>{ESC(title)}</b></div><div class="top-actions"><a href="validation.html#case-register">Case data</a>{pdf}</div></header>{preview}
<main class="main" id="main"><div class="page-head"><div class="eyebrow">{ESC(eyebrow)}</div><h1>{ESC(title)}</h1><p>{ESC(intro)}</p>{condition}</div>{content}<footer class="footer"><span>EMPOWER · Partial evaporation in geothermal Rankine cycles</span><span>{'Preview dataset' if PREVIEW else 'ThermoOpt native cycle model'} · <a href="validation.html">Methods, records and limitations</a></span></footer></main></div>
<dialog id="figure-dialog" class="zoom-dialog" aria-labelledby="zoom-title"><div class="zoom-head"><p id="zoom-title">Figure</p><button class="zoom-close" aria-label="Close figure">×</button></div><img alt=""><div class="zoom-footer"><span>Matplotlib vector figure</span><a class="original-link" href="#" target="_blank" rel="noopener">Open full-size SVG ↗</a></div></dialog></body></html>'''
    target = OUTPUT / f"{key}.html"
    target.write_text(html_text, encoding="utf-8")
    PAGES.append(dict(path=target.name,title=title,purpose=intro))


def summary_figures():
    output = (OUTPUT / "assets/report_summary") if PREVIEW else (SOURCE / "report_summary")
    output.mkdir(parents=True,exist_ok=True)
    labels = list(LABELS.values())
    standard, partial, gain = [], [], []
    for f in run.FAMILIES:
        b,p = selected(f,"superheated"),selected(f,"partial_evaporation")
        standard.append(b.specific_system_power_kW_per_kg_s)
        partial.append(p.specific_system_power_kW_per_kg_s)
        gain.append(p.relative_gain_pct)
    style=dict(plt.rcParamsDefault)
    style.update({"backend":"Agg","font.family":"sans-serif","font.sans-serif":["Segoe UI","Arial"],"font.size":11,
                         "axes.labelsize":11,"xtick.labelsize":10,"ytick.labelsize":11,"axes.spines.top":False,
                         "axes.spines.right":False,"axes.spines.left":False,"axes.grid":False,"text.usetex":False,
                         "axes.edgecolor":"#cbd9df","axes.labelcolor":"#526c79","xtick.color":"#607580",
                         "ytick.color":"#264b59","savefig.dpi":250})
    with plt.rc_context(style):
        figure,ax=plt.subplots(figsize=(7.4,4.3),layout="constrained")
        x=np.arange(4)
        ax.barh(x+.17,standard,height=.3,color="#a9bdc5",label="Standard")
        ax.barh(x-.17,partial,height=.3,color="#126e78",label="Partial evaporation")
        for pos,values in ((x+.17,standard),(x-.17,partial)):
            for y,v in zip(pos,values): ax.text(v+.35,y,f"{v:.2f}",ha="left",va="center",fontsize=10,color="#355966")
        ax.set_yticks(x,labels);ax.invert_yaxis();ax.set_xlim(0,max(partial)*1.17)
        ax.set_xlabel("Net power per source flow [kW per kg/s]")
        ax.xaxis.grid(True,color="#e8eef0",linewidth=.7);ax.set_axisbelow(True);ax.tick_params(axis="y",length=0)
        ax.legend(loc="lower left",bbox_to_anchor=(0,1.015),ncol=2,frameon=False,fontsize=10)
        run.save_figure(figure,output/"power_by_fluid");plt.close(figure)
        figure,ax=plt.subplots(figsize=(7.4,4.3),layout="constrained")
        ax.barh(x,gain,height=.48,color=["#5c939d","#3d8491","#126e78","#2f7889"])
        for y,v in zip(x,gain): ax.text(v+.16,y,f"+{v:.2f}%",va="center",fontsize=11,color="#264c5b")
        ax.set_yticks(x,labels);ax.invert_yaxis();ax.set_xlim(0,max(gain)*1.25)
        ax.set_xlabel("Gain over the same fluid’s standard cycle [%]")
        ax.xaxis.grid(True,color="#e8eef0",linewidth=.7);ax.set_axisbelow(True);ax.tick_params(axis="y",length=0)
        run.save_figure(figure,output/"gain_by_fluid");plt.close(figure)


def overview():
    best=COMPARISON.loc[COMPARISON.specific_system_power_kW_per_kg_s.idxmax()]
    largest=COMPARISON.loc[COMPARISON.relative_gain_pct.idxmax()]
    stem1="report_summary/power_by_fluid";stem2="report_summary/gain_by_fluid"
    content='<div class="intro-grid"><p class="intro-statement">At equal turbine efficiency, partial evaporation raises net output for all four fluids at the 70 °C reinjection limit.</p><ul class="fact-list"><li><strong>Compare fairly</strong><span>Both modes at 88.6%</span></li><li><strong>Test the efficiency penalty</strong><span>90 / 80 / 70 / 60%</span></li><li><strong>Explore heat recovery</strong><span>135 → 20 °C outlet limits</span></li></ul></div>'
    content+=section("The eight-case comparison","70 °C minimum source outlet · net power after all modelled pumps")
    content+='<div class="grid cols-2">'+plot(stem1,"Electricity potential per geothermal flow","Standard and partial cycles at the same 88.6% turbine efficiency. These are model net powers, before separately specified generator losses.")+plot(stem2,"What partial evaporation adds","Each partial cycle is compared with the same fluid’s standard cycle; this differs from the common n-butane reference.")+"</div>"
    bestfamily=best["family"]
    content+='<div class="text-columns"><div><h3>Highest absolute output</h3><p>'+f'{LABELS[bestfamily]} reaches {best.specific_system_power_kW_per_kg_s:.2f} kW per kg/s in the partial cycle. Its recuperated architecture contributes to this comparison.'+'</p></div><div><h3>Largest relative improvement</h3><p>'+f'{LABELS[largest["family"]]} gains {largest.relative_gain_pct:.2f}% against its own standard reference. Relative gain and absolute output answer different questions.'+'</p></div></div>'
    content+=section("Where does an expander penalty remain worthwhile?","Highest sampled minimum source outlet with positive gain versus an 88.6% standard cycle")
    rows=[]
    for f,t in TABLES.items():
        cells=[f'<a class="fluid-name" href="{f}.html">{LABELS[f]}</a>']
        for efficiency in (.8,.7,.6):
            positive=t[t.success & t["mode"].eq("partial_evaporation") & np.isclose(t.turbine_efficiency,efficiency) & t.relative_gain_pct.gt(0)]
            cells.append(f"{positive.source_min_C.max():.0f} °C" if len(positive) else "No sampled gain")
        rows.append(cells)
    content+=data_table(["Fluid","80% turbine","70% turbine","60% turbine"],rows)
    content+='<p class="table-note">The limits are sampled in 5 K steps. They are model break-even observations, not an expander performance guarantee.</p>'
    content+=section("Explore one fluid at a time","Each fluid page separates performance, heat recovery, phase state and cycle evolution.")
    content+='<div class="grid cols-2">'
    for f in run.FAMILIES:
        p=selected(f,"partial_evaporation")
        content+=f'<a class="link-card" href="{f}.html"><span class="arrow">↗</span><span class="label">{"Recuperated" if f=="novec" else "Simple"} cycle</span><h3>{LABELS[f]}</h3><p>+{p.relative_gain_pct:.2f}% at equal turbine efficiency · investigate the full sensitivity.</p></a>'
    content+='</div>'
    content+=note('The external approach is constrained at 100 heater/condenser locations. Denser checks find lower approaches between some locations. <a href="validation.html">Read the pinch-grid limits before using close rankings or break-even values.</a>')
    page("index","Where partial evaporation increases net power","Study overview","A consistent comparison of four working fluids, deeper geothermal cooling and the efficiency penalty of wet expansion.",content)


def comparison_page():
    content=section("Equal efficiency. Equal boundary conditions.","Eight optimized designs · 88.6% turbine efficiency · 70 °C minimum source outlet")
    content+='<div class="grid cols-2">'+plot("report_summary/power_by_fluid","Absolute net power","Power per unit geothermal-source flow makes the different 1 MW-sized cycles comparable.")+plot("report_summary/gain_by_fluid","Same-fluid gain","This isolates the inlet-state formulation at equal turbine efficiency.")+"</div>"
    content+=section("Compare the results","The common cross-fluid reference is standard n-butane, matching Cristina’s workbook.")
    rows=[]
    for r in COMPARISON.to_dict("records"):
        fluid=LABELS[r["family"]]
        mode="Standard" if r["mode"]=="superheated" else "Partial evaporation"
        gain=r["relative_gain_pct"]
        rows.append([f'<span class="fluid-name">{fluid}</span><small>{mode}</small>',f'{r["specific_system_power_kW_per_kg_s"]:.3f}',
                     f'{r["source_exit_C"]:.2f}',f'{100*np.clip(r["inlet_quality"],0,1):.1f}',f'{100*np.clip(r["outlet_quality"],0,1):.1f}',
                     f'{r["gain_vs_standard_butane_pct"]:+.2f}%',f'<span class="positive">{gain:+.2f}%</span>'])
    content+=data_table(["Fluid / mode","Net kW<br>per kg/s","Source out<br>°C","Inlet<br>vapor %","Outlet<br>vapor %","Gain vs<br>standard n-butane","Gain vs<br>same-fluid standard"],rows)
    content+='<p class="table-note">100% displayed vapor includes superheated gas. The raw quality values are preserved in the CSV. '+link(SOURCE/"fluid_comparison/comparison.csv","Download comparison CSV",True)+"</p>"
    content+=section("Inspect the heat-addition profile","Choose a fluid; standard and partial cycles share one diagram grid.")
    items=[]
    for f in run.FAMILIES:
        caption="Standard and partial-evaporation columns; T-s above T-Q. Both turbine efficiencies are 88.6%."
        items.append(("compare-"+f,LABELS[f],plot(f"fluid_comparison/comparison_grid_{f}",LABELS[f]+": cycle comparison",caption,True)))
    content+=tabset(items)
    content+=section("Alignment with Cristina’s results","Overlapping cases at the workbook’s 194.44 kg/s source flow")
    reference={"butane_superheated":(6795.782384193805,13.55,2.77),"pentane_superheated":(6396.,4.34,.78),"cyclopentane_superheated":(5917.,2.81,None),"pentane_partial_evaporation":(7061.,8.33,.81),"cyclopentane_partial_evaporation":(6988.,6.09,.51)}
    rows=[]
    for r in COMPARISON.to_dict("records"):
        if r["case"] not in reference:continue
        power,hp,lp=reference[r["case"]]
        rows.append([f'{LABELS[r["family"]]}<small>{"Standard" if r["mode"]=="superheated" else "Partial evaporation"}</small>',
                     f"{power:,.1f}",f'{r["power_at_194p44_kg_s_kW"]:,.1f}',f'{hp:.2f} / {r["turbine_inlet_pressure_bar"]:.2f}',
                     f'{"n.a." if lp is None else number(lp)} / {r["condenser_inlet_pressure_bar"]:.2f}'])
    content+=data_table(["Case","Workbook<br>net kWe","Model<br>net kW","Turbine inlet bar<br>workbook / model","Condenser inlet bar<br>workbook / model"],rows)
    content+=note('The workbook reports electrical power. Generator losses, auxiliary definitions, superheat and recuperation assumptions are not fully specified. Pressure and quality results are validation targets rather than imposed inputs. <a href="methods.html">See the confirmed and retained assumptions.</a>')
    if not PREVIEW:
        rows=[]
        for r in COMPARISON.to_dict("records"):
            directory=SOURCE/"fluid_comparison"/r["case"]
            label=LABELS[r["family"]]+" · "+("Standard" if r["mode"]=="superheated" else "Partial evaporation")
            workbook=directory/"optimal_solution.xlsx"
            bounds=directory/"design_variable_bounds.csv"
            if workbook.exists() or bounds.exists():
                rows.append([ESC(label),link(workbook,"Cycle workbook",True) if workbook.exists() else "Unavailable",
                             link(bounds,"Design-variable bounds",True) if bounds.exists() else "Unavailable"])
        if rows:
            content+=section("Download the comparison cases","Native cycle states, energy analysis and design-variable bounds for each saved optimum")
            content+=data_table(["Case","Excel workbook","Bounds CSV"],rows)
            content+='<p class="table-note">Bounds are reported in SI units. The fractional distance uses each variable’s full range. A lower or upper bound is flagged within 0.01% of that range.</p>'
    page("comparison","Compare the optimized designs","Equal-efficiency comparison","Separate absolute power, same-fluid improvement and agreement with the supplied workbook.",content)


def fluid_page(f):
    t=TABLES[f];b=selected(f,"superheated");p=selected(f,"partial_evaporation")
    content='<div class="metrics">'+metric("Standard net output",val(b,"specific_system_power_kW_per_kg_s"),"kW/(kg/s)","88.6% turbine · 70 °C limit")+metric("Partial net output",val(p,"specific_system_power_kW_per_kg_s"),"kW/(kg/s)","88.6% turbine · 70 °C limit")+metric("Same-fluid gain",val(p,"relative_gain_pct",2,True),"%","Equal turbine efficiency",True)+"</div>"
    stem=f"sensitivity_{f}/"
    performance='<div class="grid cols-2">'+plot(stem+"system_efficiency_contours","Net efficiency on a fixed heat reference","Net system power divided by source heat from 140 to 15 °C at 10 bar. The partial turbine efficiency varies from 60% to 90%.")+plot(stem+"relative_gain_contours","When does partial evaporation pay off?","Gain against this fluid’s 88.6%-efficient standard cycle at the matching minimum source outlet.")+"</div>"
    performance+=section("Read the sampled gains","Relative net-power gain at equal geothermal-source flow; standard reference efficiency 88.6%")
    rows=[]
    for temp in (80,70,60,50,40,20):
        rows.append([f"{temp} °C"]+[val(selected(f,"partial_evaporation",temp,e),"relative_gain_pct",2,True)+"%" for e in run.TURBINE_EFFICIENCIES])
    performance+=data_table(["Minimum source outlet","90% turbine","80% turbine","70% turbine","60% turbine"],rows)
    performance+='<p class="table-note">The sensitivity includes differences in turbine efficiency. The top summary and comparison page hold both modes at 88.6%.</p>'
    recovery=plot(stem+"baseline_performance","What the standard cycle can recover","Standard turbine efficiency 88.6%. Allowed source cooling and the outlet selected by the optimizer are different quantities.")
    b20,p20=selected(f,"superheated",20,.886),selected(f,"partial_evaporation",20,.9)
    recovery+='<div class="grid cols-2" style="margin-top:20px">'+plot(stem+"actual_exploitation_contours","Actual heat-source exploitation","Actual and nominal exploitation use 15 °C as their cold reference; the 20 °C limit represents 96% nominal exploitation.")+'<div class="panel"><h2>Permission to cool is not achieved cooling</h2><p>At the 20 °C outlet limit, the optimized source leaves at:</p><ul class="fact-list"><li><strong>Standard, 88.6% turbine</strong><span>'+val(b20,"source_exit_C")+' °C</span></li><li><strong>Partial, 90% turbine</strong><span>'+val(p20,"source_exit_C")+' °C</span></li></ul><p style="margin-top:18px">Partial evaporation can shorten the isothermal evaporation duty and improve matching to the cooling geothermal stream. Additional recovered heat and its conversion efficiency both affect net power.</p><p>All cycles are sized to 1 MW after the working-fluid pump. Their geothermal flows differ; power per source flow is the comparison basis.</p></div></div>'
    heat_ratio=p20.achieved_heat_utilization_pct/b20.achieved_heat_utilization_pct
    conversion_ratio=p20.cycle_efficiency_pct/b20.cycle_efficiency_pct
    power_ratio=p20.specific_system_power_kW_per_kg_s/b20.specific_system_power_kW_per_kg_s
    assert np.isclose(heat_ratio*conversion_ratio,power_ratio,rtol=1e-8)
    recovery+=section("Where the power gain comes from","20 °C source-outlet limit · standard turbine 88.6% · partial turbine 90%")
    recovery+=data_table(["Contribution","Change from standard"],[
        ["Recovered source heat per kilogram",number(100*(heat_ratio-1),2,True)+"%"],
        ["Net power per unit absorbed heat",number(100*(conversion_ratio-1),2,True)+"%"],
        ["Net power per source flow",number(100*(power_ratio-1),2,True)+"%"],
    ])
    recovery+='<p>The recovered-heat ratio multiplied by the heat-conversion ratio gives the net-power ratio. Both terms include the optimized state changes, and conversion includes the prescribed turbine efficiency and external pumping.</p>'
    if f=="novec":
        effectiveness=[run.read_formulation(
            SOURCE/"sensitivity_novec"/row["case"]/"optimal_solution.yaml"
        )["design_variables"]["recuperator_effectiveness"]["value"] for row in (b,p,b20,p20)]
        recovery+=section("How much recuperation is selected","Effectiveness is optimized between zero and one")
        recovery+=data_table(["Source-outlet limit","Standard effectiveness<br>88.6% turbine","Partial turbine efficiency","Partial effectiveness"],[
            ["70 °C",number(effectiveness[0],3),"88.6%",number(effectiveness[1],3)],
            ["20 °C",number(effectiveness[2],3),"90%",number(effectiveness[3],3)],
        ])
        recovery+='<p>Internal preheating reduces external heat demand and can limit deep cooling of the geothermal stream. The optimizer balances these effects as the source-outlet limit and inlet-state formulation change.</p>'
    phases='<div class="grid cols-2">'+plot(stem+"inlet_quality_contours","Expander admission","The inlet state is optimized. Allowing partial evaporation does not require a wet inlet.")+plot(stem+"outlet_quality_contours","Expander discharge","Displayed 100% includes superheated vapor. Wet admission can end either wet or dry.")+"</div>"
    counts=AUDITS[f]["expansion_classes"]
    phases+=section("Endpoint phases across the sensitivity","96 partial-evaporation cases, excluding the extra 88.6% comparison point")
    phases+=data_table(["Expansion endpoints","Number of designs"],[[ESC(k.capitalize()),str(v)] for k,v in counts.items()])
    phases+=note("Endpoint phase classes do not establish a real expander’s internal phase path or attainable efficiency. The efficiency is prescribed in this study.")
    cycles='<div class="grid cols-2">'+plot(stem+"baseline_cycle_evolution","Standard cycle · 88.6% turbine","T-s and T-Q overlays at alternating source-outlet limits. Expand for a full-size view.")+plot(stem+"partial_evaporation_cycle_evolution","Partial evaporation · 90% turbine","The same heat-duty coordinate and source-limit colour convention make the evolving heat-input profiles comparable.")+"</div>"
    cycles+=note("The overlays use different stated turbine efficiencies. Heat duties retain each design’s 1 MW cycle sizing. Source and cooling-water T-s traces use the associated working-fluid entropy coordinate, as in ThermoOpt’s native diagrams.")
    cycles+=section("Admission pressure and temperature","Equal 88.6% turbine efficiency · 70 °C source-outlet limit")
    cycles+=data_table(["Mode","Admission pressure bar","Admission temperature °C"],[
        ["Standard",val(b,"turbine_inlet_pressure_bar"),val(b,"inlet_temperature_C")],
        ["Partial evaporation",val(p,"turbine_inlet_pressure_bar"),val(p,"inlet_temperature_C")],
    ])
    cycles+='<p>Flexible admission changes the high-pressure level and temperature of heat addition along with the vapor fraction, which changes how the cycle matches the cooling geothermal stream.</p>'
    content+=tabset([("performance","Performance",performance),("recovery","Heat recovery",recovery),("phases","Phase behaviour",phases),("cycles","Cycle evolution",cycles)])
    arch="recuperated" if f=="novec" else "simple"
    page(f,LABELS[f],f"Working fluid · {arch} cycle","Follow the gain in net power, the heat-source recovery and the expander phase states across the sensitivity.",content)


def methods_page():
    content=section("Confirmed inputs and retained assumptions","Nishith’s request determines the revised boundary conditions.")
    settings=[
        ("External heater and condenser approach","4.5 K for both cycle types","Cristina, 19 October; adopted in Nishith’s 20 October request"),
        ("Cooling water","15 °C inlet; exactly 25 °C return","Cristina’s exchange and revised workbook"),
        ("Standard turbine reference","88.6% isentropic efficiency","Cristina’s workbook; user confirmed throughout"),
        ("Partial turbine sensitivity","90%, 80%, 70%, 60%","Retained sensitivity requested by the user"),
        ("Source","Water at 140 °C and 10 bar","Workbook inlet temperature; retained pressure and water model"),
        ("Minimum source outlet","70 °C comparison; 135 to 20 °C sensitivity","Retained study grid; constraint is a lower bound"),
        ("Internal recuperator","Novec only; 5 K approach; optimized effectiveness","Retained assumption; not explicitly revised in the exchange"),
        ("Standard admission and discharge","2 K inlet superheat; dry outlet","Retained study assumption"),
        ("Partial admission and discharge","Optimized inlet; no dry-outlet requirement","Flexible inlet formulation"),
        ("Pumps and pressure losses","70% pump efficiencies; 1% per active exchanger side","Retained study assumptions"),
        ("Exchanger grid","100 heater/condenser locations; 25 recuperator locations","Retained optimization grid; dense checks at 501 locations"),
    ]
    content+=data_table(["Parameter","Applied setting","Basis"],[[ESC(x) for x in row] for row in settings],"methods")
    content+=section("Optimization variables","Mass flows follow from the optimized states, heat balances and 1 MW cycle-sizing target.")
    variables=[
        ("Source outlet temperature","Optimized above the case’s minimum outlet"),
        ("Cooling-water outlet temperature","Retained as a variable and constrained to exactly 25 °C"),
        ("Working-fluid pump inlet pressure and enthalpy","Low-pressure state and at least 1 K subcooling"),
        ("Expander inlet pressure and enthalpy","High-pressure state and admission phase"),
        ("Recuperator effectiveness","Optimized between zero and one for Novec only"),
    ]
    content+=data_table(["Variable","Role"],[[ESC(x) for x in row] for row in variables],"methods")
    content+='<p>Pressure and enthalpy bounds are broad and property-based. The standard-cycle superheat and dry-outlet constraints restrict the same inlet-state variables that remain flexible in the partial-evaporation formulation. Bounds at the eight comparison optima are available on the comparison page.</p>'
    content+=section("A comparison at equal geothermal flow")
    content+='<div class="grid cols-2"><div class="panel"><h2>Fixed-reference efficiency</h2><p>Net system power divided by source flow and the water enthalpy drop from 140 to 15 °C at 10 bar.</p><div class="summary-value">'+number(run.fixed_reference_heat_J_per_kg()/1000,3)+' <span style="font:14px Segoe UI">kJ/kg</span></div><p>The heat reference is constant, so this efficiency is directly proportional to net power per unit geothermal flow.</p></div><div class="panel"><h2>Source exploitation</h2><p>Nominal exploitation uses the permitted outlet; actual exploitation uses the optimized outlet.</p><div class="code">100 × (140 − outlet temperature) / 125</div><p>The 20 °C endpoint therefore represents 96% nominal exploitation. Heat extraction and conversion efficiency are separate contributors to net power.</p></div></div>'
    content+=note("ThermoOpt’s native objective uses heat available down to each case’s minimum source outlet. At a fixed outlet limit, it selects the same optimum as fixed-reference efficiency or power per geothermal flow.")
    content+=section("Native ThermoOpt optimization","The public entry point uses the repository’s cycle and component models.")
    content+='<div class="code">import thermopt as th<br><br>cycle = th.ThermodynamicCycleOptimization(<br>&nbsp;&nbsp;&nbsp;&nbsp;"case.yaml", out_dir="native_run"<br>)<br>cycle.run_optimization()</div>'
    content+='<div class="text-columns"><div><h3>Initialization</h3><p>Versioned states under inputs/initializers supply the original study’s warm initial guesses. Every finite-heat case is optimized with native ThermoOpt. The saved v2 reference tables provide the numerical comparison for a fresh run.</p></div><div><h3>Acceptance</h3><p>The native solver must succeed. The result must satisfy the constraints, produce positive physical flows and net power, and close the cycle energy balance. Complete input and solver records are retained.</p></div></div>'
    content+=section("Source documents")
    files=[("RE_ EMPOWER Project - ORC Conference.eml","Email exchange · clarified external pinches and cooling temperatures"),("summary results_rev1.xlsx","Revised workbook · efficiency, pressure and power targets"),("summary results.xlsx","Original workbook · same reported performance values")]
    content+='<div class="panel"><ul class="download-list">'+"".join(f'<li>{link(run.REFERENCE_INPUTS/filename,filename,True)}<span>{ESC(desc)}</span></li>' for filename,desc in files)+"</ul></div>"
    content+=note("The workbook reports net electrical power. Generator losses and some auxiliary and equipment assumptions are not specified. This report aligns the stated boundary conditions and keeps those remaining differences visible.")
    page("methods","Assumptions & method","Basis of the study","Trace each applied setting to the exchange or the retained cycle model, then inspect the native optimization method.",content)


def validation_page():
    accepted=sum(int(t.success.sum()) for t in TABLES.values())
    dense=sum(a["dense_points"] for a in AUDITS.values())
    content='<div class="metrics">'+metric("Native optimized designs",accepted,"", "24 finite-heat limits × five designs × four fluids + four comparison cases",True)+metric("Dense-grid checks",dense,"designs","501 exchanger locations; diagnostic reevaluation")+metric("Zero-heat endpoints",20,"excluded","140 °C minimum outlet releases no geothermal heat")+"</div>"
    content+='<div class="provenance-strip">'+("Layout preview. The final register is populated only by native optimization records." if PREVIEW else "Every listed design was optimized by ThermoOpt’s native cycle problem.")+"</div>"
    content+=section("An auditable computation path")
    content+='<div class="grid cols-2"><div class="panel"><h2>Native model, throughout the solve</h2><p>Each case uses <code>thermopt.ThermodynamicCycleOptimization</code> and its unmodified <code>ThermodynamicCycleProblem</code>. There is no custom fitness override or alternate evaluator in the active run.</p><p>The convergence header records the native problem class. The source-code hashes, complete YAML input and solver records establish which model produced each result.</p></div><div class="panel"><h2>Reproduce one case directly</h2><p>Download its YAML input from the register below and run the public ThermoOpt entry point.</p><div class="code">import thermopt as th<br>cycle = th.ThermodynamicCycleOptimization(<br>&nbsp;&nbsp;"case.yaml", out_dir="verification_run")<br>cycle.run_optimization()</div></div></div>'
    reproduction=SOURCE/"reproducibility_comparison.json"
    if not PREVIEW and reproduction.exists():
        verification=json.loads(reproduction.read_text(encoding="utf-8"))
        status="passed" if verification["passed"] else "failed"
        power_checks=[item for item in verification["metrics"]
                      if item["metric"]=="specific_system_power_kW_per_kg_s"]
        relative_difference=max(item["maximum_relative_difference"] for item in power_checks)
        absolute_difference=max(item["maximum_absolute_difference"] for item in power_checks)
        content+=section("Comparison with the saved v2 reference","Recorded status · "+status)
        content+=f'<p>Comparison {status} for {verification["expected_finite_sweep_cases"]} finite-heat sweep cases, {verification["expected_unavailable_sweep_cases"]} unavailable endpoints and {verification["expected_comparison_cases"]} comparison cases. The largest specific-power difference is {100*relative_difference:.3e}% relative and {absolute_difference:.3e} kW per kg/s absolute.</p>'
        content+='<p>Efficiencies, powers and source mass flow use 0.01% relative tolerance. Pressures use 0.1% relative tolerance, temperatures use 0.02 K absolute tolerance and raw qualities use 0.001 absolute tolerance. '+link(reproduction,"Download comparison details and tolerances",True)+" · "+link(reproduction.with_suffix(".csv"),"Download metric summary CSV",True)+"</p>"
    kkt_fields={"cases_with_recorded_kkt","strict_kkt_passes"}
    if not PREVIEW and all(kkt_fields<=audit.keys() for audit in AUDITS.values()):
        content+=section("Native first-order diagnostics","KKT residual checks use five times the solver tolerance.")
        content+=data_table(["Fluid","Cases with native<br>KKT records","All strict KKT<br>checks passed"],[
            [LABELS[f],a["cases_with_recorded_kkt"],a["strict_kkt_passes"]] for f,a in AUDITS.items()
        ])
        content+='<p>Case acceptance follows the native solver-success, physical-flow, energy-balance and constraint checks used in the v2 reference. These additional first-order diagnostics qualify the local solutions. Global optimality remains unproven.</p>'
    content+=section("Pinch discretization remains a limitation","These are the lowest approaches among the selected dense reevaluations, not all designs.")
    rows=[]
    for f,a in AUDITS.items():
        minima=a["minimum_dense_approach_K"]
        rows.append([LABELS[f],number(minima["heater"],3),number(minima["cooler"],3),number(minima.get("recuperator"),3)])
    content+=data_table(["Fluid","Heater min K<br>target 4.5 K","Condenser min K<br>target 4.5 K","Recuperator min K<br>target 5 K"],rows)
    content+=note("Some approaches fall below the imposed grid constraint between the optimization locations. The 100-point results are retained. Precise break-even efficiencies, small gains and close fluid rankings need refinement before design conclusions.",True)
    content+='<div id="case-register">'+section("Native optimization register","Filter the cases and inspect the exact inputs, converged states and solver records.")+"</div>"
    options=''.join(f'<option value="{f}">{label}</option>' for f,label in LABELS.items())
    temperatures=''.join(f'<option value="{t}"{" selected" if t==70 else ""}>{t} °C</option>' for t in range(135,19,-5))
    content+='<div class="filters"><label>Fluid<select data-filter="family"><option value="">All fluids</option>'+options+'</select></label><label>Cycle<select data-filter="mode"><option value="">Both modes</option><option value="superheated">Standard</option><option value="partial_evaporation">Partial evaporation</option></select></label><label>Minimum source outlet<select data-filter="minimum"><option value="">All temperatures</option>'+temperatures+'</select></label><label>Turbine efficiency<select data-filter="efficiency"><option value="">All efficiencies</option>'+''.join(f'<option value="{e:g}">{100*e:g}%</option>' for e in (.886,.9,.8,.7,.6))+'</select></label><label>Search<input data-filter="search" type="search" placeholder="Fluid or case name"></label></div><p class="filter-count" id="case-count" aria-live="polite"></p>'
    rows, attrs=[],[]
    for f,t in TABLES.items():
        for row in t[t.success].sort_values(["source_min_C","mode","turbine_efficiency"],ascending=[False,True,False]).to_dict("records"):
            directory=SOURCE/f"sensitivity_{f}"/row["case"]
            saved=json.loads((directory/"result.json").read_text(encoding="utf-8"))
            if not PREVIEW:
                assert saved["row"]["optimization_backend"]=="thermopt.native"
                assert saved["row"]["problem_class"]=="thermopt.optimize_cycle.ThermodynamicCycleProblem"
                case=directory/saved["backend"]["input_file"]
            else:
                case=directory/"attempt_1/initial_guess.yaml"
            recordlinks='<div class="file-links">'+link(case,"Input YAML",True)+link(directory/"optimal_solution.yaml","Optimum",True)+link(directory/"convergence_history.txt","Solver",True)+"</div>"
            rows.append([f'<span class="fluid-name">{LABELS[f]}</span><small>{"Standard" if row["mode"]=="superheated" else "Partial evaporation"}</small>',f'{row["source_min_C"]:g} °C',f'{100*row["turbine_efficiency"]:g}%',f'{row["specific_system_power_kW_per_kg_s"]:.3f}',f'{row["max_constraint_violation"]:.1e}',recordlinks])
            attrs.append(dict(family=f,mode=row["mode"],minimum=f'{row["source_min_C"]:g}',efficiency=f'{row["turbine_efficiency"]:g}',search=f+" "+row["case"]))
    content+=data_table(["Fluid / mode","Minimum<br>source outlet","Turbine<br>efficiency","Net kW<br>per kg/s","Constraint<br>violation","Reproduction records"],rows,table_id="case-table",attributes=attrs)
    content+='<p id="no-cases" class="no-results" hidden>No cases match these filters.</p>'
    content+=section("Downloads and implementation")
    downloads=[(SOURCE/"fluid_comparison/comparison.csv","Eight-case comparison","CSV"),(run.HERE/"case_template.yaml","Shared case definition","YAML"),(SOURCE/"environment_and_sources.json","Native backend and source fingerprints","JSON"),(run.SCRIPTS/"study.py","Run orchestration","Python"),(Path(run.th.__file__).parent/"optimize_cycle.py","ThermoOpt native optimization implementation","Python"),(Path(run.th.__file__).parent/"cycles/cycle_power_recuperated.py","ThermoOpt cycle equations","Python"),(Path(run.th.__file__).parent/"components/basic_components.py","ThermoOpt component models","Python")]
    downloads.extend((SOURCE/f"sensitivity_{f}"/"sweep_results.csv",LABELS[f]+" complete sensitivity","CSV") for f in run.FAMILIES)
    content+='<div class="panel"><ul class="download-list">'+"".join(f'<li>{link(p,label,True)}<span>{kind}</span></li>' for p,label,kind in downloads)+"</ul></div>"
    content+='<p class="table-note">The native source fingerprints cover the optimization manager, cycle implementations and component models. Input files retain initial guesses; optimum files retain the selected final state.</p>'
    page("validation","Validation & data","Evidence and reproducibility","Verify the model that ran, understand the numerical limits and download the full optimization records.",content)


def build(preview=False):
    global SOURCE, OUTPUT, PREVIEW, TABLES, AUDITS, COMPARISON
    PREVIEW=preview
    SOURCE=(run.ARCHIVE/"coolprop_evaluator_2026-11-02"/"results") if preview else run.RESULTS
    OUTPUT=run.OUTPUT/("report_preview" if preview else "report")
    OUTPUT.mkdir(parents=True,exist_ok=True)
    assets=OUTPUT/"assets";assets.mkdir(exist_ok=True)
    for file in ("report.css","report.js"):
        shutil.copy2(WEB_ASSETS/file,assets/file)
    TABLES={f:pd.read_csv(SOURCE/f"sensitivity_{f}"/"sweep_results.csv") for f in run.FAMILIES}
    AUDITS={f:json.loads((SOURCE/f"sensitivity_{f}"/"validation_summary.json").read_text(encoding="utf-8")) for f in run.FAMILIES}
    COMPARISON=pd.read_csv(SOURCE/"fluid_comparison/comparison.csv")
    if not preview:
        for t in TABLES.values():
            assert t.loc[t.success,"optimization_backend"].eq("thermopt.native").all()
            assert t.loc[t.success,"optimization_called"].all()
    PAGES.clear()
    summary_figures()
    overview();comparison_page()
    for f in run.FAMILIES:fluid_page(f)
    methods_page();validation_page()
    (OUTPUT/"pages.json").write_text(json.dumps(PAGES,indent=2),encoding="utf-8")
    if not preview:
        import pack_standalone_report
        pack_standalone_report.build(OUTPUT,run.HERE/"report_partial_evaporation.html")
    print(f"Built {'preview' if preview else 'native'} report: {len(PAGES)} pages at {OUTPUT}")


if __name__=="__main__":
    build(preview="--preview-archive" in sys.argv)
