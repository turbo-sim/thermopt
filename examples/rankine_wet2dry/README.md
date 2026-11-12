# Partial evaporation: updated EMPOWER study

Simple n-butane, n-pentane and cyclopentane cycles are compared with a recuperated
Novec 649 cycle. The boundary conditions follow Nishith's 20 October 2026 request
to adopt Cristina's clarification. Her revised workbook supplies the 88.6%
turbine efficiency used for the standard reference and equal-efficiency comparison.

## Directory layout

```text
rankine_wet2dry/
  case_template.yaml                 Shared case definition
  README.md
  inputs/initializers/              Versioned warm-start states
  reference/v2/                     Authoritative v2 comparison tables
  reference/reproduction_summary.json Recorded comparison of the consolidated run
  report_partial_evaporation.html    Standalone shareable report
  report_partial_evaporation.md      Linear Markdown reference
  report_partial_evaporation.pdf     Standalone PDF reference
  scripts/                          Run and report-generation code
    run_all.py
    compare_results.py
    study.py
    refine_results.py
    validate_saved_results.py
    run_sensitivity_efficiency_and_exploitation.py
    run_fluid_comparison.py
    build_report.py
    build_html_report.py
    pack_standalone_report.py
    render_report_pdf.py
    paths.py
    web_assets/                     CSS and JavaScript sources
  output/                           Generated/local files; ignored by Git
    results/                        Native cases, tables, figures and JSON
    logs/                           Run and report-generation logs
    report/                         Intermediate multipage HTML build
    figures.html                    Local figure gallery
    reference_inputs/               Supplied email and workbooks
    archive/                        Earlier dataset and historical snapshots
    cache/                          Local caches
    layout_migration.json           Preserved-record hashes and path mapping
    standalone_manifest.json        Embedded report asset inventory
```

The case definition, code, initializers, v2 reference tables, README and final
reports can be committed. Everything
under `output/` is excluded by this directory's `.gitignore`.

## Sharing the report

**Send only `report_partial_evaporation.html`.** It contains all eight logical
pages, CSS, JavaScript, vector figures, result tables and downloadable case inputs,
optima and solver histories. Navigation, filters, figure enlargement and downloads
work without a server or companion folders. The file makes no network requests.

The navigation keeps the report divided into an overview, a comparison page, four
fluid pages, methods, and validation/data. Each fluid page has tabs for performance,
recovery, phases and cycle evolution. Related figures remain side by side on wide
screens and stack on narrow screens.

The heat-recovery tabs factor the power gain into recovered source heat and net
conversion of that heat. Novec's tab also reports optimized recuperator
effectiveness. Cycle-evolution tabs compare admission pressure and temperature at
equal turbine efficiency, and the methods page lists the optimization variables.

Each of the eight comparison cases includes `optimal_solution.xlsx` with native
cycle states and energy analysis, plus `design_variable_bounds.csv` with lower and
upper bounds, optimum values and the fractional distance to the nearest bound.
Bounds within 0.01% of the variable range are flagged. The comparison page embeds
these downloads in the standalone HTML.

The standalone document offers SVG downloads and printing of the current logical
page. PNG/PDF/individual-figure HTML exports remain in `output/results/`; the full
PDF report is also available as the separate top-level PDF file. Original email
and workbook attachments remain local under `output/reference_inputs/`; their
citations are included in the shareable report.

The Markdown reference uses local images under `output/results/`. Use the HTML or
PDF when sharing one file.

## Running

From the repository root, use the ThermoOpt environment. On Windows the local
environment provides the following command. An activated equivalent environment
can use `python` in its place.

```console
.venv/Scripts/python.exe examples/rankine_wet2dry/scripts/run_all.py
```

To rebuild only the reports from existing tables and figures:

```console
.venv/Scripts/python.exe examples/rankine_wet2dry/scripts/run_all.py --report-only
```

Paths are resolved relative to the study, so individual scripts also work from
another current directory. Console output from the full pipeline is saved under
`output/logs/`. Completed cases are reused only when their physical configuration,
native backend and source-code signature match. Use a separate results directory
for changed physical assumptions.

The full pipeline runs optimization, targeted recovery, validation, Matplotlib
figure generation, the equal-efficiency comparison, comparison against the saved
v2 tables, report assembly and PDF
rendering. PDF generation requires `reportlab` and Pillow; the installed Codex
runtime is used as a fallback when these are absent from the active environment.

To compare existing output with the versioned v2 tables independently, run

```console
.venv/Scripts/python.exe examples/rankine_wet2dry/scripts/compare_results.py
```

The comparator checks all 484 finite-heat sweep cases, 20 unavailable endpoints
and eight comparison cases. It requires matching case grids and metadata, then
checks each physical metric using `abs(rerun - reference) <= atol + rtol * abs(reference)`.
Efficiencies, powers and source mass flow use 0.01% relative tolerance. Pressures
use 0.1% relative tolerance, temperatures use 0.02 K absolute tolerance and raw
qualities use 0.001 absolute tolerance. The complete metric tolerances and observed
differences are saved in `output/results/reproducibility_comparison.json` and
`.csv`. A failed comparison stops the full pipeline before report generation.
`--reference-dir` and `--results-dir` select other table locations. Run
`compare_results.py --help` for these path options.

## Applied conditions

| Parameter | Setting |
|---|---|
| Source | Water, 140 °C, 10 bar |
| Minimum source outlet | 140, 135, ..., 20 °C; 70 °C for the eight-case comparison |
| Cooling water | 15 °C inlet; exactly 25 °C return |
| Heater / condenser approach | 4.5 K at 100 locations each |
| Novec recuperator | 5 K approach; 25 locations; effectiveness optimized |
| Standard turbine reference | 88.6% isentropic efficiency throughout |
| Partial turbine sensitivity | 90%, 80%, 70%, 60% |
| Equal-efficiency comparison | Both modes at 88.6%; 70 °C source-outlet limit |
| Standard inlet / outlet | 2 K inlet superheat; dry outlet |
| Partial inlet / outlet | Optimized inlet; no dry-outlet requirement |
| Pumps / pressure losses | 70% pump efficiencies; 1% per active exchanger side |
| Sizing | 1 MW turbine output minus working-fluid pump; external pumps deducted |

There are **484 finite-heat optimized designs**, including four extra partial
cases at 88.6% for the main comparison. Twenty zero-heat endpoints at 140 °C are
unavailable. Internal recuperation, standard superheat, pump efficiencies and
pressure losses retain the previous study's assumptions where the exchange does
not specify replacements. Generator losses and EXERGY auxiliary-power accounting
are not fully specified, so model mechanical power and workbook electrical power
are distinguished.

## Native ThermoOpt backend and reproducibility

Every reported finite-heat case was optimized through the public
`thermopt.ThermodynamicCycleOptimization` entry point using the unmodified native
`ThermodynamicCycleProblem`. The active run uses no custom fitness override or
alternate cycle evaluator. Each attempt preserves a complete input `case.yaml`
that can be reproduced directly:

```python
import thermopt as th

cycle = th.ThermodynamicCycleOptimization("case.yaml", out_dir="verification_run")
cycle.run_optimization()
```

Versioned states under `inputs/initializers/` supply the original v2 initial
guesses. These inputs let the consolidated folder reproduce the study without
another example directory. Every finite-heat case is optimized with native
ThermoOpt. Tables under `reference/v2/` preserve the authoritative numerical
reference for comparing fresh runs.

Historical local datasets are retained under `output/archive/`. Native solver
success, physical flows, power, energy balances and constraints are checked for
the active results. The full pipeline compares all resulting physical outputs
with the versioned v2 reference before generating the reports.

`reference/v2/environment_and_sources.json` preserves the original study's
environment and source fingerprints. The current run writes its package versions,
input hashes and native model-source fingerprints to
`output/results/environment_and_sources.json`. The committed
`reference/reproduction_summary.json` records the consolidation verification.

## Normalization and numerical limits

Fixed-reference efficiency uses source-water heat from **140 to 15 °C** at 10 bar,
so it is proportional to net power per source mass flow. Nominal exploitation is
$100(140-T_{\min})/125$; actual exploitation uses the achieved outlet. A 20 °C limit
therefore represents 96% nominal exploitation.

ThermoOpt's native objective uses heat available down to each case's minimum
outlet. At a fixed limit these efficiency definitions select the same optimum.
Powers can be scaled to Cristina's 194.44 kg/s source flow without changing states.

The 4.5 K external approach is imposed at the optimization grid locations.
Selected designs are reevaluated at 501 locations; some between-grid approaches
are below target. The report documents this limitation. Precise break-even values,
small gains and close fluid rankings require refinement before design conclusions.
The optimization remains local and does not guarantee a global optimum.

The validation page reports native KKT diagnostics at five times the solver
tolerance alongside accepted-case counts. Acceptance uses native solver success,
physical flows, energy balance and constraint satisfaction. The stricter
first-order diagnostics provide additional information about each local solution.
