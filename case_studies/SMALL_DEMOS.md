# Small demonstration scripts: execution check

All 14 scripts passed on 16 November 2026, including all 15 optimized cycle
runs. This profile runs the actual one/two-case launchers
with their configured iteration limits, normal reporting callbacks and final
exports. It covers 14 scripts containing 15 cycle optimizations and two
standalone property/correlation demonstrations. The machine-readable
[execution inventory](SMALL_DEMOS.json) records each final result, iterations,
objective, constraint violation, source hash and output location.

| Script (relative to this directory) | Optimizations | Configuration/profile |
| --- | ---: | --- |
| `sCO2_recuperated/run_optimization.py` | 1 | Recuperated CO2, CoolProp HEOS |
| `sCO2_recompression/run_optimization.py` | 1 | Recompression CO2, native bicubic 200x200 |
| `rankine_steam/run_optimization.py` | 1 | Water Rankine cycle |
| `rankine_organic/run_optimization.py` | 1 | Default butane ORC; commented alternatives excluded |
| `refrigeration_R134a/run_optimization.py` | 1 | R134a refrigeration |
| `refrigeration_air/run_optimization.py` | 1 | Recuperated air refrigeration |
| `heat_pump_CO2/run_optimization.py` | 1 | Default Joule-Thomson valve case |
| `heat_pump_CO2/run_comparison.py` | 2 | Valve and two-phase turbine |
| `PTES_CO2/case_baseline/run_optimization.py` | 1 | Fixed-efficiency PTES |
| `PTES_CO2/case_ideal/run_optimization.py` | 1 | Ideal PTES, including its experimental setup |
| `PTES_CO2/case_shaft_constraint/run_optimization.py` | 2 | Unconstrained and shaft-constrained PTES |
| `PTES_CO2/case_turbo_model/run_optimization.py` | 2 | Turbomachinery and matched fixed-efficiency PTES |
| `HX_development/blended_calculation.py` | 0 | 300-state blended property plot |
| `HX_development/PHE_correlations.py` | 0 | Water-state evaluation and correlation plots |

The recompression run generated a new matching table through JaxProp; it did
not substitute the 128x128 benchmark configuration for the ordinary YAML.

Run a demo from its own directory. For example, from the repository root:

```powershell
$demoPython = (Resolve-Path .\.venv\Scripts\python.exe).Path
Set-Location .\case_studies\sCO2_recompression
& $demoPython .\run_optimization.py
```

For a headless run with a fresh output location, set these before invoking the
same launcher:

```powershell
$env:MPLBACKEND = "Agg"
$env:DISABLE_PLOTS = "1"
$env:THERMOPT_RESULTS_DIR = Join-Path (Get-Location) "results\manual_check"
& $demoPython .\run_optimization.py
```

Use a new directory for a new PTES check: the two-case PTES launchers otherwise
intentionally reload their existing `simulation_data.pkl`. Verification outputs
are under each case's `results/demo_checks/<timestamp>/`. Property-table caches
remain under the case's configured `results/property_tables/` directory.

The verification used Matplotlib Agg, JAX CPU/x64 and two independent worker
processes. It disabled interactive file watchers, blocking windows and optional
movie encoding. Normal report, plot, configuration and Excel callbacks remained
enabled as configured; final saves and the launchers' pickle exports executed.
Remaining figures from plotting-only demos were saved by the verification
worker. Interactive GUI behavior and movie encoding are outside this check.

Each optimized cycle had to report native solver success, return a finite final
objective/residual vector and satisfy the configured native feasibility and
first-order checks. The inventory retains the measured final constraint
violation. Script exit alone was insufficient: initial runs which stopped with
a failed first-order check were recorded as unconverged. No physical reference
fixtures, equations, bounds or initial designs were changed.

Minimal repairs found by execution:

- The plate heat-exchanger demonstration now constructs `jaxprop.Fluid`
  directly and uses `jp.PT_INPUTS`; the removed `thermopt.Fluid` export failed.
- The heat-pump comparison uses canonical property labels and native state
  identifier metadata. This fixes the `KeyError: 's'` and restores the intended
  process selection. Its legend now follows the actual configured cycle colors.
- Launchers accept a fresh results directory and headless execution; the PTES
  baseline and ideal figures now stay inside `results/`.
- R134a uses three-point central differences. Two-point derivatives in the
  actual variable order stopped with first-order residual 0.209, despite small
  steps. Central differences converge in 23 iterations to COP 2.3928322, in
  agreement with the frozen optimum. Report/plot callbacks were independently
  excluded as the cause of the failure.
- The turbomachinery PTES YAML uses native stopping tolerance `1e-8` instead
  of `1e-6`. Its original feasibility gate `1e-6`, first-order gate `1e-4` and
  two-point derivatives remain unchanged. Central differences alone did not
  resolve its premature stopping. The tighter native stopping continues beyond
  the former 57-iteration point; the unconstrained solve converges in 82
  iterations with round-trip efficiency 0.4490831.
- The default-solver regression test now preserves YAML variable order during
  its round trip. Alphabetical sorting previously changed the numerical path
  and concealed the R134a failure. Explicit numerical overrides are reviewed
  separately from the immutable physical-input hashes.

The full regression suite passes **564 tests with four unchanged strict expected
failures**. Source/bicubic optimization and heat-pump comparison figures, plus
the PHE correlation output, were visually inspected. These are execution and
numerical checks; no between-version timing comparison is claimed.

The environment uses PySolverView 0.8.2, Pounce 0.13.0, JAX/JAXlib 0.6.2 and
SciPy 1.17.1. All 43 installed JaxProp Python package files match the official
PyPI 0.6.0 wheel, SHA256
`0a8dfb68ca4d21cf97b7483e1da776b946cec108a94d50ebb770bb9848127d3c`;
see [the release assessment](../benchmark/RELEASE_060.md). The separate local
metadata-fix candidate used in earlier performance work was not used for these
demo runs.

The following programs remain outside this small-demo profile: the six-way
solver/backend comparisons, the ten-case HX node/mesh comparison, sensitivity
sweeps, full Rankine research/report pipelines, historical figure preparation
and saved-input table exporters. Larger studies and output round-trip coverage
remain open in Phases 6/7 of the implementation checklist.
