# Thermopt regression suite

## JaxProp 0.6 compatibility checkpoint

The current code requires the released JaxProp 0.6 package. A fresh locked
installation passes **262 tests**, with four unchanged strict expected failures.
See [the release-adoption record](../benchmark/RELEASE_060.md) and
[phase plan](../z_JAXPROP_060_IMPLEMENTATION_PLAN.md).

The earlier prepared-wheel checkpoint passed 253 tests, including three complete
optimization solves and five native-bicubic integration tests. Nine additional
deterministic tests now verify benchmark checkpoint eligibility, ratios and
invalid measurements without wall-clock assertions. See
[the native-backend workflow](../development/bicubic/README.md) for interpolation
budgets and [benchmark/](../benchmark/README.md) for repeatable phase checks.

The original `cycle_baselines.json` physical designs/states/energy/objective
references and their tolerances remain unchanged by this migration. New source
reference coordinates affect initialization; those expected designs are kept
separately in `initialization_060.json`. Corrected generalized-quality reporting
outside the two-phase region and two initial R134a superheating-constraint
residuals are kept in `source_annotations_060.json`. Those 38 values were
independently checked through fresh JaxProp source calls and saturation endpoint
formulas before acceptance. `read_baselines()` applies only these explicit
versioned reporting corrections; it does not relax physical comparisons or
regenerate references while testing.

Tests use the new keyword backend spelling, explicit mask policies in already
existing JAX-transformation tests, and public state serialization for mocked
quality fields. This compatibility work does not port the cycle solver to JAX.

This suite established behavior before changing fluid-property access. The jaxprop-only migration and functional refactor were checked against those original numerical references. Current compressor and nondimensional PTES references include the deliberate Figure 4 model correction described below. JAX cycle compilation is a separate future change.

## Run the tests

Use the repository's existing development environment (`poetry install`, or its `.venv` interpreter):

```powershell
# Everything, including three complete optimizations
.\.venv\Scripts\python.exe -m pytest

# Routine development checks, excluding complete solves
.\.venv\Scripts\python.exe -m pytest -m "not slow"

# Just the complete cycle solves
.\.venv\Scripts\python.exe -m pytest -m optimization

# Just plots and graphical exports
.\.venv\Scripts\python.exe -m pytest -m plotting
```

The equivalent Poetry command is `poetry run python -m pytest`. `python tests/run_tests.py` also works from any directory and forwards arguments to pytest. `pytest.ini` restricts discovery to `tests/`, so development scripts and archived example programs are not executed as tests.

All plotting is headless (`Agg`). Files are written to isolated, automatically cleaned directories under the ignored `.test_artifacts/`; example outputs are never overwritten. The custom `tmp_path` fixture uses ordinary directory permissions because private mode-0700 temporary directories can be inaccessible in managed Windows sessions. Pytest's optional cache is disabled for the same environment. No network, REFPROP license, GUI interaction or external optimization service is needed.

Expect roughly half a minute for the fast suite and another one to two minutes for the complete solves on a similar CPU. Runtime is not an assertion and depends on machine load.

## Coverage

| Module | Checks |
|---|---|
| `test_example_cases.py` | Numerical regression at initial and nearby designs for 14 cases, all seven distinct topology implementations and their aliases; energy/heat/shaft balances; repeated-call isolation; invalid topology; explicit two-phase expansion. |
| `test_optimization.py` | Three complete SLSQP solves from original initial designs, independent final reevaluation, normalized bound feasibility, every equality/inequality residual, reference objective/COP/efficiency and fixed-design state regression at saved optima. |
| `test_components.py` | Isentropic and polytropic compression/expansion, reversible limits, real-fluid CO₂ regression, gas/flashing valves, exchanger mapping/zero duty/energy, live nondimensional compressor/turbine outlet and design calculations, friction correlation. |
| `test_heat_exchanger_saturation.py` | Existing snapping tests: evaporation/condensation, pressure losses, counter/co-current layouts, endpoints, coarse competing crossings, dense pinch reference and constraint counts. |
| `test_heat_exchanger_plotting.py` | Existing nonuniform heat-duty alignment and artist updates. |
| `test_plotting.py` | Real bounded-resolution phase diagrams; T-s and p-h plots, logarithmic axes, pinch alignment, design updates without duplicate artists, reopening figures, stream visibility, scaling helpers, PNG/SVG/PDF rendering. |
| `test_workflow.py` | Manager setters/reloads, constraint changes, current-design YAML roundtrip, Excel readback of all endpoint properties and energy fields, iteration configuration/plot/report callbacks. |
| `test_configuration_and_utilities.py` | Configuration expressions/types/errors, nested state aliases and arrays, objective signs/scalar requirements, constraint signs/normalization/order/errors, serialization. |
| `test_math.py` | Smooth extrema, absolute values and sigmoids: broadcasting, symmetry, analytical values, approximation bounds, transition slopes and numerical stability. |
| `test_property_boundary.py` | Project-wide source check for direct backend imports/raw-state escapes and direct dependency declarations. Includes executable archives when present. |
| `test_machine_property_access.py` | Explicit fluid forwarding, inlet reuse, enthalpy-only flashes, rotor-static h,s state and entropy allocation. |
| `test_study_property_access.py` | Study configuration/property calculations in isolated interpreters and cached-evaluator/native residual parity. |
| `test_turbomachinery_correlations.py` | Independent paper equations, frozen Air/CO2 designs, 85 frozen demo points and physical scaling/loss trends. |
| `test_turbomachinery_functional.py` | JIT/gradient/batching comparisons, shrouded stress endpoints, equal-temperature limits, input validation and solver errors. |
| `test_paper_figures.py` | Independent PDF-curve comparisons for component Figures 1-7, visible marker counts/provenance, corrected Figure 4 correlation and documented Equation 23 error, Figure 3 thermodynamic reconstruction, affine extraction geometry and all PNG/SVG exports. |

Frozen cycle cases include polytropic simple CO₂, recuperated/recompression CO₂, steam Rankine, subcritical/partial-evaporation/transcritical organic cycles, an explicit wet-inlet organic expansion, air/R134a refrigeration, CO₂ heat-pump valve/turbine variants, and both PTES models including nondimensional machinery.

The original `test_example_cases.py` launched an interactive example program with reporting and animation side effects in its example directory. Its replacement exercises numerical cycle behavior directly, while dedicated tests cover optimization, plotting and exports independently.

## What the numerical references mean

`data/cases/` contains frozen YAML inputs copied from repository examples; ordinary example edits do not silently change a test's input. `orc_wet_expansion.yaml` is an intentional additional variant with inlet quality $0.7$ at the frozen expansion pressure. The source PEORC example's current initial design is dry, so its name alone did not provide wet-expansion coverage.

`data/cycle_baselines.json` records:

- Design ordering, initial and fixed nearby design vectors.
- The objective and **all** equality/inequality residuals, preserving order and signs.
- Every energy-analysis metric, component endpoint $p,h,T,s,\rho,Q$, mass flows and work/duty quantities.
- Complete exchanger temperature profiles and pinch-difference arrays.
- Nondimensional machinery scalar design/performance outputs.

The references were captured before production changes using thermopt 0.2.5, jaxprop 0.5.7/HEOS, CoolProp 7.2.0, NumPy 2.4.4, SciPy 1.17.1 and pysolver_view 0.7.1. The JSON records versions and the thermopt commit. It does not claim every original initial design is feasible or optimal. Physics tests complement these characterization records rather than treating existing output as proof of correctness.

`assert_snapshot` in `regression_support.py` uses a relative tolerance of $2\times10^{-7}$ and field-specific absolute allowances, including $0.05$ Pa, $0.02$ J/kg, $2\times10^{-5}$ K, $2\times10^{-4}$ J/(kg K), $10^{-6}$ kg/s, and $0.1$ W for near-zero energy balances. Dimensionless metrics retain an absolute allowance of $10^{-7}$. The initial design-vector check is tighter to catch changing special-point initialization.

Optimization references use explicit SLSQP options shared with the tests. Tests require solver success and normalized constraint violations below $2\times10^{-6}$, then compare objective-related metrics with relative tolerance $5\times10^{-5}$ and absolute tolerance $2\times10^{-6}$. They do not require identical iteration counts, inactive design variables or convergence histories: multiple designs can represent equivalent optima. Initial guesses may be infeasible, so an objective need not improve monotonically while feasibility is restored.

These are current-HEOS regression tolerances. A later bicubic/perfect-gas/JAX migration needs an explicit model-specific error policy; do not broaden all tolerances merely to make a new backend pass. Perfect-gas approximation error and bicubic interpolation error are separate from implementation regression.

## Known defects recorded as expected failures

The following tests use `xfail(strict=True)` with the expected exception class. A future fix produces an unexpected pass and fails the suite until its marker is removed. Other tests remain active; there are no broad component/topology skips.

| Test | Existing defect |
|---|---|
| `test_cocurrent_total_heat_is_nonzero_and_balanced` | Total exchanger heat accounting takes the first hot-side array entry, which is zero for co-current ordering. |
| `test_recompression_reported_energy_balance_closes` | Reported balance includes external auxiliary pump work against working-fluid heat duties; the nonzero report equals the auxiliary powers. Working-loop and individual exchanger balances are checked separately. |
| `test_documented_logcosh_absolute_selector_is_accepted` | `smooth_abs` documents `logcosh` but accepts `logarithmic` instead. |
| `test_logarithmic_absolute_remains_finite_at_default_smoothing_scale` | Direct `cosh(x/epsilon)` overflows for ordinary arguments at the default smoothing width. |

The turbine coefficient test follows the published Figure 6 and author implementation (0.50). Table 6's 0.30 is a documented typo; these four checks pass, including speeds away from the default.

The expected Agg warning from the existing `plt.pause` call may appear during plotting tests. It does not require interactive display.

## Updating references deliberately

The approved compressor velocity correction uses $w/U=0.38+3\phi-7.8\phi^2$ from Figure 4, replacing the erroneous printed Equation 23. References were deliberately refreshed for the 50 compressor sweep designs, two high-pressure compressor designs, the default component outlet, and the initial/nearby `ptes_turbo` cycle evaluations. Model-correction metadata records this choice. Turbine-only references, the other 13 cycle cases, and optimization references are unchanged. Existing tolerances and independently extracted publication curves are retained. A dedicated test checks the corrected velocity in both actual and reference Reynolds numbers.

Tests **never** regenerate their expected outputs. The recording tools protect existing files unless `--overwrite` is explicitly supplied:

```powershell
.\.venv\Scripts\python.exe tests\record_baselines.py --overwrite
.\.venv\Scripts\python.exe tests\record_optimization_baselines.py --overwrite
```

Inspect diffs and version changes before accepting new references. The cycle recorder preserves the previous reference and exits with an error if any case cannot be captured. The optimization recorder writes completed solves incrementally so an interrupted later solve does not lose prior work; it is a developer recording tool, not part of pytest. Keep fixture changes separate from implementation changes when reviewing a migration.

## Property-access audit and remaining scope

[PROPERTY_ACCESS_AUDIT.md](PROPERTY_ACCESS_AUDIT.md) retains the pre-migration inventory. Executable property accesses now use jaxprop, including the nondimensional machinery, studies and benchmarks. CoolProp remains an implementation dependency of jaxprop's HEOS backend. Historical outputs and recorded backend/version descriptions are preserved.

The machinery API is functional: `solve_compressor(fluid, inlet, p_out, mass_flow, design)` and `solve_turbine(...)` use the supplied jaxprop model, reuse the inlet and solve only in outlet enthalpy. The shrouded stress-key defect is fixed and now has passing interpolation/endpoint checks. Separate performance functions support JAX arithmetic; `test_turbomachinery_functional.py` checks JIT, derivatives, batching, equal-temperature limits and invalid inputs. The SciPy root solvers remain host operations. See the [functional API and JAX boundary](../development/turbomachinery/README.md).

This is a regression foundation, not a claim of 100% line coverage or thermodynamic certification. It does not test full optimization for every topology, REFPROP/IPOPT/SNOPT, interactive live-update loops, video encoders/animations, every plotting aesthetic, every correlation option or optional JAX backends. The first backend changes should run the complete suite; dedicated cross-backend and derivative checks can then be added without changing these HEOS references.
