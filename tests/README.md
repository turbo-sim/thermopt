# Thermopt regression suite

This suite establishes current behavior **before** changing fluid-property access or introducing JAX cycle evaluation. It extends the existing `tests/` directory. Production models and dependencies are unchanged.

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

## Five existing defects recorded as expected failures

The following tests use `xfail(strict=True)` with the expected exception class. A future fix produces an unexpected pass and fails the suite until its marker is removed. Other tests remain active; there are no broad component/topology skips.

| Test | Existing defect |
|---|---|
| `test_cocurrent_total_heat_is_nonzero_and_balanced` | Total exchanger heat accounting takes the first hot-side array entry, which is zero for co-current ordering. |
| `test_shrouded_compressor_can_evaluate_an_outlet` | Shrouded-compressor stress calculations reference missing `K_stress_shrouded_phi1/phi2` reference keys. |
| `test_recompression_reported_energy_balance_closes` | Reported balance includes external auxiliary pump work against working-fluid heat duties; the nonzero report equals the auxiliary powers. Working-loop and individual exchanger balances are checked separately. |
| `test_documented_logcosh_absolute_selector_is_accepted` | `smooth_abs` documents `logcosh` but accepts `logarithmic` instead. |
| `test_logarithmic_absolute_remains_finite_at_default_smoothing_scale` | Direct `cosh(x/epsilon)` overflows for ordinary arguments at the default smoothing width. |

The expected Agg warning from the existing `plt.pause` call may appear during plotting tests. It does not require interactive display.

## Updating references deliberately

Tests **never** regenerate their expected outputs. The recording tools protect existing files unless `--overwrite` is explicitly supplied:

```powershell
.\.venv\Scripts\python.exe tests\record_baselines.py --overwrite
.\.venv\Scripts\python.exe tests\record_optimization_baselines.py --overwrite
```

Inspect diffs and version changes before accepting new references. The cycle recorder preserves the previous reference and exits with an error if any case cannot be captured. The optimization recorder writes completed solves incrementally so an interrupted later solve does not lose prior work; it is a developer recording tool, not part of pytest. Keep fixture changes separate from implementation changes when reviewing a migration.

## Property-access audit and remaining scope

[PROPERTY_ACCESS_AUDIT.md](PROPERTY_ACCESS_AUDIT.md) inventories direct imports, indirect bypasses, obsolete helpers and dependency declarations across the repository. Cleanup is deferred as requested. New test property evaluations use jaxprop; production nondimensional machinery still performs its existing internal CoolProp calls while being tested. CoolProp remains an implementation dependency of jaxprop's HEOS backend, even when thermopt imports only jaxprop.

This is a regression foundation, not a claim of 100% line coverage or thermodynamic certification. It does not test full optimization for every topology, REFPROP/IPOPT/SNOPT, interactive live-update loops, video encoders/animations, every plotting aesthetic, every correlation option or optional JAX backends. The first backend changes should run the complete suite; dedicated cross-backend and derivative checks can then be added without changing these HEOS references.
