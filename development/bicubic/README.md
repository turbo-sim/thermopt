# JaxProp 0.6 compatibility and native bicubic cases

This workflow keeps ThermOpt's eager NumPy/SciPy implementation. It selects
actual JaxProp evaluators through YAML; there is no application fluid adapter,
replacement thermodynamic solver, outer JIT or silent source fallback.

The completed [verification report](VERIFICATION_060.md) and
[portable measurements](verification_060.json) record the prepared-wheel tests,
optimized designs and synchronized timings.

## Select the property model

Keep the existing example's equations, solver and design bounds, and replace
its `working_fluid` section with:

```yaml
working_fluid:
  name: CO2
  method: bicubic
  backend: HEOS
  table_options:
    h_min: 110000
    h_max: 1300000
    p_min: 8000000
    p_max: 39000000
    N_h: 200
    N_p: 200
    N_p_sat: 8
    coarse_step: 1
    gradient_method: central
    table_dir: development/bicubic/results/tables
```

ThermOpt passes these settings to `jaxprop.Fluid`. JaxProp generates or validates
and loads the table. Components continue calling `fluid.get_state(pair, a, b)`.
Each optimization problem owns reusable configured models, including separate
working/heating/cooling stream identities. Changing `method` to `coolprop` and
removing `table_options` restores source evaluation.

The 0.6 interface uses keyword `backend`; it no longer accepts the old
`Fluid("Air", "HEOS")` positional spelling. ThermOpt uses the constructor
vocabulary directly: constructor-level `exceptions` is rejected even when true.
Normal cycle calls use JaxProp's raising policy. An explicit `errors="mask"` or
`errors="clip"` belongs on an individual state call only where its consuming
algorithm has defined behavior for failed or adjusted states.

Phase 2 migrates 36 maintained YAML inputs (108 fluid declarations), preserving
all physical values, constraints and solver options. Historical initializer and
figure-result YAML files remain dated data. The current wet-to-dry study reads
only design-variable values from its 484 saved initializers and uses the new
template for actual fluid settings. Old standalone configurations must be
migrated before replay; no runtime legacy translation is provided.

Table bounds must cover every physical state and finite-difference/optimizer
probe. The CO2 rectangle above covers the validated supercritical cases, not
every possible CO2 cycle. Failure/out-of-domain requests raise; they do not
clip, extrapolate or switch backends. Setup-only saturation/dilute reference
queries and source-EOS phase plots use a separate configured JaxProp source,
so they need not lie within the numerical table rectangle.

Optional source diagnostics are requested only on models supporting them.
Unsupported/disabled optional spreadsheet fields are blank. A constraint that
actually needs an unsupported property still fails explicitly. Saturation-node
snapping requires a saturation-capable model; ordinary equilibrium bicubic
models reject it. Phase-changing Rankine/ORC cases remain verified with the
CoolProp route; successful single-phase cases do not establish universal
interchangeability around phase boundaries.

## Run the comparisons

Install the released 0.6.0 package into the environment used for ThermOpt.
Confirm `jaxprop.__version__` and `jaxprop.__file__` before running. The package
dependency requires `^0.6.0` and the refreshed lock selects 0.6.0. The earlier
prepared-wheel results remain historical evidence; the
[migration plan](../../z_JAXPROP_060_IMPLEMENTATION_PLAN.md) records verification
against the actual PyPI artifact and a fresh locked environment.

```powershell
python -c "import jaxprop; print(jaxprop.__version__, jaxprop.__file__)"
python -m pytest
python development/bicubic/validate_060.py --case sco2_recuperated --grids 100 200 --optimize
python development/bicubic/validate_060.py --case sco2_recuperated --grids 200 --all-fluids --optimize --output development/bicubic/results/compatibility_060/native_all --table-directory development/bicubic/results/compatibility_060/native/tables
python development/bicubic/validate_060.py --case sco2_simple_polytropic --grids 200
python development/bicubic/validate_060.py --case sco2_recompression --grids 200
python development/bicubic/validate_060.py --case refrigeration_air --grids 100 --all-fluids
```

The runner saves ordinary complete `config.yaml` files under each case/backend
directory. They run with `ThermodynamicCycleOptimization` directly and can be
used by the existing example launcher/`--config` option. Reported solve settings
retain each example's own stopping tolerances and iteration budgets; only
interactive/reporting callbacks are disabled for the comparison.

Source and table cases use identical physical designs and source-based setup
references. The interpolation budgets are fixed separately: maximum 0.3 K
temperature/profile/pinch error, 2% endpoint density error, 0.02 normalized
residual error and 0.0005 absolute efficiency/COP error. These allowances do
not replace the much tighter source regression tolerances. Optimization results
are reevaluated with both models at the returned design, preserving numerical
and source feasibility differences in the report. Timings alternate complete
eager source/table evaluations with explicit synchronization; construction/table
generation and complete solve costs are reported separately.

The old `bicubic_adapter.py`, `compare_backends.py` and earlier results describe
the historical 0.5 experiment. Use `validate_060.py` for the native 0.6 workflow;
it imports neither the old adapter nor its monkey patches.

## Implementation and verification checklist

- [x] Preserve the existing dirty worktree and record its initial patch/source
  hashes before the compatibility changes. Leave JaxProp and Barotropy unchanged.
- [x] Install the prepared 0.6 wheel in an isolated dependency directory outside
  ThermOpt, leaving its old `.venv` JaxProp installation intact. Record the
  actual module path, wheel hash and dependencies with the verification report.
- [x] Migrate construction, explicit source reference setup and source-EOS
  plotting; reuse real model objects during repeated cycle evaluations.
- [x] Request optional component diagnostics according to model capabilities,
  preserve scalar ODE state shapes, and export unavailable optional fields
  honestly. Retain current NumPy/SciPy equations and solver tolerances.
- [x] Pass all 14 frozen cycle cases at initial/nearby physical designs and the
  existing optimization/plot/export/functional-machine suite: 253 passed and
  four unchanged expected failures, including five new bicubic integration tests.
- [x] Retain the original physical fixtures and comparison tolerances. Record
  corrected 0.6 initialization separately in `tests/data/initialization_060.json`.
  Record 36 corrected generalized-quality values and two R134a superheating
  residuals in `tests/data/source_annotations_060.json`, verified independently
  from fresh source calls/saturation endpoints. No objective or physical-state
  reference was replaced.
- [x] Complete native 100/200-grid recuperated-CO2 optimization comparisons,
  including same-design source cross-evaluation and all-stream tables.
- [x] Complete additional simple/polytropic/recompression/refrigeration native
  table comparisons and record numerical agreement and eager execution costs.
- [x] Verify native bicubic phase/cycle/pinch plotting and numeric Excel cells,
  including blank unsupported optional fields. Correct JAX scalar values being
  written as text at the export boundary; no new work enters property calls.
- [x] Normalize YAML's NumPy scalar setup options to Python metadata once,
  avoiding static-array warnings in configured native models. Direct raw JaxProp
  construction with NumPy static options still emits a harmless upstream warning
  in two component tests; central library normalization is an optional follow-up.
- [x] Refresh `poetry.lock` against published JaxProp 0.6 and verify a fresh
  locked main-dependency installation. Only JaxProp changes version, from 0.5.7
  to 0.6.0. Correction to the earlier audit: the incoming lock contained 0.5.7;
  it was outdated, rather than lacking the JaxProp package.
- [ ] Optional follow-up: migrate/archive the old exploratory `HX_development`
  scripts still referring to `thermopt.Fluid`/`states_to_dict`; these are outside
  the tested case/optimization workflows. The archived recompression figure
  script also has a separate Barotropy dependency and external result inputs.
- [ ] Existing independent defects remain tracked by strict expected failures:
  co-current total heat reporting, recompression's reported auxiliary-power
  balance, the documented logcosh selector, and overflow in that smoothing form.

Full raw logs/configurations/tables are ignored under `results/compatibility_060`.
The compact final verification record summarizes the retained results and
their scope without treating unexecuted scripts as passing examples.
