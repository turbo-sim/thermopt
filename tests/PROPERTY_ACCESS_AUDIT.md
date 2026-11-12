# Fluid-property access audit

Date: 12 November 2026. This inventory records the current repository before the requested property-access migration. No production import, fluid implementation or dependency was changed by this audit.

## Intended dependency boundary

Thermopt should declare `jaxprop` as its fluid-property dependency and obtain property states through `jaxprop.Fluid`, `FluidJAX`, `FluidPerfectGas` or `FluidBicubic`. Direct `CoolProp`/`coolpropx` use in thermopt, examples and development scripts remains to be migrated after the regression suite is established.

**Removing direct thermopt imports does not mean uninstalling CoolProp from the environment.** The installed and locked jaxprop 0.5.7 package declares CoolProp as a dependency; its HEOS-backed `Fluid` and callback-backed `FluidJAX` need that implementation internally. The intended change is ownership of the API and dependency declaration: thermopt depends on jaxprop, and jaxprop manages its backends.

## Installation declarations

| Location | Finding | Later action |
|---|---|---|
| `pyproject.toml:24` | Active direct dependency is `jaxprop = "^0.5.7"`; no active direct `CoolProp` or `coolpropx` dependency. | Preserve jaxprop as the direct dependency. |
| `pyproject.toml:27` | Commented development-path example for `coolpropx`. | Remove obsolete comment during the migration cleanup. |
| `poetry.lock:493` | Locked `coolprop` package, version 7.2.0. | Retain while required transitively; do not manually delete lock entries. |
| `poetry.lock:1003–1016` | Locked jaxprop 0.5.7 explicitly depends on `CoolProp >=7.0.0,<8.0.0`. | Expected transitive dependency. |
| `poetry.lock:2404` | `pysolver-view` lists an optional `thermodynamics` extra containing CoolProp. | This extra declaration alone does not establish that thermopt requests it. |
| `environment.yaml` | Declares Python, pip and Poetry; no direct fluid-property package declaration. | No fluid-dependency edit needed here. |
| `.github/workflows/` | Workflows install dependencies through Poetry; no separate CoolProp installation command found. | Keep dependency installation centralized. |

No additional explicit CoolProp/coolpropx installation command was found in the README or installation documentation.

## Production property bypasses

| Location | Current access | Migration consideration |
|---|---|---|
| `thermopt/components/basic_components.py:9` | `import CoolProp.CoolProp as cp`. | Replace the nondimensional branches' raw-state interface after baseline tests exist. |
| `thermopt/components/basic_components.py:361–369` | Compressor branch constructs `cp.AbstractState("HEOS", fluid.name)` and passes it to `set_CoolProp_fluid`. | Must preserve the supplied fluid/backend choice and compressor outputs. |
| `thermopt/components/basic_components.py:469–477` | Turbine branch constructs an equivalent raw state. | Must preserve turbine outlet, transport properties and model outputs. |
| `thermopt/components/turbomachinery_nondimensional.py:18` | Direct CoolProp import. | Replace constants and mutable state calls through a jaxprop-facing adapter/API. |
| `thermopt/components/turbomachinery_nondimensional.py:41–124` | `_Machine` uses raw `update`, `T`, `p`, `hmass`, `smass`, `rhomass`, `viscosity` and `speed_sound` calls to evaluate inlet, isentropic outlet and iterative outlet states. | Covers more than an import rename; regression tests should check both enthalpy- and temperature-based outlet solving. |

Cycle modules and `thermopt/optimize_cycle.py` already import `jaxprop` and instantiate its `Fluid` class. Most basic-component property calls already use the passed jaxprop fluid object. The direct property bypass in production is concentrated in nondimensional turbomachinery.

Existing docstrings in `basic_components.py` still name `CoolProp.AbstractState`/CoolProp states (for example lines 128, 138, 245 and 530). Update these type descriptions when the API is migrated; narrative mentions of a backend are not direct imports.

## Examples, tests and development code

| Location | Finding |
|---|---|
| `tests/utilities.py:4` | Obsolete `import coolpropx as cpx`; the import is not used in the helper body. |
| `tests/test_heat_exchanger_saturation.py:189` | Existing test reaches through `Fluid.abstract_state.p_critical()` rather than using jaxprop's critical-point state. Deferred with the other property-interface changes. |
| `examples/rankine_transcritical_vs_trilateral/run_sensitivity.py:21` | Direct CoolProp import; `PropsSI` calls initialize states, bounds, saturation departures and fluid-specific seeds (examples at lines 55–56, 72–106 and 149–208). |
| `examples/rankine_transcritical_vs_trilateral/fast_evaluator.py:9–17` | Gets `CP` from `run_sensitivity` and caches raw `AbstractState` objects. Further `PropsSI` calls evaluate exchanger temperatures and isentropic/saturation states (lines 33, 92–94 and 110–112). A scan for import statements alone would miss this bypass. |
| `examples/rankine_wet2dry/scripts/study.py:15` | Direct CoolProp import; `PropsSI` computes water heat availability and initial property bounds/states (lines 51–75). |
| `development/jax_assessment/benchmark_backends.py:27–28` | Direct CoolProp imports used for the prior JAX assessment's native-HEOS and native-CoolProp-bicubic benchmark references; raw states constructed at lines 137–138. This script is explicitly a backend comparison, but it is still outside the requested jaxprop-only access boundary. |
| `examples/rankine_wet2dry/output/archive/coolprop_evaluator_2026-11-02/study.py:14` | Archived direct CoolProp import and historical `PropsSI` evaluator. |
| `examples/rankine_wet2dry/output/archive/coolprop_evaluator_2026-11-02/fast_evaluator.py:7` | Archived direct import, raw `AbstractState` and `PropsSI` calls. |

The archived files were inventoried but left intact. Migration should explicitly decide whether historical executable snapshots are retained as historical evidence or updated along with active examples. Existing numerical benchmark JSON and report prose refer to CoolProp as the model used when they were generated; changing those historical descriptions would misrepresent the measurements.

The new regression tests access property states through jaxprop. They intentionally do not enforce a zero-CoolProp-import assertion yet: the current production bypasses are the behavior being captured before conversion.

## Scope and verification

The scan covered the entire thermopt repository's Python source, available notebook files, project/environment/lock declarations, workflows and text documentation, including untracked development code and Python files in archived outputs. No notebooks were present in the inspected checkout. `.git`, `.venv`, Python/pytest caches and generated binary artifacts were excluded; generated HTML and numerical output tables were not treated as source or installation declarations.

Both direct imports and indirect API use (`AbstractState`, `PropsSI`, `set_CoolProp_fluid`, and aliased property-module access) were checked. Dependencies' internal imports were not audited as thermopt code. This is a dated inventory: the test-suite work and unrelated working-tree changes may add or move files afterward.

After establishing the baseline suite, migrate active production accesses first, then active examples/development helpers, and finally obsolete test helpers/docstrings. Rerun component, cycle, optimization and plotting regression tests at each step. Keep the lockfile's transitive CoolProp dependency unless jaxprop itself changes its requirements.
