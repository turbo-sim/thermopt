# Using the JaxProp 0.6 interface in ThermOpt

This checkout requires JaxProp 0.6 and deliberately rejects the old ThermOpt
fluid schema. Install the locked dependencies with `poetry install`. ThermOpt's
cycle and component calculations remain eager NumPy/SciPy. A future JAX port is
separate work; configuration/provenance checks belong at setup boundaries.

## Fluid configuration

Use explicit native construction settings:

```yaml
working_fluid:
  name: CO2
  method: coolprop
  backend: HEOS
```

For native tables, use `method: bicubic` and supply the domain/resolution in
`table_options`. See the validated inputs in the [case-study guide](case_studies/README.md).
ThermOpt constructs actual JaxProp evaluators; all state calculations use their
public state interface. There are no fluid adapters or implicit source fallbacks.
Perfect gas and tables must support the properties/input pairs a case requires;
phase-change initialization and plotting use explicit source capabilities.

Constructor-level `exceptions` and obsolete fluid fields are rejected. Solver
trials use the shared normalized-input error boundary so controlled SLSQP/Pounce
can reject recognized invalid thermodynamic states. Direct physical evaluation
retains the native JaxProp exception. Unsupported capabilities remain explicit.
The error boundary never clips a failed optimizer trial into a valid state.

Stream models and reference data belong to the problem. Unchanged configuration
reloads reuse models; a real model/source-reference change rebuilds them at setup.
The former global `get_property_fluid` helper has been removed. Study utilities
must receive or own their native models rather than sharing a mutable global one.

## Solvers, outputs and saved studies

Standard SLSQP YAMLs now select PySolverView's controlled initial step. Physical
inputs and per-case acceptance tolerances remain explicit. Deliberate SciPy,
IPOPT and Pounce comparisons are separate programs under `benchmark/` with
case-local launchers. Ordinary demos are listed separately from research sweeps.

Generated case artifacts belong under each case's `results/`. Older research
`output/` trees remain historical inputs, with no rewritten signatures. Rankine
study caches now include core/study code, dependencies and source-reference
identity. Changed provenance invalidates reuse; choose a fresh output root after
migration. This work does not claim that full scientific grids were regenerated.

Manager snapshots, iteration YAML and `optimal_solution.yaml` contain complete
`solver_options` and `problem_formulation` and can initialize a manager directly.
The problem-level writer emits the formulation for `ThermodynamicCycleProblem`.
Both retain variable order, floating-point precision and construction recipes.
Study readers accept these two documented payloads without translating old fluid
settings. Excel cells contain native numbers; unavailable properties remain blank.

## Verification status

The [Phase 7/8 assessment](benchmark/PHASE7_8.md) records local Python 3.11/3.13
tests, bounded case runs, a complete release numerical profile, installed-wheel
checks and matched study timings. CI is prepared for Windows/Ubuntu with both
Python versions; its hosted runs have not yet been executed. Ordinary CI has
numerical gates; manual benchmark artifacts support controlled timing review.

The official locked JaxProp artifact remains 0.6.0. The separately verified local
metadata fix `52c98b8` is still an upstream release follow-up. Old absolute timing
flags remain open because their artifact/solver/protocol context differs. A
protocol-3 baseline candidate does not close those comparisons. Direct JaxProp
NumPy-valued static-option warnings have a [local upstream repro](benchmark/JAXPROP_STATIC_OPTIONS.md).
ThermOpt normalizes these options once during setup.
