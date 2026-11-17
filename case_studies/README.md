# ThermOpt case studies

This directory contains the cycle cases and research programs formerly under
`examples/`. Each case retains its original subdirectory layout, YAML inputs,
scripts, reference material and existing local results. The existing
`partial_evaporation_orc_offdesign_prompt.md` is also preserved here.

Available case families include:

- `sCO2_recuperated/`, `sCO2_recompression/`, `sCO2_comparison/`.
- `rankine_steam/`, `rankine_organic/`, `rankine_wet2dry/`,
  `rankine_transcritical_vs_trilateral/`.
- `heat_pump_CO2/`, `refrigeration_air/`, `refrigeration_R134a/`.
- `PTES_CO2/`, `demo_hx_node_snap/`, `HX_development/`.

Read the individual case instructions before running it. Scientific sweeps,
historical figure preparation and input-dependent programs have different
requirements; relocating them does not establish that every full study runs.
Keep generated outputs, tables and logs inside the case's `results/` directory.

The [small-demo execution report](SMALL_DEMOS.md) lists the one/two-case
launchers verified on 16 November 2026, their fresh outputs, numerical fixes
and headless invocation. Scientific sweeps and benchmark suites have a separate
verification scope.

The [execution inventory](EXECUTION_INVENTORY.json) classifies all 47 Python
programs, with commands, working directories, required inputs, budgets and
acceptance criteria. It distinguishes ordinary demonstrations from full
scientific sweeps, saved-input tools and historical figure preparation. The
[Phase 7/8 assessment](../benchmark/PHASE7_8.md) records the additional bounded
runs; full research grids have not been rerun as part of this migration.

## Selecting a native evaluator

The following launchers accept `THERMOPT_CONFIG`. Run from the case directory;
the variable selects a YAML file without editing the default input.

| Directory | Source EOS input | Native bicubic input |
| --- | --- | --- |
| `sCO2_recuperated` | `case_sCO2_recuperated.yaml` | `case_sCO2_recuperated_bicubic.yaml` |
| `sCO2_comparison` | `case_sCO2_simple.yaml` | `case_sCO2_simple_bicubic.yaml` |
| `sCO2_recompression` | `case_sCO2_recompression_coolprop.yaml` | `case_sCO2_recompression.yaml` (default) |
| `refrigeration_air` | `case_air_refrigeration_recuperated.yaml` | `case_air_refrigeration_recuperated_bicubic.yaml` |

For example, from `case_studies/sCO2_recuperated` in an activated environment:

```powershell
$env:THERMOPT_CONFIG = "case_sCO2_recuperated_bicubic.yaml"
$env:DISABLE_PLOTS = "1"
python run_optimization.py
Remove-Item Env:THERMOPT_CONFIG
```

These four bicubic inputs use 200x200 working-fluid tables and enable JAX x64 in
the launchers. Their domains cover the validated single-phase operating cases;
they are not a universal table for arbitrary initial guesses or other fluids.
The external streams retain their configured source EOS. First construction can
generate tables; later runs reuse JaxProp-validated files in `results/property_tables`.
Changing an evaluator must retain suitable table bounds and supported input pairs.
Phase-change Rankine/HX/PTES studies continue to use their source models.

All native evaluators are ordinary JaxProp objects with `get_state(...)` calls.
There is no application-side bicubic adapter or fallback to a different EOS.

Standard SLSQP cases and study templates use PySolverView's controlled backend:

```yaml
solver_options:
  library: pysolver
  method: slsqp
  update_on: iteration
  extra_options:
    first_step_fraction: 0.2
    max_step_fraction: 1.0
    step_growth_factor: 1.5
    feasibility_tol: 1e-6
    optimality_tol: 1e-4
```

Each case retains its physical inputs, scaling and iteration budget. Feasibility
tolerances are explicit (some cases use `1e-8` or `1e-12` rather than the
illustrative `1e-6`). Native stopping can be stricter: the PTES turbomachinery
case uses `1e-8` native stopping with its original `1e-6` feasibility gate.
R134a uses three-point derivatives to avoid false small-step termination near
saturation; the other standard cases keep their configured derivative method.
The optimality threshold is `1e-4` for finite-difference KKT checks. Accepted
iterations drive reporting. Study runners inherit these settings from their
templates. The named IPOPT configuration and deliberate SciPy/Pounce solver
comparisons remain explicit alternatives, rather than default SLSQP runs.

Study cache signatures include solver settings. Preserve results produced with
older settings as historical runs and choose a fresh results location when
rerunning those sweeps; do not edit saved signatures to reuse them as new solves.

The default-solver migration passed 564 tests with four unchanged strict expected
failures, including complete controlled solves of recuperated sCO2, R134a
refrigeration and a butane ORC. A separate two-iteration check exercised all six
intentional SciPy/controlled comparison configurations; those bounded runs are
not claimed as converged optimizations or full research sweeps.

For the verified six-way recompression solver/backend comparison, run from the
repository root after installing the locked dependencies:

```powershell
.\.venv\Scripts\python.exe case_studies/sCO2_recompression/compare_solvers.py --no-show
```

This compares controlled SLSQP, Pounce interior point and Pounce SQP using
CoolProp HEOS and a runner-defined 128x128 native bicubic table. It preserves
the input YAML and saves figures, histories and an ASCII timing comparison in
`sCO2_recompression/results/solver_comparison/`.

The implementation lives in `benchmark/recompression_solvers.py`; the case-local
command remains a thin launcher. Likewise, `sCO2_recuperated/compare_slsqp.py`
launches `benchmark/recuperated_slsqp.py`. These comprehensive comparisons are
separate from the ordinary demonstration profile.

Most `run_optimization.py` launchers expect their case directory as the working
directory. Standalone comparison runners resolve their own input paths; their
case README documents the supported invocation.

Manager-saved `optimal_solution.yaml` and iteration YAML files contain a complete
configuration, including solver settings and the evaluated design, and can be
passed directly to `ThermodynamicCycleOptimization`. Calling
`problem.save_current_configuration(...)` instead saves a formulation for
`ThermodynamicCycleProblem`. Both preserve native fluid recipes and variable
order. See the [Phase 6 output assessment](../benchmark/PHASE6.md).

Tests keep their frozen physical inputs under `tests/data/cases/`. Historical
benchmark JSON, checkpoints and audits retain the original `examples/` paths
and numerical hashes. Those records describe the run that produced them; active
source code and current instructions use `case_studies/`.

The relocation was verified on 16 November 2026: all 660 tracked files were
recognized as renames, all original local files were retained, and physical
YAML/data/figure bytes were unchanged. The full suite passed 543 tests with four
unchanged strict expected failures. The maintained recuperated source/bicubic
smoke profile also passed its numerical checks from the new layout; this is
path/workflow verification, not a new between-phase performance acceptance.
