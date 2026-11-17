# ThermOpt phase benchmarks

Maintained cycle programs now live under `case_studies/`. Active workload
builders and test entry points use those locations. Dated reports, checkpoints
and frozen comparison trees retain their original paths and hashes; they are
historical evidence, not rewritten results. The path relocation and updated
solver dependency change the current workload/environment identity, so earlier
timing checkpoints are not automatically comparable to a fresh run.

The [published 0.6 adoption record](RELEASE_060.md) describes the verified
installation and tests. The [first phase checkpoint](checkpoints/pypi_060_phase.json)
is available as a reference for a matching subsequent phase run.

The [Phase 2 assessment](PHASE2.md) retains three matched controlled cohorts,
including all individual flags and the original unrestricted observations.
Its [new checkpoint](checkpoints/phase2.json) records CPU affinity explicitly;
use its matching protocol for subsequent phase comparisons.

The [Phase 3 assessment](PHASE3.md) records completed state-interface verification,
matched complete-consumer measurements and a confirmed model-metadata cache
confounder. Its strict historical timing gate remains open; Phase 2 remains the
accepted checkpoint.

The [metadata-fix follow-up](METADATA_FIX.md) verifies the local JaxProp fix and
records an [artifact-specific checkpoint](checkpoints/metadata_fix.json) for the
next phase. Its historical flags and targeted follow-up remain visible; use the
recorded environment/artifact when comparing it.

The [Phase 4 assessment](PHASE4.md) verifies component physics across the three
native models, complete-return polytropic speedups, the full test suite and the
200x200 release profile. Its [checkpoint](checkpoints/phase4.json) uses the same
local JaxProp artifact. First-use compilation observations and older absolute
timing flags remain documented follow-ups.

The [Phase 5 assessment](PHASE5.md) verifies problem-owned model/reference reuse,
atomic reloads and first-use host-solver compilation. Its
[checkpoint](checkpoints/phase5.json) starts benchmark protocol 2 with the same
unpublished local JaxProp artifact. Matched complete-cycle and balanced component
gates pass; earlier raw flags and one-time setup/compilation costs remain recorded.

This is the canonical benchmark entry point for the
[JaxProp 0.6 migration plan](../z_JAXPROP_060_IMPLEMENTATION_PLAN.md).

The focused [Phase 6 reporting check](PHASE6.md) separates complete cycle
evaluation from YAML/Excel/plot/KKT-report costs:

```powershell
python benchmark/output_workflows.py --output benchmark/results/output_check
```

It compares the normalized-input exception boundary with its inherited
PySolverView route on the same configured native model. It is not a new release
benchmark or a substitute for the full phase/release profiles below.

ThermOpt stays eager NumPy/SciPy. The maintained workload implementation is
`benchmark/native_backends.py`. The historical development command delegates
to it; dated reports retain their original source identities. Current solver
comparisons live in `recompression_solvers.py` and `recuperated_slsqp.py`, with
thin case-local launchers preserving their documented commands.

The [Phase 7/8 assessment](PHASE7_8.md) records study-owned models, five additional
complete case runs, local Python 3.11/3.13 verification and installed-wheel QA.
Hosted CI and older absolute timing gates remain explicitly open.

Use the environment containing the released JaxProp 0.6 package, with no
checkout or `PYTHONPATH` substitution:

```powershell
python -m pytest
python benchmark/run.py --profile smoke
python benchmark/run.py --profile phase
python benchmark/run.py --profile release
```

- `smoke`: recuperated CO2, one 100x100 working-fluid table, no full solve.
- `phase`: five source/native-table case or stream configurations at 100x100,
  including simple/polytropic CO2, recompression and nitrogen refrigeration.
- `release`: workloads at 200x200 plus full source/bicubic recuperated optimizations,
  source cross-evaluation, feasibility and objective agreement.

`state_paths.py` complements these profiles with complete source PTES consumers,
including nondimensional radial machinery and its returned state paths:

```powershell
python benchmark/state_paths.py --case ptes --cpu-affinity 0xffff --output benchmark/results/ptes_reference.json
python benchmark/state_paths.py --case ptes_turbo --cpu-affinity 0xffff --output benchmark/results/ptes_turbo_reference.json
```

Use new output paths for each run. It retains 12 synchronized full-cycle samples,
setup cost, complete numerical snapshots, CPU affinity and dependency/source/input
identities. Compare matching before/after records; ordinary tests have no timing
assertions. Like the main runner, this helper currently shares frozen workload
definitions with the regression suite; Phase 7 will consolidate them.

For small interface changes, alternate both source implementations in one process:

```powershell
python benchmark/compare_versions.py --reference-root PATH_TO_REFERENCE_TREE --output benchmark/results/paired_run --cpu-affinity 0xffff
```

This loads real source packages in separate namespaces and retains complete-cycle
outputs, raw samples and an executable harness snapshot. `--workload` restricts
cases, `--first-version current` counterbalances construction/warmup order, and
`--identical-source-control` uses reference code for both independent instances.
`--shared-models` compares the existing functional cycle API with identical native
model arguments to isolate application changes from metadata-identity cache costs.
That mode excludes the problem object's bookkeeping; keep full `run.py` and
`state_paths.py` measurements as lifecycle/end-to-end evidence too. No fluid
backend or numerical implementation is patched by these helpers.

For full component returns across CoolProp, perfect-gas and bicubic models:

```powershell
python benchmark/compare_components.py --reference-root PATH_TO_REFERENCE_TREE --cpu-affinity 0xffff --output benchmark/results/components_before_first
python benchmark/compare_components.py --reference-root PATH_TO_REFERENCE_TREE --first-version current --cpu-affinity 0xffff --output benchmark/results/components_after_first
```

This covers both machines with isentropic, polytropic and nondimensional
efficiencies, plus exchanger and valve controls: 24 workloads. Each version
owns independent reusable native models. It compares every returned property,
optional NaN, validity mask, status and identity outside timing, retaining raw
samples for three balanced cohorts. Output data use the public state serializer;
backend models are never replaced. New output directories are required.

Model setup and the first component invocation are recorded separately. That
first invocation follows inlet evaluation and earlier workloads in the same
process; it is **not** an isolated cold-JIT measurement. Reverse construction
order and retain both observations. New batch shapes can incur one-time JIT
compilation even when sustained evaluation is faster.

First-time property-table generation can take several minutes. The shared
`benchmark/results/tables/` cache is reused only after JaxProp validation.
Construction/table setup, warmed evaluations and complete solves remain
separate measurements. Never compare fresh generation against cached loading.

`lifecycle.py` measures complete configuration reloads and cycle evaluations
across source, bicubic and an independently specified perfect-gas Air cycle:

```powershell
python benchmark/lifecycle.py --reference-root PATH_TO_REFERENCE_TREE --cpu-affinity 0xffff --output benchmark/results/lifecycle
```

It includes every source-provenance/rendering operation in reload setup time,
records first construction/evaluation separately, and checks complete native
state data through the public serializer. Unchanged and nonfluid configuration
updates retain the current problem's stream-model identities. It does not request
quality from models that do not support it.

`first_use.py` separates first-use compilation from sustained component timings:

```powershell
python benchmark/first_use.py --reference-root PATH_TO_REFERENCE_TREE --cpu-affinity 0xffff --output benchmark/results/first_use
```

Each version/model/machine/warmup/cohort runs in its own process, using installed
JaxProp, CPU/x64, hash seed zero and disabled persistent compilation caching.
Import, construction, inlet evaluation, explicit warmup, the first complete
component and twelve warm returns are recorded separately. Default paths contain
10 compressor or 50 expander nodes. All returned fields, masks and identities
are compared; no property projection or replacement backend is used.

- `none`: only the inlet PT evaluation precedes the component.
- `scalar`: explicitly prepare HP/PS scalar representations used by both source
  versions (Python and NumPy scalar types); no reporting batch is warmed.
- `path`: additionally prepare a complete HP batch at the reporting shape.

This is protocol 2. Earlier protocol-1 investigations warmed fewer scalar
representations; their scalar/path results are retained as diagnostics and are
not interchangeable checkpoints. `--log-compiles` records a separate diagnostic
run with compiler logging overhead. Explicit warmup moves compilation into setup;
it does not eliminate that work or happen automatically in the fluid constructor.

When assessing independent-model component comparisons, use both construction
orders with equal weight. Run each order twice in fresh processes if a flag
changes direction with the order. Retain individual flags and use the median of
those four per-run medians for a balanced aggregate; the 1.10 budget is unchanged.
An unchanged valve/exchanger control helps detect native JIT metadata/cache
identity effects. Do not treat an order-dependent control as a measured code
speedup, or share mutable models globally to mask its cost.

Each run writes a new directory containing `manifest.json`, full workload JSON,
effective YAML inputs and logs. Raw synchronized timing samples and scientific
comparison results are in each workload's summary. Existing runs are never
overwritten. All generated data under `results/` is ignored by Git.

To compare the next phase against an accepted checkpoint, use the same profile,
grid, sample count, machine and numerical environment:

```powershell
python benchmark/run.py --profile phase --baseline benchmark/results/ACCEPTED_RUN/manifest.json
```

The runner rejects incompatible protocols, effective physical/solver inputs or
environments. Main-runner protocol 3 retains protocol 2's `PYTHONHASHSEED=0` for
each workload subprocess and adds the consolidated native workload source.
It removes the obsolete `exceptions` normalization. Protocol-1/2 checkpoints
are incompatible; retain them as historical evidence.
Both source versions must use the same harness when capturing a new matched
baseline. This improves reproducibility without proving that hash randomization
caused any prior timing flag.

CPU affinity is observed and included in that signature. On Windows
or Linux, an optional mask constrains only the benchmark process and its children,
before numerical libraries create worker threads:

```powershell
python benchmark/run.py --profile phase --cpu-affinity 0xffff
```

Choose a nonempty subset of the processors available to your process. Both sides
of a comparison must use the same mask and protocol; do not compare a newly pinned
run with an old unrestricted checkpoint. Affinity does not reserve those CPUs or
eliminate shared-machine load. It changes no machine-wide power settings.

The runner excludes the table storage directory from the physical input
signature; source and cache identities are otherwise validated by JaxProp.
Numerical source hashes may change between phases and are always recorded.

The Phase 2 records used a historical workload builder that removed the old
`exceptions` field. That normalization is absent from the current builder and
ThermOpt rejects the obsolete field. The Phase 7/8 release run establishes a
protocol-3 baseline candidate; it does not establish a speedup against an
incompatible protocol-2 artifact, solver or dependency environment.

The provisional warm budget is 1.10 times the preceding accepted median for
each matching model/workload. Exit code 2 means timing flags require assessment,
not that a noisy sample proves an intrinsic regression. Repeat flagged cases
in at least three fresh cohorts and inspect unchanged controls, raw samples and
host conditions. Do not discard slow samples, relax limits silently or replace
full-state consumers with selected-property timings. Setup/compilation/full-solve
costs need their separate recorded assessment; this first runner does not yet
automate their between-phase gates.

Ordinary correctness tests contain no timing assertions. Hosted shared-runner
CI timing is advisory; controlled scheduled/manual runs carry performance
acceptance. The initial smoke checkpoint validates the infrastructure and the
published installation, not full phase/release performance acceptance.

## Study-model ownership and CI

For a matched full-evaluation comparison of the Rankine study evaluator, prepare
an independent reference checkout containing the study and ThermOpt sources:

```powershell
python benchmark/study_models.py --reference-root PATH_TO_REFERENCE_TREE --samples 60 --output benchmark/results/study_models_run
```

Four fresh processes run reference/current/current/reference. Every complete
residual vector must agree before medians are compared. Raw samples, source
identities, dependency versions, setup costs and the 1.10 flag are retained.
This measures the actual study `fitness`, including native full-state calls;
it is separate from the general cycle release profile.

Ordinary CI runs the locked tests and a bounded numerical smoke profile on
Windows/Ubuntu and Python 3.11/3.13. The manual benchmark workflow records the
selected profile and uploads artifacts. Shared-runner stopwatch observations
do not replace matched controlled cohorts or automatically close timing flags.
