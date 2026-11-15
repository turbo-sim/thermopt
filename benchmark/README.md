# ThermOpt phase benchmarks

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

This is the canonical benchmark entry point for the
[JaxProp 0.6 migration plan](../z_JAXPROP_060_IMPLEMENTATION_PLAN.md).
ThermOpt stays eager NumPy/SciPy. Current workload code is reused from
`development/bicubic/validate_060.py`; consolidating shared implementations here
is a later checklist item, not a second copy of the thermodynamic equations.

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
environments. CPU affinity is observed and included in that signature. On Windows
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

For the Phase 2 comparison, the unchanged historical workload builder removes
the old `exceptions` field when constructing its comparison inputs. This keeps
the pre/post scientific workload identical; it is not the production
configuration path. ThermOpt's actual constructor rejects that field. Removing
this obsolete benchmark normalization is part of Phase 7's harness consolidation
and requires an explicitly checked new protocol/checkpoint.

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
