# ThermOpt phase benchmarks

The [published 0.6 adoption record](RELEASE_060.md) describes the verified
installation and tests. The [first phase checkpoint](checkpoints/pypi_060_phase.json)
is available as a reference for a matching subsequent phase run.

The [Phase 2 assessment](PHASE2.md) retains three matched controlled cohorts,
including all individual flags and the original unrestricted observations.
Its [new checkpoint](checkpoints/phase2.json) records CPU affinity explicitly;
use its matching protocol for subsequent phase comparisons.

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
