# Phase 2: direct JaxProp fluid configuration

The reference is commit `8f56ccc92d712a17f8216e70b781ee4bb04c6e22`, using
the published JaxProp 0.6.0 package. This is a ThermOpt phase comparison; both
sides use the same released JaxProp, locked dependencies and eager CPU/x64
execution. It is not a comparison with JaxProp 0.5.

## Configuration changes

- All 36 maintained cycle YAML files now explicitly select `method: coolprop`:
  108 stream declarations across examples, templates, assignments, documentation
  and regression inputs. Only the obsolete `exceptions` field is replaced.
- `fluid_from_configuration` delegates to `jaxprop.Fluid` without translating
  old fields. JaxProp's typed errors explain obsolete/unknown constructor
  options. Error policies belong on state calls; normal solver calls raise.
- Static NumPy scalars and read-only nested mappings are normalized once when
  constructing an evaluator. Configuration is not mutated. Source, perfect-gas
  and bicubic models remain actual JaxProp objects.
- Historical saved solutions remain data. The current wet/dry study is tested
  with both independently generated and saved design-variable seeds; executable
  archives are not a compatibility route for old constructor settings.

The numerical equations, property-call paths, solver settings, scientific
references and tolerances are unchanged. Existing tests explicitly verify that
configured cycle evaluation does not reconstruct fluid models. The complete
document signatures in `tests/data/fluid_configuration_phase2.json` protect all
physical, solver, constraint and output settings during this schema migration.
Future JAX design guidance is retained in [AGENTS.md](../AGENTS.md); this phase
adds no mutable numerical object or cycle JAX execution.

## Verification

The configuration tests cover native evaluator classes, identifiers, unchanged
inputs, read-only mappings, scalar normalization, obsolete/unknown fields and
state-call error policy. All 36 maintained YAML cases have schema and physical
signature checks. Eight additional deterministic benchmark-affinity checks
validate masks and observation without stopwatch assertions.

The final clean locked-environment run passes **319 tests**, with **four existing
strict expected failures**, zero failures and zero ordinary skips. All three
complete optimization solves remain part of that run. No numerical golden files
were regenerated. The report retains the full test log and JUnit XML hashes.

## Complete-cycle timing assessment

The first unrestricted before/after cohort had six flags, with increases of
15–41% in the affected workloads and improvements in others. Those observations
are retained in [PHASE2.json](PHASE2.json).

The repeat experiment uses the same new affinity-aware runner on both frozen
reference and current source trees, leaving workload code and frozen physical
inputs unchanged. This protocol change is explicit: the original runner and
original measurements remain retained. Both processes use mask `0xffff`, inherited
by subsequently spawned workers. This constrains allowed CPUs, without reserving
them or changing machine-wide power settings. Affinity is part of the comparison
signature and is checked again after the run.

Three fresh matched pairs use balanced order: pre/post, post/pre, pre/post. Each
cohort measures 12 alternating, synchronized complete cycle evaluations for each
backend, after warmup. The result below is the ratio of the three post-change
cohort medians' median to the corresponding pre-change median. All samples and
per-pair flags are retained; the original 1.10 limit is unchanged.

| Complete cycle / stream selection | Source change | Bicubic change |
|---|---:|---:|
| Recuperated CO2, working fluid tabulated | +2.2% | +0.7% |
| Recuperated CO2, all streams tabulated | -1.3% | -1.8% |
| Simple CO2, polytropic machinery | +7.1% | +5.4% |
| Recompression CO2 | +5.7% | +7.1% |
| Nitrogen refrigeration, all streams tabulated | +0.1% | -3.7% |

All ten aggregate comparisons meet the warm budget. Individual controlled pairs
still have five, zero and two flags respectively. In particular, both polytropic
variants exceed 10% in two individual pairs; the aggregate changes are 7.1% and
5.4%. This is evidence of acceptable aggregate cost on this profile, not proof
that every run or every machine has identical performance. Unchanged reference
code also varies between cohorts. Affinity alone has not been established as
the cause of the original shifts, and no intrinsic speedup is claimed.

All eight cohorts (the original two and the six repeats) have exactly identical
source and bicubic physical snapshots before/after. Every workload also passes
the existing source-versus-table interpolation budgets. Cached setup costs are
retained separately and are not substituted for warm evaluation measurements.
This setup-only phase does not change component/cycle numerics, so the complete
optimization benchmark profile was not repeated; full optimization correctness
is covered by the full suite.

## Retained checkpoint and rerunning

[PHASE2.json](PHASE2.json) retains raw samples, individual medians and flags,
setup timings, numerical agreement, snapshot hashes, source/protocol/dependency
identities and the aggregation rule. Detailed generated configs/logs and full
test XML remain in ignored `benchmark/results/phase2/`.

The [Phase 2 checkpoint](checkpoints/phase2.json) contains the current three-cohort
median timings and their constituent observations. On this matching machine and
locked environment, the next phase can use:

```powershell
python benchmark/run.py --profile phase --samples 12 --cpu-affinity 0xffff --baseline benchmark/checkpoints/phase2.json
```

Other hardware, affinity, versions or protocol hashes require a new matched
reference. The earlier unrestricted checkpoint remains historical evidence and
is incompatible with this affinity-aware protocol. Windows affinity was
exercised here; the Linux affinity branch and other platforms were not executed.
