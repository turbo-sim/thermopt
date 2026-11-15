# Phase 3: native state access and setup requirements

Follow-up: the metadata fix is committed and verified through a separate package.
See [METADATA_FIX.md](METADATA_FIX.md). This original report retains the earlier
open-gate observations; the follow-up preserves the historical flags too.

Reference commit: `b211fa084b800e970d4bc86f594c02bfc7fef385` (Phase 2).
Both implementations use installed PyPI JaxProp 0.6.0 and the same locked
CPU/x64 environment. ThermOpt remains eager NumPy/SciPy.

**Implementation and numerical verification are complete. The strict historical
checkpoint performance gate remains open.** The controlled comparisons below
isolate the application change, while preserving three historical bicubic flags.
No clean Phase 3 checkpoint has been promoted and no budget has been widened.

## Implemented behavior

- Production components, cycle equations and export lookups use canonical
  property fields. Display labels and public input aliases remain separate.
- Nondimensional machinery now returns a native batched immutable state, like
  isentropic and polytropic machinery. Endpoints remain scalar, with original
  provenance, statuses, numerical values and solver behavior.
- Problem setup checks required fields/input pairs and component-state property
  references in constraint variables, targets and objectives. Unsupported
  requirements fail explicitly before reference preparation. Exchanger diagnostics
  that are supported but disabled remain explicit errors.
- Expression traversal accepts read-only mappings and reads guarded attributes
  once. Typed JaxProp property failures retain provider/index context and gain
  an expression note, instead of being wrapped as generic expression errors.
- Generic state serialization uses public lossless data, preserving enabled
  fields, identity, status and masks. Optional undefined values remain unavailable.
- Canonical plot axes preserve auxiliary streams. The blending example uses
  native construction/constants, public state stacking and an explicit smoothing
  import, with its original illustrative coefficients and physical assumptions.

See the [consumer inventory](../tests/STATE_ACCESS_060.md). No new fluid adapter,
numerical class, property-call wrapper, implicit clipping or source fallback is
introduced. The [repository guidance](../AGENTS.md) preserves the later JAX
direction and records static metadata equality as a performance concern.

## Verification and source fidelity

The final clean locked suite passes **396 tests**, with **four unchanged strict
expected failures**, zero failures and zero ordinary skips. It includes the
three standard optimization solves, existing JIT/gradient/batching checks and
77 new state-contract tests. These cover all three native model routes, singleton
and batch paths, all 14 case stream identities, all topology aliases, optional
property failures, source quality meaning, serialization, setup-only validation
and a complete 300-point headless blending example producing both figures.

Nine component/cycle modules have identical ASTs after the documented canonical
field spelling and nondimensional path-container transformations. This comparison
omits docstrings only in the audit representation; source explanations remain.
The math module and all 21 captured scientific-reference/input files are
byte-identical to the reference. No physical golden values or existing numerical
tolerances were regenerated. Synthetic numeric-mapping fixtures adopt canonical
keys without changing their numbers or expected equations.

The five-mode 200x200 release profile also passes. Working-fluid-only and
all-stream recuperated CO2 source/bicubic optimizations converge, pass feasibility
checks, agree at source optima, and pass independent source cross-evaluation.
Table generation is retained separately from warm calls; the new 200x200 CO2
table takes about 205 s to generate. Full-solve costs are current observations,
not a matched before/after optimization speed claim.

## Fresh-process timing observations

Three fresh before/after pairs use the unchanged phase harness, 12 synchronized
samples per backend and affinity `0xffff`, with balanced pair order. All six
cohorts have exactly identical physical snapshots, also matching Phase 2 hashes.

Against the stored Phase 2 three-cohort median, the current three-cohort median
still flags:

| Bicubic complete cycle | Ratio to stored checkpoint |
|---|---:|
| Recuperated CO2, all streams tabulated | 1.107 |
| Recompression CO2 | 1.136 |
| Nitrogen refrigeration | 1.228 |

Against the freshly rerun reference aggregate, all ten comparisons meet 1.10;
the largest increase is 5.1%. Individual fresh pairs have zero, zero and six
flags. Unchanged reference timing also shifts, so the old absolute flags alone
do not isolate application overhead. They remain retained and unresolved as a
strict historical checkpoint gate.

Separate complete source PTES measurements give +1.9% for the ordinary model
and +0.8% for nondimensional PTES in one matched fresh pair.

## Cache-identity investigation

Alternating two implementations in one process initially leaves a bicubic
polytropic flag of +13.1%. Profiling 20 complete cycles per variant finds the
same **4,580 bicubic requests**. The disadvantaged model additionally incurs
**13,480 dataclass equality calls** and **67,400 Mapping equality calls** in
installed JaxProp/JAX dispatch. The setup validator is outside these evaluations.

Two controls establish a dependency/cache-identity confounder:

- Two separately constructed models running **identical committed source** show
  a 28.8% apparent slowdown for the second variant.
- Reversing construction/initial warmup order makes current/reference 0.752 for
  the polytropic case. This is a cache-order observation, not an application
  speedup claim.

The cost arises when equivalent independent model metadata is compared during
compiled dispatch. It is pre-existing behavior in released JaxProp 0.6.0.
Phase 5 records the upstream investigation: efficient immutable metadata equality
or interning needs its own tests. Installed JaxProp is not patched, and mutable
evaluators are not globally shared as a production workaround.

## Controlled application comparison

The public functional cycle interface already accepts explicit native models.
Supplying the same model arguments to both implementations removes that identity
bias while retaining complete cycle output. Three bursts of 12 alternating calls
per version check full returned data, property availability, statuses and
identities before/after, and synchronize all outputs. Path-container normalization
for comparison uses public state stacking outside timings.

| Complete functional cycle | Source change | Bicubic change |
|---|---:|---:|
| Recuperated CO2, working-fluid selection | +0.4% | -5.5% |
| Recuperated CO2, all-stream selection | +1.6% | -2.0% |
| Simple CO2, polytropic machinery | -0.3% | -1.3% |
| Recompression CO2 | +0.2% | -3.3% |
| Nitrogen refrigeration | +0.1% | -3.1% |
| Source PTES | +0.3% | — |
| Source nondimensional PTES | +4.0% | — |

All twelve comparisons meet the unchanged 1.10 limit. These observations isolate
the changed numerical consumer code; they do not replace fresh-process problem
lifecycle measurements or erase the three historical flags. Full `fitness`
measurements, independent-model controls and setup costs remain retained.

## Evidence and next work

[PHASE3.json](PHASE3.json) retains raw samples, source/input/dependency identities,
test log/XML hashes, source-fidelity checks, all fresh-process flags, optimization
metrics, cache-order/identical-source controls and profiling details. Detailed
configs, logs, profiles and executable benchmark snapshots remain under ignored
`benchmark/results/phase3/`.

[state_paths.py](state_paths.py) and [compare_versions.py](compare_versions.py)
remain runnable benchmark helpers. The main phase runner is unchanged, so stored
checkpoint flags stay comparable rather than disappearing behind a new protocol.
The accepted checkpoint remains Phase 2.

Next performance work is to resolve the released-library metadata/cache issue
and reassess the strict historical gate with its controls. Additional component
physics/domain checks are Phase 4; broader lifecycle behavior is Phase 5.
Other platforms, licensed REFPROP and every standalone study remain unverified
by this phase's checks.
