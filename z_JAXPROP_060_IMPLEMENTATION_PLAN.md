# ThermOpt migration to the JaxProp 0.6 interface

## Objective and scope

Use the released JaxProp 0.6 interface consistently throughout ThermOpt and its
maintained examples. Backward compatibility with JaxProp 0.5 and old ThermOpt
fluid settings is not a requirement. Keep cycle/component equations in eager
NumPy/SciPy: this plan does not port ThermOpt to JAX. Existing isolated native
JAX machinery tests remain regression coverage, not a new execution strategy.

Each phase must leave a reviewable implementation, tests, measured application
costs and a checked status. Preserve scientific explanations and useful
docstrings. Do setup/configuration checks once; do not add repeated provider,
reference or metadata validation to property calculations.

The user's future JAX direction is persistent repository guidance in
[AGENTS.md](AGENTS.md). Prefer explicit functional numerical interfaces and
predictable data structures, keeping mutable orchestration/I/O separate.
Assess dataclasses deliberately; frozen records alone do not establish JAX
transformability. This is preparation for a later sweep, not permission to
port the current cycle implementation now.

## Starting evidence

- [x] Prepared-wheel compatibility checkpoint: 253 tests passed, zero failures
  or skips, four explicitly tracked expected failures. Coverage includes 14
  cycle configurations at initial/nearby physical designs, three complete
  optimizations, plots/exports and five native-bicubic integration tests.
- [x] Native bicubic comparison: recuperated CO2 with 100/200 grids, all three
  streams tabulated, simple/polytropic CO2, recompression and nitrogen
  refrigeration. See [verification](development/bicubic/VERIFICATION_060.md).
- [x] Preserve original physical operating-point fixtures and their tolerances.
  Separate 0.6 records describe corrected initialization, 36 generalized-quality
  annotations and two R134a source-diagnostic residuals. Never regenerate
  expected outputs merely because a refactor makes a test fail.
- [x] Live official PyPI metadata confirms JaxProp **0.6.0**, uploaded on
  **15 November 2026 at 13:14 UTC**. A git tag or local version constant alone
  is not evidence of publication.
- [x] Close the released-artifact/dependency gate below. Prepared-wheel evidence
  remains dated evidence; do not relabel it as a published-wheel test run.

The incoming worktree contains substantial unrelated uncommitted changes.
Preserve them and keep phase changes separately reviewable; no automatic bulk
commit or reference rewrite is part of this plan.

## Per-phase acceptance protocol

1. Record the starting source/lock/dependency identities and the specific
   behavior being changed. Add missing meaningful regression cases first.
2. Implement one coherent phase. Use actual JaxProp evaluator objects and
   ordinary `get_state(pair, first, second, ...)` calls; no fluid adapters,
   replacement thermodynamic equations or hidden backend fallback.
3. Run targeted tests, then `python -m pytest` before closing the phase. Preserve
   strict expected failures with their reasons; an unexpected pass needs review.
4. Run `python benchmark/run.py --profile phase`, compare against the preceding
   accepted benchmark manifest, and run `--profile release` for component/cycle
   numerical changes and release boundaries. Separate model/table setup,
   compilation/warmup, complete warm evaluation and full optimization costs.
5. Investigate timing flags before acceptance. Initial warm-performance budget:
   at most 10% slower median per matching workload/model than the preceding
   accepted checkpoint. A noisy one-cohort flag requires repeats; use at least
   three fresh cohorts before calling a material change repeatable. Do not
   silently widen budgets, discard slow samples or substitute property-only
   kernels for complete cycle evaluations. Document any explicitly accepted cost.
6. Update this checklist and link exact test/benchmark evidence, limitations and
   newly discovered issues. Commit a completed phase only when requested, with
   a descriptive message explaining behavior and verification.

Tests and timings serve different purposes. Correctness tests must not contain
fragile wall-clock assertions. A timing run must first pass its numerical
checks. Compare like-for-like hardware, Python, numerical dependencies, solver
options, precision, table domain/resolution and physical inputs. If those
conditions change, declare the comparison incompatible or design a separate
matched experiment; never label it a speedup automatically.

## Phase 0 — Establish maintained verification and benchmark entry points

**Status: foundation implemented and verified on 15 November 2026.** Canonical
directories are `tests/` and `benchmark/` (singular). The initial checkpoint
records one cohort; it does not claim between-phase performance acceptance.

- [x] Keep existing source/component/cycle/optimization/plot/export tests and
  their independently verified reference corrections.
- [x] Add `benchmark/run.py` with smoke, phase and release profiles, immutable
  run manifests, raw samples, complete-cycle synchronization and optional
  comparison against a previous accepted manifest.
- [x] Reuse the already validated native workloads during this transition,
  without duplicating equations. Moving shared definitions out of the
  version-specific development script and test helpers remains a Phase 7 item.
- [x] Capture a first benchmark on the published-wheel environment and freeze
  it as the starting checkpoint. Match source/lock identities in every run.
- [x] Add nine deterministic checkpoint-protocol tests: known ratios/budget
  boundary, incompatible or numerically failed references, missing workloads,
  and invalid/nonfinite measurements. No stopwatch assertions enter pytest.
- [x] Assign remaining coverage gaps to Phases 2–6: configuration changes,
  invalid options, unsupported diagnostics, failed/out-of-domain states and
  model-specific saturation operations. Those later tests remain unimplemented.

Acceptance: documented commands actually run, preserve logs/results and fail
clearly on numerical errors, incompatible references or unresolved timing flags.
The historical 0.5 adapter stays outside the maintained execution path.

## Phase 1 — Adopt the published dependency and reproducible installation

**Status: complete locally on 15 November 2026.** Do not maintain a parallel
JaxProp 0.5 runtime route. See [release adoption evidence](benchmark/RELEASE_060.md).

- [x] Declare `jaxprop = "^0.6.0"` in `pyproject.toml`.
- [x] Download/hash-check the actual PyPI wheel; compare its Python source
  against the prepared/tested artifact without treating differing ZIP hashes
  or line endings as numerical differences.
- [x] Refresh `poetry.lock` using the published package, preserving already
  locked compatible dependencies where possible. Explain required changes.
- [x] Update the normal ThermOpt environment and editable package metadata;
  prove imports come from installed `site-packages` with no `PYTHONPATH` shadow.
- [x] Run dependency checks, the complete suite and a benchmark phase profile
  against that installation. Record versions and artifact hashes.
- [x] Verify the lock with CI's Poetry 2.2.1 and perform a clean locked
  main-dependency installation, separate from the development environment.
  Hosted CI and other Python/platform profiles remain Phase 8 work and are
  not claimed as executed here.
- [x] Update active installation/example instructions to the released package;
  preserve earlier preparation reports as historical evidence.

Acceptance: normal installed-package imports work and the dependency/lock check
passes. A working `PYTHONPATH` override alone does not close this phase.

## Phase 2 — Make fluid configuration use the new interface directly

**Status: planned.** The temporary `exceptions: True` translation is explicitly
scheduled for removal, not retained as the desired public interface.

- [ ] Migrate maintained example/case YAML and corresponding test input schemas
  to explicit `name`, `method`, `backend`, `identifier`, `table_options` and
  applicable `model_options`. Remove constructor-level `exceptions` settings.
- [ ] Remove the legacy translation in `fluid_from_configuration`; reject
  obsolete/unknown fields with useful guidance. Do not accept old positional
  backend spelling or silently reinterpret constructor policies.
- [ ] Keep error policy on state calls. Normal ThermOpt solver calls raise on
  invalid states; mask/clip are explicit choices only where the algorithm has
  defined semantics and tests. Never hide a failed optimizer trial with clipping.
- [ ] Normalize static YAML scalar settings once and preserve native model
  objects. No state-call wrappers or backend-specific fluid subclasses.
- [ ] Cover source, bicubic and perfect-gas construction, malformed settings,
  input immutability and identifiers. Update tests to the new schema without
  changing their physical inputs or expected numerical outputs.

Acceptance: maintained cases use the new vocabulary, no old-interface shim is
needed, complete source tests pass and the phase benchmark meets its budget.

## Phase 3 — Standardize state access and capability requirements

**Status: planned; basic capability-aware diagnostics are already implemented.**

- [ ] Inventory core versus optional fields used by each component, constraint,
  study, plot and exporter. Use canonical property names in maintained code;
  keep friendly display labels separate from field access.
- [ ] Resolve required capabilities during component/problem setup. Request
  optional source diagnostics only when supported and needed; a constraint
  requiring an unsupported property must fail clearly, never disappear.
- [ ] Use the public immutable-state/index/stack/serialization interfaces.
  Remove dataclass/private-field mutation and obsolete `states_to_dict` helpers.
- [ ] Separate unsupported, disabled and locally unavailable properties. Preserve
  statuses, optional NaNs and honest blank report cells; never invent zeros.
- [ ] Test scalar/singleton/batch shapes, immutability, source/table identities,
  optional-property failures, source diagnostics and source-model-specific
  quality meaning outside the two-phase region.

Acceptance: same complete state/property behavior through supported routes,
typed informative failures elsewhere, and no new per-property validation work.

## Phase 4 — Verify eager component physics across applicable models

**Status: planned; source regressions and initial bicubic checks exist.**

- [ ] Test compressor/expander isentropic and polytropic paths, heat exchangers,
  valves and functional nondimensional machinery using actual supplied models.
- [ ] Preserve SciPy ODE/root solving, scalar states, energy identities,
  endpoint recovery, finite physical outputs and input/model reuse. Keep the
  numerical solver tolerances and acceptance checks unless a defect is proven.
- [ ] Expand bicubic component coverage at realistic table boundaries and
  applicable HP/PT/PS/HS routes. Require explicit failure for unavailable paths.
- [ ] Add perfect-gas component and applicable eager cycle checks against its
  own ideal-gas identities; do not require it to reproduce real-fluid physics.
- [ ] Verify saturation-node placement only with a saturation-capable model.
  Test real two-phase source cases and explicit rejection for unsupported
  PG/table configurations. Do not infer universal phase interchangeability.
- [ ] Benchmark complete component returns and cycles, not exposed individual
  property projections. Keep diagnostics/output costs visible.

Acceptance: independent physical tests, full source regressions, model-specific
cross-checks and matched complete-consumer performance pass.

## Phase 5 — Make cycle/model lifecycle reliable and inexpensive

**Status: planned; configured model reuse and separate source setup exist.**

- [ ] Verify every maintained topology owns reusable stream models with correct
  identifiers. Repeated fitness calls must not rebuild/load tables or recheck
  source metadata.
- [ ] Define invalidation when configuration, model method, table domain,
  reference convention or backend changes; rebuilding once on a real setup
  change is correct. Unchanged configuration reloads should avoid needless work.
- [ ] Keep source-only critical/triple/saturation/dilute initialization and phase
  plotting explicit, outside a table's evaluation domain. Handle pseudo-pure
  references honestly and reject models lacking the needed setup physics.
- [ ] Cover independent problems, interleaved designs, configuration reloads,
  copy/worker lifecycle where used and source-reference changes. Do not cache
  mutable provider identity globally without defined invalidation.
- [ ] Confirm initial designs, bounds, objectives and residual ordering across
  the 14 cycle cases; rerun three standard optimizations plus native bicubic
  recuperated optimization and source cross-evaluation.

Acceptance: lifecycle tests pass, repeated evaluations contain no construction
or table loading, and setup/warm/whole-solve timing costs are recorded separately.

## Phase 6 — Finish plotting, exports and normal example workflows

**Status: planned; native phase plotting and numeric Excel cells are verified.**

- [ ] Cover source and bicubic cycle/process/pinch diagrams, reused artists,
  closed figures, optional-field gaps and source-EOS backgrounds explicitly.
- [ ] Standardize numeric Excel/CSV/JSON outputs for NumPy/JAX scalar storage;
  preserve units, labels, identities and unavailable values without new flashes.
- [ ] Test saved configurations and reloaded cases under the new schema;
  exported inspection metadata must not be mistaken for a universal fluid recipe.
- [ ] Verify normal callbacks, reports and headless launchers in addition to
  callback-disabled benchmark runs. Time reporting separately from core fitness.

Acceptance: graphical/data round trips reflect actual evaluated states and the
ordinary user workflow works with the same native objects as the timed cases.

## Phase 7 — Audit and migrate maintained examples/case studies

**Status: planned.** Tests of frozen inputs do not prove every standalone script
works; maintain an explicit execution inventory.

- [ ] Catalogue runnable examples, scientific studies, input-dependent figure
  scripts and historical research material. Record commands, inputs, models,
  expected outputs, execution budget and numerical acceptance criteria.
- [ ] Make minimal new-interface edits in maintained programs. Replace old
  `thermopt.Fluid`, positional backend, `exceptions`, state mutation and state
  collection calls. Preserve useful explanations and physical assumptions.
- [ ] Give the four validated native table studies convenient ordinary YAML
  configurations with justified table domains; source-only phase cases remain
  explicit. Reuse validated caches and configured models.
- [ ] Run bounded representative studies/headless figures and full declared
  scientific profiles separately. Label missing inputs, research limits and
  unexecuted full sweeps honestly. Larger redesigns get their own checklist item.
- [ ] Retire/archive the old 0.5 bicubic adapter and monkey-patching experiment
  after its useful provenance is retained and the native workflow replaces it.
- [ ] Consolidate active comparison code under `benchmark/`, retaining dated
  evidence under `development/`. Keep tests independent of result generation.

Acceptance: every maintained program has a documented successful profile or a
clear external-input/research classification; no blanket "all examples pass".

## Phase 8 — Integrate gates into CI and close the migration

**Status: planned.**

- [ ] Run locked source tests on supported Python versions/platforms and bounded
  native bicubic/component/example checks. Document optional REFPROP selection.
- [ ] Keep timing gates in controlled scheduled/manual jobs, with artifacts and
  repeated flags; do not fail ordinary shared-host CI on one noisy stopwatch.
- [ ] Run the full benchmark release profile and compare against accepted
  checkpoints. Require source cross-evaluation, energy balances, feasibility,
  physical metrics and iteration/call counts for complete optimizations.
- [ ] Update package/API/example docs and migration notes to the supported 0.6
  interface. Remove obsolete compatibility branches and installation caveats
  only when the corresponding checks actually pass.
- [ ] Record final test totals, expected failures, dependency/source/artifact
  hashes, unresolved costs and explicit acceptance. Review each completed phase
  as a coherent change; no release/push occurs merely because a plan is checked.

Acceptance: installed-release workflows, package tests, maintained examples,
numerical agreement and performance gates pass or have explicit scoped decisions.

## Separate follow-ups, not prerequisites for the interface migration

- [ ] Existing co-current total-heat reporting defect.
- [ ] Existing recompression auxiliary-power balance-reporting defect.
- [ ] Existing documented logcosh selector mismatch and smoothing overflow.
- [ ] Exhaustive bicubic table-domain accuracy/resolution studies for additional
  fluids/phase regions, with independent references and declared tolerances.
- [ ] Any future ThermOpt JAX cycle port, differentiable host solver replacement
  or status-packing redesign requires a separate user decision and plan.

## Evidence log

| Checkpoint | Status and evidence |
| --- | --- |
| Prepared 0.6 compatibility | 253 passed, four existing expected failures; [report](development/bicubic/VERIFICATION_060.md), [measurements](development/bicubic/verification_060.json). |
| Published 0.6 dependency/environment | Complete locally: PyPI artifact/hash and 43 source files checked; only JaxProp changes locked version; Poetry 2.2.1 accepts the lock; clean locked environment passes 262 tests plus four unchanged expected failures. [Evidence](benchmark/RELEASE_060.md). |
| Phase baseline benchmarks | Five source/bicubic case or stream modes recorded with 12 samples per backend. [Checkpoint](benchmark/checkpoints/pypi_060_phase.json). One starting cohort, not a between-phase performance acceptance claim. |
