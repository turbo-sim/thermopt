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

**Status: complete locally on 15 November 2026.** The `exceptions` translation is
removed. See [verification and timing evidence](benchmark/PHASE2.md): 319 tests
pass, four existing strict expected failures remain, and all ten three-cohort
aggregate warm comparisons meet the unchanged 1.10 budget. Original and individual
repeat timing flags remain recorded; this is scoped to the tested CPU/x64 profile.

- [x] Migrate maintained example/case YAML and corresponding test input schemas
  to explicit `name`, `method`, `backend`, `identifier`, `table_options` and
  applicable `model_options`. Remove constructor-level `exceptions` settings.
- [x] Remove the legacy translation in `fluid_from_configuration`; reject
  obsolete/unknown fields with useful guidance. Do not accept old positional
  backend spelling or silently reinterpret constructor policies.
- [x] Keep error policy on state calls. Normal ThermOpt solver calls raise on
  invalid states; mask/clip are explicit choices only where the algorithm has
  defined semantics and tests. Never hide a failed optimizer trial with clipping.
- [x] Normalize static YAML scalar settings once and preserve native model
  objects. No state-call wrappers or backend-specific fluid subclasses.
- [x] Cover source, bicubic and perfect-gas construction, malformed settings,
  input immutability and identifiers. Update tests to the new schema without
  changing their physical inputs or expected numerical outputs.
- [x] Add optional process CPU-affinity control/readback to the benchmark protocol
  and its comparison signature. Repeat pre/post with the same updated harness;
  preserve earlier unrestricted observations and declare the protocol transition.
  Windows behavior is exercised; Linux/platform coverage remains Phase 8 work.

Acceptance: maintained cases use the new vocabulary, no old-interface shim is
needed, complete source tests pass and the phase benchmark meets its budget.

## Phase 3 — Standardize state access and capability requirements

**Status: implementation and matched metadata-fix verification complete;
historical absolute timing closure remains separately tracked.** The full suite passes 396 tests plus four
unchanged strict expected failures. All twelve controlled functional complete-cycle
comparisons meet 1.10 (largest increase 4.0%); three historical-checkpoint bicubic
flags remain recorded. Equivalent-model cache identity is a confirmed dependency
confounder, now fixed in JaxProp commit `52c98b8`. The packaged fix passes the
standalone and paired-application gates; its one fresh nitrogen flag does not
repeat in the targeted three-cohort follow-up. Phase 4 proceeds on the user's
instruction while old absolute flags remain Phase 8 verification work. See
[the follow-up](benchmark/METADATA_FIX.md). No clean historical checkpoint is claimed.
See the [inventory](tests/STATE_ACCESS_060.md) and [assessment](benchmark/PHASE3.md).

- [x] Inventory core versus optional fields used by each component, constraint,
  study, plot and exporter. Use canonical property names in maintained code;
  keep friendly display labels separate from field access.
- [x] Resolve required capabilities during component/problem setup. Request
  optional source diagnostics only when supported and needed; a constraint
  requiring an unsupported property must fail clearly, never disappear.
- [x] Use the public immutable-state/index/stack/serialization interfaces.
  Remove dataclass/private-field mutation and obsolete `states_to_dict` helpers.
- [x] Separate unsupported, disabled and locally unavailable properties. Preserve
  statuses, optional NaNs and honest blank report cells; never invent zeros.
- [x] Test scalar/singleton/batch shapes, immutability, source/table identities,
  optional-property failures, source diagnostics and source-model-specific
  quality meaning outside the two-phase region.
- [x] Remove duplicate guarded property evaluation in expression traversal and
  preserve typed JaxProp failures plus expression context. Accept read-only
  mappings without adding metadata checks to state calculation functions.
- [x] Restore canonical-axis auxiliary plotting and replace the active blending
  example's obsolete state-collection call; verify its isolated headless run.
- [ ] Close the historical absolute checkpoint gate in final verification
  (carried to Phase 8 while the matched fix is verified and Phase 4 proceeds).
  Current flags: all-stream recuperated bicubic +10.7%, recompression bicubic
  +13.6%, refrigeration bicubic +22.8%. Preserve the fresh-process and cache-order
  controls; controlled shared-native-model comparisons isolate application changes
  within 4.1%, rather than establishing that those absolute flags have vanished.

Acceptance: same complete state/property behavior through supported routes,
typed informative failures elsewhere, and no new per-property validation work.

## Phase 4 — Verify eager component physics across applicable models

**Status: complete locally on 15 November 2026 for the recorded CPU/x64 profile.**
The full suite passes 461 tests plus the four unchanged strict expected failures;
65 new cases cover the supplied native models and independent physical identities.
Three phase timing cohorts, twelve paired cycle comparisons and both balanced
24-workload component runs pass their warm limits. All five 200x200 release
profiles and both recuperated optimization modes pass their numerical checks.
The JaxProp metadata fix is still an unpublished local wheel; first-use batch
compilation is a Phase 5 follow-up and earlier absolute timing observations remain
Phase 8 work. See [verification and timings](benchmark/PHASE4.md),
[raw evidence](benchmark/PHASE4.json) and the artifact-specific
[checkpoint](benchmark/checkpoints/phase4.json).

- [x] Test compressor/expander isentropic and polytropic paths, heat exchangers,
  valves and functional nondimensional machinery using actual supplied models.
- [x] Preserve SciPy ODE/root solving, scalar states, energy identities,
  endpoint recovery, finite physical outputs and input/model reuse. Keep the
  numerical solver tolerances and acceptance checks unless a defect is proven.
- [x] Expand bicubic component coverage at realistic table boundaries and
  applicable HP/PT/PS/HS routes. Require explicit failure for unavailable paths.
- [x] Add perfect-gas component and applicable eager cycle checks against its
  own ideal-gas identities; do not require it to reproduce real-fluid physics.
- [x] Verify saturation-node placement only with a saturation-capable model.
  Test real two-phase source cases and explicit rejection for unsupported
  PG/table configurations. Do not infer universal phase interchangeability.
- [x] Benchmark complete component returns and cycles, not exposed individual
  property projections. Keep diagnostics/output costs visible.
- [x] Assess whether optional source diagnostics can be deferred from polytropic
  ODE trial states to reporting states, preserving identical returned diagnostics,
  generalized-quality meaning and solver behavior. Measure before changing it.
  The verified implementation defers diagnostics and reconstructs reporting with
  one native full-state batch, retaining finite diagnostic values and optional NaNs.
- [x] Reject failed SciPy integrations explicitly instead of returning partial
  paths as complete components. Add a pre-fix failing regression; check solver
  success once after integration, outside thermodynamic state evaluation.

Acceptance: independent physical tests, full source regressions, model-specific
cross-checks and matched complete-consumer performance pass.

## Phase 5 — Make cycle/model lifecycle reliable and inexpensive

**Status: complete locally on 16 November 2026 for the recorded CPU/x64 profile.**
Phase 4 is committed as `d573638`. The full suite passes 534 tests plus the four
unchanged strict expected failures, including 73 new lifecycle/boundary/cache
cases and the available REFPROP route. Models and source
references are problem-owned; unchanged reloads reuse them after setup-only source
provenance checks. Manager reloads retain their problem identity, rebuild the
solver and roll back rejected updates. Existing workflow identity assertions
are deliberately updated to this documented lifecycle; numerical expectations
and physical reference fixtures are preserved. Three matched full profiles,
twelve paired cycle comparisons and four balanced component runs pass their
aggregate warm budgets. Earlier individual/historical flags remain visible.
Initial setup and new reporting-shape compilation costs remain explicit; this
does not close the historical absolute timing gate. See the
[assessment](benchmark/PHASE5.md), [evidence](benchmark/PHASE5.json) and
[protocol-2 checkpoint](benchmark/checkpoints/phase5.json).

- [x] Verify every maintained topology owns reusable stream models with correct
  identifiers. Repeated fitness calls must not rebuild/load tables or recheck
  source metadata.
- [x] Define invalidation when configuration, model method, table domain,
  reference convention or backend changes; rebuilding once on a real setup
  change is correct. Unchanged configuration reloads should avoid needless work.
- [x] Keep source-only critical/triple/saturation/dilute initialization and phase
  plotting explicit, outside a table's evaluation domain. Handle pseudo-pure
  references honestly and reject models lacking the needed setup physics.
- [x] Cover independent problems, interleaved designs, configuration reloads,
  copy/worker lifecycle where used and source-reference changes. Do not cache
  mutable provider identity globally without defined invalidation.
- [x] Define standalone polytropic output-count and zero-pressure-span contracts.
  With `num_steps=1`, the requested ODE output contains only the inlet, so it
  cannot supply a distinct outlet. Reject unsupported requests clearly at the
  component boundary, with endpoint/failure tests and no new property-call checks.
- [x] Measure whether unreturned isentropic reference states can omit optional
  source diagnostics when only their enthalpy contributes to work. Preserve
  all public inlet/outlet/path annotations; do not change real-fluid equations.
- [x] Isolate first-use compilation of new polytropic batch shapes in fresh
  processes. Phase 4's warm gains coexist with order-sensitive first calls;
  record model setup, inlet warmup, shape compilation and sustained evaluation
  separately. Assess explicit orchestration warmup if useful, without hidden
  flashes or changing the fluid constructor's behavior.
  Normalize SciPy scalar representations once at its boundary to avoid redundant
  JIT signatures. Prepared scalar/reporting shapes are covered by two deterministic
  compilation-reuse tests. A new bicubic batch shape still compiles on first use;
  scalar-only warmup can leave a higher first reporting-call cost than legacy
  scalar stacking. Explicit batch warmup moves this work into setup, not away.
- [x] Confirm initial designs, bounds, objectives and residual ordering across
  the 14 cycle cases; rerun three standard optimizations plus native bicubic
  recuperated optimization and source cross-evaluation.
- [x] Assess and fix the released JaxProp 0.6 bicubic static-metadata equality/cache cost
  with equivalent independent model objects. Phase 3 profiling finds repeated
  dataclass/FrozenMapping equality in JIT lookup; an identical-source control
  is 28.8% slower for the second model, and reversing cache order reverses the
  apparent application regression. Preserve numerical/model independence and
  report the dependency issue upstream; do not patch installed JaxProp or share
  mutable evaluators globally as a ThermOpt workaround. Consider immutable
  metadata interning or efficient equality in a separately tested JaxProp change.
  The aligned equality fix is committed as `52c98b8` with 17 new tests, 3,046
  distinct library tests and packaged ThermOpt verification. The local artifact
  and all residual observations are recorded in `benchmark/METADATA_FIX.md`.
- [x] Record worker hash seed zero in benchmark protocol 2, including schema,
  environment and comparison identity. Verify 19 focused protocol/audit tests
  and capture three matched full profiles with the identical harness on both
  source trees. Protocol-1 checkpoints are historical and incompatible with
  this new baseline; their raw flags are retained. The fixed-seed follow-up
  does not establish hash randomization as the sole cause of previous shifts.

Acceptance: lifecycle tests pass, repeated evaluations contain no construction
or table loading, and setup/warm/whole-solve timing costs are recorded separately.

## Phase 6 — Finish plotting, exports and normal example workflows

**Status: complete locally on 16 November 2026 for the recorded CPU/x64,
headless native-model profile.** The suite passes 598 tests plus the four
unchanged strict expected failures, with 34 new output/boundary cases. Both
ordinary CO2 launchers retain their final objectives and iteration counts.
The actual 200x200 bicubic solution reloads directly through the manager;
source and bicubic snapshot round trips are independently tested. The source
launcher's default does not request YAML exports. A 60-sample balanced check
meets the 1.10 warm limit (ratios 1.012 source, 1.043 bicubic); the initial
12-sample bicubic flag of 1.113 is retained. Reporting costs are separate.
Physical fitness is AST-identical to `bbdfaa3`; no state/property hot-loop checks
were added. See [assessment](benchmark/PHASE6.md) and [evidence](benchmark/PHASE6.json).

- [x] Cover source and bicubic cycle/process/pinch diagrams, reused artists,
  closed figures, optional-field gaps and source-EOS backgrounds explicitly.
- [x] Standardize numeric Excel/CSV/JSON outputs for NumPy/JAX scalar storage;
  preserve units, labels, identities and unavailable values without new flashes.
  The existing conversion helper supports concrete JAX arrays and full precision;
  numeric/boolean/missing-value round trips cover CSV/JSON consumers without
  adding an unsolicited public cycle writer. Excel JAX energy cells are numeric.
- [x] Test saved configurations and reloaded cases under the new schema;
  exported inspection metadata must not be mistaken for a universal fluid recipe.
  Manager/iteration/final YAML is rerunnable; the problem-level writer retains
  formulation-only output. Preserve small tolerances and variable order.
- [x] Verify normal callbacks, reports and headless launchers in addition to
  callback-disabled benchmark runs. Time reporting separately from core fitness.
  Restore accepted designs at output callbacks when probes displaced them;
  save-only plotting creates its required figure. Final exports use x_final.
- [x] Execute the 14 small one/two-case demonstration scripts with fresh outputs
  and their configured report/plot/export callbacks. All 15 optimized cycle runs
  pass native success, feasibility and first-order gates; two property/correlation
  demonstrations also pass. The profile includes the ordinary 200x200 bicubic
  recompression YAML and the ideal/turbomachinery/shaft-constrained PTES cases.
  Interactive watchers and optional movie encoding remain outside this headless
  check. See [the execution inventory](case_studies/SMALL_DEMOS.md).
- [x] Define a consistent optimizer boundary for recognized JaxProp invalid-state
  failures when using controlled SLSQP/Pounce. The shared normalized-input
  boundary preserves native causes/context; direct physical fitness retains
  native errors. Remove the two local comparison subclasses. Unsupported
  capabilities and programming errors propagate; no per-property validation.
- [x] Assess whether maintained constraints need additional exchanger diagnostics.
  None currently requires a new profile opt-in. Existing setup tests reject
  disabled/unsupported diagnostics before reference preparation; retain this
  explicit behavior instead of adding source flashes. A future study requiring
  them must define and measure an explicit setup opt-in.

Acceptance: graphical/data round trips reflect actual evaluated states and the
ordinary user workflow works with the same native objects as the timed cases.

## Phase 7 — Audit and migrate maintained case studies

**Status: in progress.** The small demonstration profile is verified; full
scientific studies and input-dependent tools retain separate open checks.
Tests of frozen inputs do not prove every standalone script works; maintain an
explicit execution inventory.

- [x] Consolidate the former examples directory into `case_studies/`, retaining
  its layout, tracked research inputs and local results. Update active runner,
  test, benchmark and documentation paths; keep historical hashes/records intact.
  This directory migration does not close the remaining study execution audit.
  Verification: 660 tracked renames, all original local files retained, 543 tests
  passed plus four unchanged strict expected failures, and the maintained
  source/bicubic smoke profile passed from the new layout. Explicitly remove
  stale table options when the benchmark selects CoolProp from a bicubic YAML;
  two regression cases cover those source/table routes.

- [ ] Catalogue runnable examples, scientific studies, input-dependent figure
  scripts and historical research material. Record commands, inputs, models,
  expected outputs, execution budget and numerical acceptance criteria.
- [x] Inventory and run the small demonstration subset, preserving initial
  failures and final successful outcomes in `case_studies/SMALL_DEMOS.json`.
  Fix obsolete PHE fluid construction, heat-pump metadata/property-label access
  and plot legends, plus launcher output/headless paths. The complete run found
  two numerical stopping problems: R134a now uses central three-point derivatives;
  turbomachinery PTES uses tighter native stopping (`1e-8`). Original feasibility
  and first-order acceptance gates remain unchanged. Preserve YAML variable order
  in default-solver regression tests; alphabetical dumping had concealed R134a's
  failure. Physical hashes and reference optima remain unchanged. All 43 installed
  JaxProp Python files match the official PyPI 0.6.0 artifact for this profile.
  The full regression suite passes 564 tests plus four unchanged strict expected
  failures. This does not verify full sweeps, interactive watchers or movie encoding.
- [ ] Make minimal new-interface edits in maintained programs. Replace old
  `thermopt.Fluid`, positional backend, `exceptions`, state mutation and state
  collection calls. Preserve useful explanations and physical assumptions.
  All 18 standard SLSQP case/template configurations now explicitly select the
  controlled PySolverView backend with a 0.2 initial fraction, 1.5 growth, and
  accepted-iteration reporting. Preserve per-case tolerances and physical inputs;
  keep named IPOPT and solver-comparison alternatives explicit. Clear inherited
  backend-specific controls before running the deliberate SciPy comparison.
  Verified 564 passing tests plus four unchanged strict expected failures,
  including three complete default-controlled cycle solves and a bounded
  six-configuration SciPy/controlled comparison. The broader research sweeps
  remain outside this verification scope.
- [ ] Give the four validated native table studies convenient ordinary YAML
  configurations with justified table domains; source-only phase cases remain
  explicit. Reuse validated caches and configured models.
  Recompression now has an explicit CO2 bicubic/controlled-SLSQP YAML and a
  runner-defined 128x128 comparison of SLSQP, Pounce interior point and Pounce
  SQP. See `case_studies/sCO2_recompression/README.md`; the broader inventory remains
  open. PySolverView 0.8.2 with Pounce 0.13.0 is a new dependency context:
  previous phase timings remain historical, not automatically comparable.
- [ ] Run bounded representative studies/headless figures and full declared
  scientific profiles separately. Label missing inputs, research limits and
  unexecuted full sweeps honestly. Larger redesigns get their own checklist item.
- [ ] Verify interactive file-watcher and optional movie-encoding workflows
  separately from the completed headless profiles, with declared dependencies.
- [ ] Retire/archive the old 0.5 bicubic adapter and monkey-patching experiment
  after its useful provenance is retained and the native workflow replaces it.
- [ ] Remove the historical benchmark input's `exceptions` normalization when
  consolidating the harness, verifying a new protocol/checkpoint explicitly.
  Phase 2 keeps the old benchmark builder unchanged for comparable measurements.
- [ ] Review study cache identity separately from physical-input schema. Existing
  caches hash complete configurations and correctly reject changed schema/API
  provenance; use a new results directory after migration. Do not rewrite old
  cache signatures or relabel saved results as current-library evaluations.
- [ ] Consolidate active comparison code under `benchmark/`, retaining dated
  evidence under `development/`. Keep tests independent of result generation.
- [ ] Replace the legacy `get_property_fluid` process-global cache in the two
  Rankine study families with explicitly study-owned models and source-change
  handling. Phase 5 removes that helper from configured cycle/reference setup;
  retaining its existing study behavior avoids adding construction or provenance
  probes to those studies' hot loops before their Phase 7 migration.
- [ ] Scope the property-access AST inventory to project-owned source, excluding
  ignored generated benchmark artifacts and frozen comparison trees. Preserve
  coverage of real case_studies/studies. Phase 4's retained reference tree currently
  causes the same historical example SyntaxWarning to be reported twice.

Acceptance: every maintained program has a documented successful profile or a
clear external-input/research classification; no blanket "all examples pass".

## Phase 8 — Integrate gates into CI and close the migration

**Status: planned.**

- [ ] Run locked source tests on supported Python versions/platforms and bounded
  native bicubic/component/example checks. Document optional REFPROP selection.
- [ ] Keep timing gates in controlled scheduled/manual jobs, with artifacts and
  repeated flags; do not fail ordinary shared-host CI on one noisy stopwatch.
- [ ] Exercise optional affinity observation/control on Linux and unsupported
  platforms, including inherited worker masks. Phase 2 verifies Windows on a
  single processor group; do not infer multi-group or cross-platform coverage.
- [ ] Run the full benchmark release profile and compare against accepted
  checkpoints. Require source cross-evaluation, energy balances, feasibility,
  physical metrics and iteration/call counts for complete optimizations.
- [ ] Reassess the carried absolute Phase 2/3 timing flags with stable matched
  artifact/context controls. Keep the original flags and successful metadata-fix
  follow-up separate; matched fix acceptance does not erase earlier ratios.
- [ ] Once separately released, update the dependency/lock to the JaxProp patch
  containing `52c98b8`. The verified local wheel is not an official PyPI artifact.
- [ ] Report the existing direct-JaxProp construction warning for NumPy-valued
  static table options upstream. ThermOpt's configuration boundary already
  normalizes these once; direct `jp.Fluid` tests still reproduce the Equinox
  warning. Prefer constructor normalization in JaxProp, with warning/regression
  coverage, over any checks in thermodynamic calculations.
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
- [ ] In that later sweep, consider preparing expression/residual functions and
  fixed result schemas at setup rather than evaluating configuration strings
  inside numerical calls. Benchmark the design before replacing existing logic.
- [ ] Assess JaxProp complete-state concatenation/reversal in a separate library
  optimization. Phase 5 profiling attributes about 40% of an all-stream cycle
  profile to per-field collection operations; the matched versions perform the
  same 1,728 property calls over 72 cycles. Consider whole-PyTree batching/fusion,
  preserving identities, masks, statuses, optional fields and the NumPy path.
  Require library regressions and cold/warm complete-consumer measurements;
  do not replace full states with property projections.

## Evidence log

| Checkpoint | Status and evidence |
| --- | --- |
| Prepared 0.6 compatibility | 253 passed, four existing expected failures; [report](development/bicubic/VERIFICATION_060.md), [measurements](development/bicubic/verification_060.json). |
| Published 0.6 dependency/environment | Complete locally: PyPI artifact/hash and 43 source files checked; only JaxProp changes locked version; Poetry 2.2.1 accepts the lock; clean locked environment passes 262 tests plus four unchanged expected failures. [Evidence](benchmark/RELEASE_060.md). |
| Phase baseline benchmarks | Five source/bicubic case or stream modes recorded with 12 samples per backend. [Checkpoint](benchmark/checkpoints/pypi_060_phase.json). One starting cohort, not a between-phase performance acceptance claim. |
| Phase 2 configuration | Complete locally: 36 YAML inputs / 108 stream declarations migrated; native objects and read-only configuration supported; obsolete fields rejected. 319 passed plus four unchanged strict expected failures. All ten aggregate warm comparisons meet 1.10 after three matched cohorts; individual flags remain retained. [Assessment](benchmark/PHASE2.md), [checkpoint](benchmark/checkpoints/phase2.json). |
| Phase 3 state access | Implementation and numerical verification complete: 396 passed plus four unchanged expected failures. Twelve controlled functional consumer comparisons pass 1.10; three historical bicubic flags leave the strict gate open. Metadata/cache-order controls and an identical-source false positive are retained. [Assessment](benchmark/PHASE3.md). |
| Packaged metadata fix | JaxProp `52c98b8` passes 3,046 distinct library tests; its verified local wheel passes ThermOpt's full suite and paired application gates. Historical observations remain separate. [Assessment](benchmark/METADATA_FIX.md), [artifact checkpoint](benchmark/checkpoints/metadata_fix.json). |
| Phase 4 component physics | 461 passed, four unchanged expected failures; 65 new native-model/physical cases. Three phase cohorts, twelve paired cycle cases and both balanced 24-component runs pass warm limits. Five 200x200 release profiles and four complete recuperated solves pass numerical/source cross-evaluation checks. First-use shape compilation and historical observations remain assigned follow-ups. [Assessment](benchmark/PHASE4.md), [checkpoint](benchmark/checkpoints/phase4.json). |
| Phase 5 lifecycle | 534 passed, four unchanged strict expected failures; 73 added cases plus 19 focused protocol/audit checks. Unchanged reloads are about 35%/85%/50% faster for source/bicubic/perfect gas. Three matched protocol-2 profiles, twelve paired cycle cases, four balanced 24-component runs and the release profile pass scoped gates. Initial setup, first reporting-shape compilation and older flags remain explicit. [Assessment](benchmark/PHASE5.md), [checkpoint](benchmark/checkpoints/phase5.json). |
| Small demonstration workflows | 14 scripts and all 15 optimized cycle runs pass with fresh outputs and normal configured callbacks; two property/correlation scripts also pass. Actual variable order is preserved in default-solver tests. Minimal API/plot repairs and two explicit numerical stopping fixes retain physical inputs and acceptance gates. Official PyPI JaxProp 0.6.0 sources verified; 564 tests plus four unchanged strict expected failures. [Execution inventory](case_studies/SMALL_DEMOS.md). No broad sweep or performance acceptance claim. |
| Phase 6 output boundaries | 598 tests plus four unchanged strict expected failures; 34 added cases. Rerunnable manager YAML, numeric bicubic energy cells, accepted-design callbacks, native phase backgrounds and shared solver-domain handling pass. Both ordinary CO2 demos retain objective/iteration results; the saved 200x200 bicubic solution reloads. Normalized-input 60-sample ratios 1.012/1.043 meet 1.10; initial 1.113 bicubic flag remains recorded. [Assessment](benchmark/PHASE6.md), [evidence](benchmark/PHASE6.json). |
