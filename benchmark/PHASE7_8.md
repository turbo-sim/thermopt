# Phase 7 case studies and Phase 8 local verification

Phase 6 is committed as `8341a03` (`Make native-model reporting reproducible and
share solver-domain recovery`). Phase 7 is complete for its declared bounded
migration profile. Phase 8's CI and local verification are implemented; external
release, hosted-platform and historical timing gates remain open. The machine-readable
[evidence](PHASE7_8.json) preserves identities, results, raw study timing samples
and the limitations below. No release or push was performed.

## Case-study migration

The [execution inventory](../case_studies/EXECUTION_INVENTORY.json) covers all 47
project-owned Python programs: 15 small demonstrations, six parameter sweeps,
four scientific drivers, ten saved-results analyses, two saved-input exporters,
two solver benchmarks, one mesh study, four helpers and three historical figure
programs. Commands, inputs, working directories, budgets and acceptance criteria
are explicit. Classifying a research program does not establish a successful
execution of its entire grid.

The two Rankine study families now use explicitly owned native source models.
The former global `get_property_fluid` cache is removed. The transcritical study's
bounded water caches belong to its problem and reset when that water model is
replaced. Repeated numerical evaluations neither construct models nor probe source
provenance. Setup/report operations receive or create their own native models.
Existing study equations, explanations and candidate reevaluation remain intact.

Study cache signatures now cover configuration/template, core/study code,
dependencies, active EOS reference conventions and solver/retry acceptance settings.
Wet/dry saved initializers are fingerprinted by content. Current design values
remain intentionally excluded for continuation. Old records are rejected on a
provenance mismatch, never rewritten as current results. Shared readers accept
complete manager snapshots and formulation-only exports without adapting old
fluid settings. Recoverable PySolverView `DomainError` trials are counted as failed
attempts instead of aborting a multi-start study. Typed JaxProp configuration and
unsupported-request errors bypass the studies' broad legacy numerical-error
handlers, so another initial guess cannot conceal an unsupported model choice.

New wet/dry outputs go under `results/studies`, `results/logs` and `results/report`;
historical `output` records and supplied inputs remain in place. PTES table tools
accept explicit input/output roots and are exercised against native exports.
The two scientific report pipelines and their full grids were not regenerated.

Four convenient YAML variants expose the already validated native choices. Three
add 200x200 working-fluid tables to source-backed cases; recompression gains an
explicit source input alongside its default bicubic input. Physical settings and
solver controls match their reviewed counterparts. Runners accept `THERMOPT_CONFIG`
and enable x64; no application adapter or backend fallback is introduced.

## Additional complete launcher runs

The actual launchers ran with normal iteration limits, callbacks, plots and exports.
Only table storage paths were redirected to the existing validated benchmark cache.
Measured recipe hashes are retained. A final trailing-whitespace/blank-EOF cleanup
of the four new YAML files preserves their parsed configurations exactly; their
committed hashes are recorded separately in the JSON. Each run has fresh output below
its case's `results/phase7_checks/20261117_084539`.

| Run | Evaluator | Iterations | Objective (minimization convention) | Maximum constraint violation |
| --- | --- | ---: | ---: | ---: |
| Recuperated CO2 | Bicubic 200x200 | 23 | -0.382689250185 | 5.00e-13 |
| Simple/polytropic CO2 | Bicubic 200x200 | 25 | -0.220887502724 | 0 |
| Recompression CO2 | CoolProp HEOS | 29 | -0.413819874551 | 7.91e-9 |
| Nitrogen refrigeration | Bicubic 200x200 | 18 | -0.696864201436 | 1.87e-7 |
| PTES system-scale, first 1 MW point | Source EOS, shared-shaft machinery | 108 | -0.448488996554 | 2.05e-10 |

All five passed native success, finiteness and the configured feasibility gates.
PTES's `THERMOPT_SMOKE=1` selects its first existing power point while retaining
the full solve; the complete 21-point sweep was not executed. The prior 14-script,
15-solve [small-demo checkpoint](../case_studies/SMALL_DEMOS.md) remains dated
evidence and is not relabeled as a fresh run of every launcher.

## Numerical and performance evidence

All component/cycle module ASTs and the core physical `fitness` AST match Phase 6.
Original `tests/data` physical inputs/reference fixtures are unchanged. Study
ownership changes pass native numerical parity and lifecycle tests, including no
model construction/source probes during repeated evaluations and changed-reference
cache invalidation.

`benchmark/study_models.py` compares the real full Rankine study `fitness` against
the Phase 6 tree in four fresh processes, reference/current/current/reference,
with 60 samples each. Complete residual vectors agree within 1e-12.

| Full study evaluation | Phase 6 median | Current median | Ratio |
| --- | ---: | ---: | ---: |
| Native butane/source-water cycle | 37.738 ms | 36.996 ms | 0.980 |

This shows comparable sustained cost for the ownership change. It is one balanced
observation, recorded before final output-only core edits and setup/provenance
additions; the numerical evaluator remains unchanged. It is not a general
final-source performance comparison or evidence that old flags disappeared.
Raw processes and samples are retained under
`benchmark/results/phase7_study_timing_complete/`.

The final complete release profile uses installed official JaxProp 0.6.0, CPU/x64,
200x200 tables and 12 samples per model/workload. All five profiles pass source
cross-evaluation and physical/numerical gates. Initial-point maximum temperature
error is below 7.8e-5 K across them. Both recuperated table modes also pass complete
source/table optimizations: source 24 iterations, bicubic 23, objectives agreeing
to approximately 1.2e-9. Raw call counts, optimum cross-checks, balances and settings
remain in the workload summaries.

| Complete warm cycle | CoolProp HEOS | Bicubic |
| --- | ---: | ---: |
| Recuperated CO2, working stream tabulated | 44.26 ms | 42.84 ms |
| Recuperated CO2, all streams tabulated | 57.63 ms | 34.48 ms |
| Simple/polytropic CO2 | 177.01 ms | 83.37 ms |
| Recompression CO2 | 67.77 ms | 52.72 ms |
| Nitrogen refrigeration | 45.41 ms | 26.92 ms |

These are backend observations in this environment. They do not compare current
ThermOpt against Phase 6. Earlier same-phase release observations are retained in
the JSON: source/table working-stream recuperation was 33.59/23.43 ms, and
recompression 53.09/53.88 ms. This variation includes both source and table controls;
neither candidate alone establishes a persistent cost or universal bicubic speedup.
Setup, generation/loading and complete solves remain separate measurements.

The canonical implementation now lives in `benchmark/native_backends.py`; the
obsolete `exceptions` stripping is removed. Main-runner protocol 3 records that
change. Solver comparison implementations also live in `benchmark/`, with thin
case-local launchers. The old adapter experiment is archived byte for byte with
verified original hashes. Dated evidence is unchanged. The final release manifest
is `benchmark/results/phase7_8_release_verified/manifest.json`; source, harness
and lock hashes stayed fixed during the run and still match current numerical
source. It is a **baseline candidate**, not an accepted comparison to incompatible
protocol-2/local-patch/older-solver checkpoints. Historical absolute flags remain.

## Local installation, tests and optional workflows

- Python 3.13.2: **623 passed, four unchanged strict expected failures**.
- Clean locked main dependencies on Python 3.11.9: **622 passed, one optional
  REFPROP skip, four unchanged strict expected failures**. The existing local
  Python 3.13 environment does exercise REFPROP. The clean environment does not
  imply a missing proprietary library is a package failure.
- Final focused checks after comparison relocation and additional restart/cache
  assertions passed (31, 21 and 14 cases respectively; overlapping counts).
  The complete runs precede these final changes, which are covered by the focused
  checks. Original strict expected failures and their reasons are preserved.
- A built ThermOpt 0.2.5 wheel was installed outside the repository import path
  in the isolated Python 3.11 environment. All 23 installed source files match
  this tree, with native source and bicubic compressor evaluations passing.
  This is packaging QA, not a new package release.
- All 43 installed JaxProp Python files match the preserved official 0.6.0
  artifact. Live PyPI metadata checked on 17 November still reports 0.6.0; the
  local metadata fix `52c98b8` has not been released. Dependency/lock unchanged.
- Live watcher tests use a real input thread with headless plots: changed designs,
  requested reports, initial rendering, EOF exit and post-watcher pickling pass.
  Input threads are excluded from serialization. Ordinary launchers reload manager
  configuration after interactive edits so the solver receives those changes.
  Visible GUI/window painting was not manually observed.
- Tiny GIF and MP4 files are actually encoded/decoded for frame order/count/rate.
  GIF seconds are converted to the Pillow writer's milliseconds; manager `fps`
  now applies to both formats. These changes are confined to I/O boundaries.

Full pytest logs/JUnit are under `benchmark/results/phase7/python313_final.*` and
`benchmark/results/phase8/python311_final.*`. The wheel receipt and hash are under
`benchmark/results/phase8/wheel_qa_verified.json`; identities are also preserved in this
report's JSON so essential results do not depend on committing generated folders.

## Phase 8 remaining gates

CI is prepared for Windows/Ubuntu and Python 3.11/3.13: locked production
dependencies, full tests, bounded native numerical smoke and retained artifacts.
A manual workflow records selectable benchmark profiles without imposing an
unrelated shared-runner stopwatch baseline. Workflow policy and Linux/unsupported
affinity branches have unit coverage. Hosted jobs and actual Linux kernel
affinity/inheritance have not run; no push was authorized or performed.

The [local upstream repro](JAXPROP_STATIC_OPTIONS.md) captures the direct-JaxProp
NumPy-valued static-option warning. ThermOpt already normalizes at setup. The
report has not been submitted externally, and the upstream constructor fix plus
published metadata patch remain pending. Historical timing closure needs a stable
matched artifact/context comparison. Full research sweeps, supplied-input reports
and visible GUI checks retain their explicit separate scope. The four existing
physics/reporting/math expected failures remain independent follow-ups.
