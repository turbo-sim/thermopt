# Phase 4: eager component physics and complete-return performance

Phase 4 passes the scoped CPU/x64 correctness and warm-performance gate. The reference is ThermOpt `6819d34`; both versions use the same verified installed JaxProp wheel containing `52c98b8`. ThermOpt remains eager NumPy/SciPy.

The dependency is an **unpublished local wheel**, SHA256 `bd9d3cdf7ee73ffb096a615f1ed34abd4ced141a18105f09e29c256101e40d71`. Its 43 installed Python files are checked against the artifact. This does not claim a new PyPI release or change the dependency lock.

## Changes and physical verification

- Polytropic ODE trial calls request complete core states without optional source reporting diagnostics. One native full-state batch reconstructs the returned path with its original diagnostics and quality convention.
- Failed SciPy integrations now raise an informative error instead of returning a partial path as a completed component. The new failure test fails before the fix.
- Added 65 cases covering all three native models: machinery energy/efficiency, reversible limits, functional radial solvers, exchanger/valve balances, HP/PT/PS/HS round trips, table edges, explicit unsupported saturation and two complete perfect-gas cycle variants checked against ideal-gas first-law identities.
- Full suite: **461 passed, four unchanged strict expected failures**, zero failures. Existing source/flashing/wet-expansion, saturation snapping and three optimization regressions pass without rewriting references or broadening their tolerances.
- Source audit preserves all function signatures/docstrings, work equations, ODE equations, RK45 settings and six unrelated component/helper functions. All 21 frozen numerical input/reference files are unchanged.

The Air table's lowest-enthalpy corner is trusted gas at low pressure but wet and rejected at high pressure. A rectangular table domain does not imply that every cell or phase operation is available.

## Warm timings

Three fresh 12-sample phase cohorts pass all ten unchanged 1.10 limits against the accepted metadata-fix artifact checkpoint. Twelve same-process complete-cycle comparisons retain exactly identical physical snapshots; their largest warm increase is 3.6%. The complete polytropic cycle is 22.0% faster with CoolProp and 43.3% faster with bicubic.

The component helper covers 24 full-return workloads, using independent reusable native models and both construction orders. Each fresh process retains three balanced cohorts of 12 samples. Both runs pass every warm limit. Complete output fields, optional NaNs, masks, statuses and identities agree (numeric rtol $2\times10^{-10}$, atol $10^{-7}$; discrete metadata/status exact).

| Model | Compressor time reduction | Expander time reduction |
| --- | ---: | ---: |
| CoolProp | 24.0–24.1% | 47.6–49.5% |
| Perfect gas | 51.5–56.4% | 86.6–87.5% |
| Bicubic | 38.7–52.7% | 77.6–82.1% |

These are default complete paths (10 compressor / 50 expander nodes), including state construction and reporting. Scalar re-flashing and per-field stacking disappear from reporting; the source ODE also avoids unused saturation diagnostics. No backend adapter, property projection or suppressed validation is used. Apparent speedups in unchanged controls are not attributed to a code change.

A preliminary four-sample harness smoke recorded 6 control flags. Neither sustained balanced run repeats them. Two earlier harness attempts failed while comparing string/heterogeneous reference metadata; those attempts precede timing acceptance and remain in the execution log. The corrected helper compares nested metadata structurally and numerical arrays numerically.

## First-use costs and limitations

First component calls are recorded separately and remain order-sensitive. For example, bicubic polytropic expansion takes about 206–244 ms on first invocation versus 44–87 ms for the reference in these processes, while warm calls are 7.7–8.8 ms versus 39.5–43.2 ms. New batch shapes can compile once. These observations follow inlet/prior-workload warmup and do not isolate fresh-process compilation; a cold-JIT no-regression claim is not made. Fresh isolated compilation measurement is recorded in Phase 5, together with singleton/zero-span component contracts and potentially unused isentropic-reference diagnostics.

Earlier absolute Phase 2/3 timing flags remain visible in the [metadata follow-up](METADATA_FIX.md) and final Phase 8 verification. Passing this matched phase does not erase them or establish universal table accuracy.

## Release profile

All five 200×200 case/stream profiles pass their existing source/table numerical budgets. Source and bicubic recuperated optimizations pass for working-fluid-only and all-stream table modes, including independent source cross-evaluation, feasibility and objective agreement. Iteration/call counts and complete-solve times are recorded separately; these are verification evidence, not a matched between-phase optimization speedup claim.

| Mode | Model | Solve [s] | Iterations | Calls | Constraint violation |
| --- | --- | ---: | ---: | ---: | ---: |
| recuperated_working | source | 5.642 | 22 | 168 | 1.71e-14 |
| recuperated_working | bicubic | 4.393 | 24 | 184 | 8.53e-15 |
| recuperated_all | source | 5.745 | 22 | 168 | 1.71e-14 |
| recuperated_all | bicubic | 3.000 | 22 | 168 | 5.68e-15 |

Full raw samples, environment/input/source/harness identities, numerical comparisons, solve results and retained preliminary observations are in [PHASE4.json](PHASE4.json). The [checkpoint](checkpoints/phase4.json) is scoped to the recorded local wheel, CPU/x64 environment and protocol. Reproduce with the commands in [README.md](README.md), using new output directories. Raw local artifacts remain under `benchmark/results/phase4/`.
