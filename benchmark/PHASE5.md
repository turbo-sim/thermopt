# Phase 5: reusable model lifecycle and host-solver performance

Phase 5 passes the scoped matched warm gate. The numerical reference is ThermOpt Phase 4 commit `d573638`. Both implementations use the verified installed JaxProp wheel containing `52c98b8`; ThermOpt remains eager NumPy/SciPy. JaxProp source code is unchanged in this phase.

## Behavior and correctness

- Models belong to a problem. Identical reloads reuse native evaluators and immutable reference points. A changed fluid, backend, method, domain, calibration or actual source convention rebuilds affected roles. Source provenance is checked only at setup.
- Reference evaluators are independently owned; configured cycles no longer use the process-global study cache. Pseudo-pure lower endpoints retain their actual source interpretation. Source reference queries remain outside interpolation domains.
- Configuration is copied and prepared before being applied. Rejected model, constraint, expression, initial-evaluation and batch updates retain the last usable problem/solver. Successful reloads invalidate old results; plot/Excel entry points explicitly require a new fitness evaluation. Saved configurations reflect the newly applied settings.
- Manager reloads now retain `manager.problem` identity and replace `manager.solver` to respect new dimensions. Discard external references to an old solver after reload. The two identity assertions in workflow tests deliberately change; their physical expectations remain intact.
- Polytropic output counts require integers >=2. Zero pressure spans return constant full-state paths without integration. Unreturned isentropic reference states omit unused reporting diagnostics; public inlet/outlet/path diagnostics remain available.
- Host-solver inputs use consistent NumPy scalar representations. SciPy passes its initial pressure as a Python float and later pressures as NumPy scalars; normalize that boundary to avoid duplicate JIT signatures. Density arithmetic for the host ODE uses a host scalar. All property requests still return native complete states.
- **534 tests passed; four unchanged strict expected failures.** There are 73 new lifecycle/boundary/cache cases, including HEOS/available REFPROP, an actual external reference-convention change, interleaving, copy/pickle restore and scalar cache reuse. The two cache tests assert compilation reuse, without stopwatch limits.
- Source audit: nine cycle modules and six unrelated component/helper functions are unchanged. ODE/work equations, RK45 arguments, result schemas and the fitness body are preserved; all 21 physical input/reference files are unchanged. No golden values or numerical tolerances are regenerated.

The test-only reference control imports exactly `set_reference_state` from CoolProp to change external source configuration; all thermodynamic evaluation remains in JaxProp. The property-access audit permits only that exact import in that test module.

## Reload timings

Three balanced cohorts of twelve samples per version include all provenance/rendering work in reload timing, followed by a complete synchronized cycle evaluation. All returned states, statuses, optional NaNs, masks, model identities and physical outputs are checked outside timings. Different backends are compared with their own baseline, including an explicitly specified ideal-gas Air cycle.

| Model / update | Reload reduction | Reload + complete evaluation reduction |
| --- | ---: | ---: |
| coolprop/unchanged | 34.9% | 17.5% |
| coolprop/nonfluid | 35.4% | 19.3% |
| bicubic/unchanged | 84.9% | 68.3% |
| bicubic/nonfluid | 85.1% | 67.3% |
| perfect_gas/unchanged | 50.1% | 26.0% |
| perfect_gas/nonfluid | 48.8% | 26.0% |

Initial construction remains a separate cost. In the recorded source case setup increases about 32→48 ms, and perfect-gas setup about 26→34 ms, including the new provenance/ownership work. Same-process first-evaluation records are cache-order dependent and are not used as cold-JIT speedup evidence.

## Complete warm consumers and retained flags

All twelve paired complete-cycle comparisons pass 1.10; their largest ratio is 1.0176. Physical snapshots remain exactly identical. Component returns cover 24 workloads with independent native objects. Both construction orders are run twice; the aggregate is the median of those four per-run medians, with unchanged 1.10 limits.

Individual order-dependent flags remain recorded, including the unchanged valve control. Reversing model construction/warmup order changes those ratios; the balanced assessment avoids treating cache population as an application-code change. No mutable evaluators are shared globally to mask the cost.

| Balanced component workload | Current/reference time |
| --- | ---: |
| coolprop/compression_process/polytropic | 0.9842 |
| coolprop/expansion_process/polytropic | 0.9924 |
| perfect_gas/compression_process/polytropic | 0.9386 |
| perfect_gas/expansion_process/polytropic | 0.9963 |
| bicubic/compression_process/polytropic | 0.8686 |
| bicubic/expansion_process/polytropic | 0.9236 |

The initial random-seed fresh-process profiles included a persistent all-stream bicubic flag. A separately declared three-cohort/36-sample fixed-seed repeat gives 0.9964×, with no flags. A concurrent JaxProp task was observed; that observation does not establish it or hash randomization as the sole cause. Independent complete-cycle profiles show the same 1,728 property calls and 23,112 NumPy conversions over 72 evaluations, taking about 2.194 versus 2.206 seconds with profiling enabled.

**Benchmark protocol 2** therefore fixes worker `PYTHONHASHSEED=0` and records it in schema, environment and comparison identity. Both numerical source trees use the identical updated harness. Three matched full profiles establish the new checkpoint; protocol-1 checkpoints are explicitly incompatible. Old manifests, raw samples, flags and the harness transition remain in the evidence.

| Protocol-2 complete cycle | Current/reference time |
| --- | ---: |
| recuperated_working/100/source | 0.9619 |
| recuperated_working/100/bicubic | 0.9985 |
| recuperated_all/100/source | 1.0440 |
| recuperated_all/100/bicubic | 1.0011 |
| simple_polytropic/100/source | 1.0414 |
| simple_polytropic/100/bicubic | 1.0412 |
| recompression/100/source | 1.0485 |
| recompression/100/bicubic | 1.0307 |
| refrigeration/100/source | 0.9832 |
| refrigeration/100/bicubic | 1.0114 |

## First-use compilation

Every model/machine/warmup/cohort/version is a fresh CPU/x64 process with persistent compilation caching disabled. Import, model construction, inlet PT evaluation, explicit warmup, first complete component and twelve sustained calls are separate. Three cohorts compare the final code against the older scalar-reporting source; another matched experiment isolates the host-boundary change.

| Backend/process, no optional warmup | Before host fix [ms] | Current [ms] |
| --- | ---: | ---: |
| perfect_gas/compression/none | 383.7 | 307.4 |
| perfect_gas/expansion/none | 373.6 | 304.3 |
| bicubic/compression/none | 1539.9 | 1231.2 |
| bicubic/expansion/none | 1412.3 | 1306.7 |

Compiler logs identify the redundant scalar representation; the regression tests verify that prepared scalar/reporting shapes need no extra `get_state` compilation inside the SciPy process. New reporting shapes still compile once. When scalar kernels are already warm, bicubic reporting batches can have a larger first call than the legacy scalar/stack path:

| Backend/process, scalar-only warmup | Legacy [ms] | Current [ms] |
| --- | ---: | ---: |
| perfect_gas/compression/scalar | 151.5 | 174.6 |
| perfect_gas/expansion/scalar | 200.4 | 175.1 |
| bicubic/compression/scalar | 130.7 | 266.6 |
| bicubic/expansion/scalar | 201.7 | 283.8 |

Explicit full-state batch warmup moves shape compilation into setup; it does not remove the work or run automatically in fluid construction. The first-use runner has its own protocol 2, warming all scalar representations used by final/legacy code. Earlier protocol-1 investigations are retained but their scalar/path warmup results are not interchangeable checkpoints.

## Release checks and next work

The protocol-2 release profile passes all five 200×200 working-grid case/stream profiles and source/bicubic recuperated optimizations in working-only and all-stream modes. Independent source cross-evaluation, feasibility and objective agreement retain the existing budgets. Whole-solve times, iterations and calls are recorded separately; no unmatched optimization speedup is claimed.

Remaining work is explicit in the implementation checklist: plotting/export/example workflows, study-local replacement of the legacy shared helper, published adoption of the metadata fix and final historical/platform checks. Profiling also identifies per-field JaxProp state concatenation/reversal as a future measured optimization; changing those operations belongs in a separately tested library change.

The dependency remains the **unpublished local JaxProp wheel**, SHA256 `bd9d3cdf7ee73ffb096a615f1ed34abd4ced141a18105f09e29c256101e40d71`. All 43 installed source files match it; this phase does not claim publication or modify the dependency lock.

Exact test/benchmark/source/environment/artifact identities, every raw component sample and the retained diagnostic stages are in [PHASE5.json](PHASE5.json). Local full worker outputs/logs remain in `benchmark/results/phase5/`. Commands are documented in [README.md](README.md).
