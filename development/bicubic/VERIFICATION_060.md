# ThermOpt verification with JaxProp 0.6.0

**Subsequent release update, 15 November 2026:** JaxProp 0.6.0 is now published
and ThermOpt's lock/environment have been updated. See the
[migration plan](../../z_JAXPROP_060_IMPLEMENTATION_PLAN.md) for the separate
published-wheel checks. Measurements below retain their original prepared-wheel
provenance and are not relabeled as new PyPI measurements.

The prepared 0.6.0 wheel works with the eager ThermOpt cycle implementation after the small interface changes recorded in the [workflow/checklist](README.md). Native bicubic models work without application fluid adapters, numerical fallback or a JAX cycle port.

**Final suite: 253 passed, zero failed/skipped, four unchanged expected failures.** All seven final validation jobs pass on package source `a0abd437d975e65a9103f0e12bd7bb2e7285606d17d79b426b7ef5d8b95552a6`. The installed package matches all 43 Python files in wheel `jaxprop-0.6.0-py3-none-any.whl` (SHA256 `5031e6bd7e00dc745fd5998e836e55666d0e243d2cfe8dc228b84a278dd374f0`).

The suite covers 14 cycle cases at initial/nearby frozen physical designs, three full source-backend optimizations, headless plots/exports, functional machinery and five added bicubic integration tests. Original physical-state/energy/objective references and their strict tolerances remain. Separate versioned records capture corrected setup references, 36 generalized-quality annotations and two initial R134a diagnostic residuals, independently verified from source states/saturation formulas.

## Recuperated CO2 optimization

| Model | Efficiency, % | Iterations / calls | Warm solve, s |
| --- | ---: | ---: | ---: |
| CoolProp HEOS | 38.268925134 | 22 / 168 | 5.660 |
| Bicubic working fluid, 100x100 | 38.268916117 | 24 / 191 | 4.219 |
| Bicubic working fluid, 200x200 | 38.268925019 | 24 / 184 | 3.936 |
| Bicubic all three streams, CO2 200x200 | 38.268925018 | 23 / 176 | 2.578 |

Each model converges with the example's original solver settings. The returned table optima also pass source evaluation: normalized source constraint violation is at most $2.6\times10^{-10}$. The 200x200 maximum temperature/profile/pinch discrepancy at matched optimized designs is about $3.2\times10^{-5}$ K; the 100x100 discrepancy is about $3.2\times10^{-4}$ K. These are sampled initial/optimized-design checks, not a whole-domain accuracy certificate.

The all-stream source solve is timed separately at 5.186 s. Timings vary by workload context. Source and table optima have different iteration/call counts; solve time is not a property-kernel-only benchmark.

## Complete eager cycle evaluations

| Case | CoolProp, ms | Bicubic, ms | Table / source | Maximum temperature error at initial design, K |
| --- | ---: | ---: | ---: | ---: |
| Recuperated CO2, working fluid | 34.746 | 24.690 | 0.711 | 1.99e-06 |
| Recuperated CO2, all streams | 30.498 | 16.714 | 0.548 | 4.56e-06 |
| Simple/polytropic CO2 | 129.689 | 78.099 | 0.602 | 2.05e-06 |
| Recompression CO2 | 35.636 | 24.274 | 0.681 | 7.77e-05 |
| Nitrogen refrigeration, all streams | 22.863 | 15.034 | 0.658 | 3.15e-06 |

All table comparisons pass the separately declared interpolation allowances. Timings use 12 alternating synchronized samples for each backend after warmup on a shared CPU; raw samples are retained. No other local validation job runs concurrently with these timers. Other user/agent activity is uncontrolled.

Fresh CO2 table generation took 66.27 s (100x100) and 268.98 s (200x200) in the first run; these costs are separate from the warm timings. Validated cache reuse is automatic. Scalar YAML options are normalized once at model setup, outside every repeated state call.

## Normal workflow and limits

The native bicubic cycle renders its source-EOS phase background and cycle/pinch curves successfully. Excel export now writes native scalar properties as numeric cells; unsupported optional diagnostic columns remain blank. A new regression catches the prior JAX-scalar-as-text behavior.

CoolProp phase-changing Rankine/ORC cases pass the source suite. A single bicubic table is not promised to replace every phase/diagnostic route. Unsupported operations and out-of-domain states raise explicitly, with no hidden fallback or clipping.

At this measurement checkpoint, the package dependency targeted `^0.6.0` while the lock still selected 0.5.7; the lock was subsequently refreshed against the published package. The earlier description of a missing JaxProp lock entry was incorrect. Legacy exploratory HX/figure scripts and four unrelated expected failures remain in the checklist. No REFPROP, GPU, exhaustive table-domain study or full parametric scientific sweep is claimed.

The existing dirty worktree is preserved; no commit was created and JaxProp/Barotropy source was not edited. [Portable data and complete timings](verification_060.json); full raw logs, YAMLs, histories and generated tables remain in ignored `results/compatibility_060/`.
