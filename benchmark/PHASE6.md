# Phase 6: plots, configurations and output boundaries

The starting checkpoint is `bbdfaa3`, following solver-default commit `2d8a05c`.
Thermodynamic equations and the physical-input/reference fixtures are unchanged.
This phase adds reporting and solver-interface checks around the native objects.

## Changes and verified behavior

- Manager snapshots, iteration YAML and `optimal_solution.yaml` contain both
  `solver_options` and the evaluated `problem_formulation`. They can initialize
  `ThermodynamicCycleOptimization` directly. The problem-level writer retains its
  formulation-only contract. Fluid recipes remain construction inputs; critical
  points, evaluated reference metadata and table contents are not substituted.
- Configuration saving preserves variable order and full float precision. The
  common conversion helper accepts concrete NumPy/JAX arrays and supports
  `precision=None`; small solver tolerances such as `1e-14` remain nonzero.
  Numeric/boolean types and missing values survive YAML/JSON/CSV conversion.
  This does not introduce a separate public cycle CSV/JSON writer.
- Excel energy-analysis cells now contain native numbers for bicubic/JAX
  results, matching the state-property cells. Units, column labels and endpoint
  identities remain unchanged. Unsupported diagnostics stay blank; output never
  triggers additional thermodynamic property evaluations.
- Final result saving refreshes the solver's final design before writing.
  Iteration callbacks restore their supplied physical design only when a probe
  displaced it; later observers reuse that evaluation. Report diagnostics restore
  it on exit. These checks occur at output boundaries, outside repeated fitness
  and state calls. Save-only plot callbacks now automatically create a figure.
- Source and native bicubic process/pinch diagrams are checked numerically,
  including artist reuse and closed-figure recreation. Bicubic saturation
  backgrounds use JaxProp's source EOS outside the interpolation table domain.
- The normalized-input solver boundary translates recognized JaxProp invalid
  states/calculation failures into PySolverView `DomainError`, preserving the
  native cause and context. Controlled SLSQP/Pounce can reject invalid trials.
  Unsupported requests, configuration errors and ordinary programming defects
  propagate. Direct physical `fitness` calls retain native JaxProp errors. The
  comparison scripts use this shared boundary and no longer define adapters.

Pre-fix failures are retained under `benchmark/results/small_demos/phase6_*before.log`.
The configuration/Excel/serialization regressions and four callback failures
were demonstrated before their fixes. Synthetic error-routing tests are
supplemented by actual invalid source and bicubic cycle trials.

## Reproduction and scope

```powershell
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe benchmark/output_workflows.py --output benchmark/results/phase6_manual
```

The focused timing harness compares complete normalized cycle evaluations
through the inherited pre-Phase-6 PySolverView route and the new exception
boundary, alternating order on one configured native model. This isolates the
boundary cost from independent-model metadata/cache order. It separately times
YAML, Excel, plot updates, the first KKT report and cached reports. Report elapsed
times include output-only design synchronization; solver counters do not count
those direct refreshes. This is not a comparison of different JaxProp releases.
The default 32x32 table bounds this I/O test; the ordinary recompression launcher
is separately verified with its configured 200x200 table.

The final full suite, native CO2 launcher results, snapshot reload checks and
bounded timing observations are recorded in [PHASE6.json](PHASE6.json).
The full suite passes **598 tests with four unchanged strict expected failures**,
including 34 new checks. The normal source and 200x200 bicubic CO2 launchers
retain their exact objectives and 24/29 iteration counts. The saved bicubic
solution reloads with the same objective and numeric Excel energy cells. The
source launcher's default does not write YAML; its manager snapshot contract
is covered by the independent native-source tests.

| Complete normalized cycle route | Inherited | Current | Ratio |
| --- | ---: | ---: | ---: |
| CoolProp HEOS | 43.068 ms | 43.579 ms | 1.012 |
| Native bicubic | 26.917 ms | 28.076 ms | 1.043 |

These are the 60-sample balanced follow-up medians, both within the declared
1.10 limit. The initial 12-sample check recorded 0.915/1.113; its bicubic flag
is retained in the evidence rather than erased. This bounded check does not
establish a persistent 1.2%/4.3% regression or close older absolute timing gates.

The earlier [14-demo inventory](../case_studies/SMALL_DEMOS.md) remains a dated
checkpoint rather than being relabeled as a run of this new source version.

No current maintained default needs extra exchanger-profile diagnostics in
constraints. Existing setup checks reject disabled/unsupported requests before
reference preparation. New diagnostic opt-in remains a future study-driven
decision, rather than silently adding source flashes to ordinary profiles.
Interactive file watchers, movie encoding, full scientific sweeps and historical
figure pipelines retain their explicit Phase 7/8 verification scope.
