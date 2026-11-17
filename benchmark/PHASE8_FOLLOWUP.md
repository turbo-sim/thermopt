# Phase 8: matched timing reassessment

The Phase 7/8 implementation is committed as `746f146`:
**Migrate study workflows to owned native fluids and add release verification**.
This follow-up completes the targeted matched reassessment locally. All six
workloads pass the unchanged 1.10 warm limit in every worker. Phase 8 as a whole
remains open for the external and historical gates below.

The [portable evidence](PHASE8_FOLLOWUP.json) retains all timed samples, source,
input, dependency and artifact identities, individual flags and controls. No
ThermOpt thermodynamic equations, numerical source, physical reference fixtures
or dependency lock changed in this follow-up. No performance workaround or new
validation was added to property calculations.

## What was measured

Compare current ThermOpt `746f146` against the preserved Phase 2 source
`b211fa084b800e970d4bc86f594c02bfc7fef385`, in the same installed environment:
official JaxProp 0.6.0, PySolverView 0.8.2, CPU/x64, hash seed zero and affinity
`0xffff`. All 43 installed JaxProp Python files match the published artifact.
Live official PyPI metadata still reports 0.6.0; the local metadata patch is not
substituted into this comparison or adopted into the lock.

The workloads are the three formerly flagged bicubic complete cycles and their
CoolProp controls: all-stream recuperated CO2, recompression CO2 and nitrogen
refrigeration. Working-fluid tables use 100x100 grids, matching the targeted
historical workload resolution. Supporting stream tables retain the canonical
60x8 settings. Table generation/loading and initial compilation are outside the
warm timings and remain separately recorded.

Eight fresh processes run the existing `compare_versions.py` implementation:

- Four production-style comparisons, with independent native models: both model
  construction/warmup orders, repeated twice.
- Two identical-source controls: both labels run the preserved Phase 2 source
  with separate native models, once per construction order.
- Two functional controls: both source versions receive the same explicitly
  supplied native models through their public functional cycle API, once per
  order. This is an experiment, not global model sharing in production.

Each process retains three bursts of 36 complete evaluations per version and
workload. This gives **10,368 timed complete cycle returns** across 48
worker/workload pairs. All physical snapshots agree exactly across workers;
the functional control additionally checks complete returned state data,
identities, statuses and optional NaNs. Every worker retains source/affinity
stability and no individual worker/workload has a timing flag.

## Results

| Complete problem evaluation | Phase 2 rerun | Current | Current/reference |
| --- | ---: | ---: | ---: |
| Recuperated CO2, all streams / CoolProp | 33.367 ms | 32.992 ms | 0.989 |
| Recuperated CO2, all streams / bicubic | 18.166 ms | 19.327 ms | 1.064 |
| Recompression CO2 / CoolProp | 36.270 ms | 35.918 ms | 0.990 |
| Recompression CO2 / bicubic | 25.441 ms | 25.320 ms | 0.995 |
| Nitrogen refrigeration / CoolProp | 24.292 ms | 24.098 ms | 0.992 |
| Nitrogen refrigeration / bicubic | 16.417 ms | 16.689 ms | 1.017 |

Each version's aggregate is the median of four equally weighted per-worker
medians; the table gives their ratio. The JSON also reports the median of the
paired worker ratios. For all-stream bicubic recuperation these statistics are
1.064 and 1.005 respectively: they are different estimators under varying host
conditions, and both are retained. Neither samples nor controls were discarded
or subtracted, and the 1.10 limit was not widened.

| Workload | Identical-source control ratio | Explicit shared-model functional ratio |
| --- | ---: | ---: |
| Recuperated CO2 / CoolProp | 1.006 | 1.000 |
| Recuperated CO2 / bicubic | 0.990 | 1.016 |
| Recompression CO2 / CoolProp | 1.002 | 0.997 |
| Recompression CO2 / bicubic | 1.000 | 0.998 |
| Nitrogen refrigeration / CoolProp | 1.009 | 0.997 |
| Nitrogen refrigeration / bicubic | 0.992 | 0.998 |

Construction order remains visible: identical-code all-stream bicubic controls
give 1.065 in one order and 0.930 in the other. Current/reference workers likewise
change direction. The functional control's largest aggregate increase is 1.6%.
These observations are consistent with the previously documented independent-model
cache-identity confounder and host variation. They do not establish a new 6.4%
intrinsic cycle-code regression or a precise causal split between those effects.
No numerical-loop optimization is justified by this evidence.

## Historical absolute observations remain separate

The original Phase 2 stored medians have a different protocol, solver environment
and observation time. They cannot supply formal acceptance for the new matched
experiment. The following advisory ratios are preserved rather than dropped:

| Bicubic cycle | Fresh Phase 2 rerun / stored median | Current / stored median |
| --- | ---: | ---: |
| All-stream recuperated CO2 | 1.092 | **1.162** |
| Recompression CO2 | 0.946 | 0.941 |
| Nitrogen refrigeration | 1.023 | 1.040 |

The fresh matched gate passes. The all-stream historical absolute flag remains
open at +16.2%; it is not waived or relabeled as resolved. The older raw flags,
metadata-fix experiments and protocol-specific checkpoints remain unchanged.
This is a targeted reassessment, not a new full release/optimization benchmark
or an unsupported comparison of protocol-3 manifests to protocol-1 clocks.

## Other verification and remaining work

The four new native YAML variants passed their final semantic regression checks
before commit. Their whitespace-only cleanup preserves the evaluated recipes,
with measured and committed byte hashes recorded separately in `PHASE7_8.json`.
The prior 623/622-test correctness checkpoints remain dated evidence; the full
suite was not rerun for this evidence-only follow-up.

A real archive defect was corrected: with `core.autocrlf=true`, Git checkout
filters changed all three retired scripts' recorded byte hashes. A scoped
`*.py -text` attribute now preserves them. Before/after Git-filter checks verify
all three original hashes; script contents are unchanged. This affects archive
reproducibility, not thermodynamic execution.

The visible GUI attempt opened the native 200x200 bicubic Tk plot and produced its
report/pickle artifacts under the case's `results/`. Computer Use then returned
`Computer Use app approval timed out`; recovery found no target window. The
watcher had ended before the planned file edit, so the diagnostic harness's
edit assertion failed. This is an incomplete visual/reload check, not a passing
GUI result or a demonstrated thermodynamic defect. Existing headless live-control,
EOF, report, serialization and actual GIF/MP4 tests remain the applicable verified
coverage. No other Windows application was operated.

Remaining Phase 8 gates are the historical all-stream absolute timing observation,
hosted Windows/Ubuntu CI, actual Linux affinity, publication/adoption of the JaxProp
metadata patch, external submission/fix of the static-option warning, and visible
GUI inspection. No push, external issue submission or package publication occurred.

Raw workers, executable helper snapshots, source archive and artifact receipts
remain under `benchmark/results/phase8_followup_20261117/`. Reproduce with the
existing [paired helper](compare_versions.py), both construction orders, the six
workload selectors and the `--identical-source-control` / `--shared-models` controls
described above. Use fresh output directories and retain all repeats.
