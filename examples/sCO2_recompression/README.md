# Recompression sCO2 cycle

The normal YAML uses native JaxProp bicubic interpolation for CO2 and
PySolverView's controlled SLSQP (`library: pysolver`, `method: slsqp`). The first
step is capped at 0.2 of the proposed direction; the cap grows by 1.5 toward 1.
Heating nitrogen and cooling water retain their CoolProp HEOS models.

Install the updated locked dependencies from the repository root:

```powershell
poetry install
```

PySolverView 0.8.2's `pounce` extra supplies Pounce 0.13.0. Both are included in
ThermOpt's dependency declaration and lock. The YAML's normal CO2 table is
200x200, with $150\leq h\leq1300$ kJ/kg and $8\leq p\leq42$ MPa. The pressure
range includes upstream pressure-drop margins; its lower enthalpy avoids an
invalid source-EOS corner at high pressure. This domain is for the configured
cycle, not a claim that arbitrary trial designs or phase operations are valid.

For the interactive single-run example, launch from this directory:

```powershell
poetry run python run_optimization.py
```

The initial cycle plot watches the YAML. Press Enter in the terminal to start
optimization, or close the figure. Its callbacks save the normal cycle outputs.
Generated cycle outputs and the property-table cache live under `results/`.

To compare all three algorithms with both CoolProp and bicubic, run from the
repository root:

```powershell
poetry run python examples/sCO2_recompression/compare_solvers.py
```

For a headless run:

```powershell
poetry run python examples/sCO2_recompression/compare_solvers.py --no-show
```

The comparison includes controlled SLSQP, Pounce interior point (`ipm`) and
Pounce active-set SQP (`sqp`), once with CoolProp HEOS and once with bicubic.
The runner selects the working-fluid evaluator in memory, using a **128x128
bicubic grid** and an ordinary native CoolProp model without table options;
it never writes the input YAML. Auxiliary fluids stay on CoolProp in both cases.
Physical initial designs, bounds, constraints, scaling and finite differences
are identical. Both launchers enable 64-bit JAX table arithmetic. Pounce uses
the common feasibility/optimality assessment for early KKT stopping and its own
line search; SLSQP step-fraction options are
not sent to Pounce. Pounce SQP's inner-QP feasibility and optimality tolerances
are `1e-7`: its default `1e-9` settings stalled at the initial QP in the validation
run. The outer feasibility (`1e-6`) and optimality (`1e-4`) limits are unchanged.
`--max-iterations N` provides a bounded investigation;
`--output PATH`, `--config PATH` and `--verbose` are also available.

All generated files are grouped under `results/` beside the runner. The
comparison writes to `results/solver_comparison/`:

- `solver_comparison.png`, `.svg`, `.pdf`: efficiency, normalized constraint
  violation and independently assessed first-order optimality.
- `comparison.json`: versions, YAML hash, grid and all solver summaries.
- `timing_comparison.txt`: the final ASCII comparison also printed by the runner.
- `coolprop/configuration.json`, `bicubic/configuration.json`: effective recipes.
- Each `<backend>/<solver>/` directory: `history.csv`, full `history.json`, `summary.json`
  and `backend_events.json`.
- `tables/`: the comparison's cached property tables.
- `validation/`: retained test evidence and diagnostic logs.

Histories preserve the initial point and returned final point. SLSQP records
accepted steps. Pounce records native observations whose coordinates come from
its last objective evaluation; these are not certified accepted steps.
Restoration observations without an original-model point remain in raw backend
events. No best-so-far envelope replaces a trajectory. Exact zeros are preserved
on symmetric logarithmic feasibility/optimality axes; unavailable values are
gaps. Stopped runs are labelled as stopped.

The same KKT assessment is applied independently at each recorded point. Those
postprocessing model calls are reported separately and excluded from solve
times and solve-call totals. The final ASCII table reports setup and complete
solve seconds, iteration/call counts, efficiency and convergence status for all
six combinations. Setup includes model construction, table loading/generation
and one complete initial-cycle evaluation to prepare native JIT kernels. Solve
time includes the solver's stopping and final KKT checks and waits for the
returned cycle; extra history diagnostics, file output and plotting are outside
it. These are single observed runs, not repeated benchmark medians. Repeat with
balanced ordering before drawing general speed conclusions. Recognized invalid thermodynamic
trials are translated to `DomainError` at the comparison's optimizer boundary,
so line searches can reject them. No clipping, extrapolation, alternate fluid
backend or fabricated objective is used.

## Local verification, 16 November 2026

The initial bicubic-only 128x128 comparison converged for all three algorithms
(retained under `validation/bicubic_only_before_backend_comparison/`):

| Algorithm | Reported iterations | Final cycle efficiency | Maximum normalized violation |
| --- | ---: | ---: | ---: |
| Controlled SLSQP | 29 | 41.38198% | $7.44\times10^{-9}$ |
| Pounce interior point | 70 | 41.14839% | 0 |
| Pounce SQP | 21 | 41.38198% | $5.90\times10^{-10}$ |

These are native iteration counts, with the different observation semantics
described above. Both Pounce runs stopped after the common KKT assessment;
the interior-point endpoint has first-order optimality $9.12\times10^{-5}$,
versus approximately $1.38\times10^{-6}$ for the other two. Meeting the same
first-order limits does not guarantee identical endpoints in a nonlinear cycle.
The YAML hash was unchanged throughout. The original failed SQP run with default
QP tolerances is retained locally as `comparison_default_qp.json`.

The updated dependency environment passed **540 tests with four unchanged strict
expected failures**. A separate two-iteration headless manager run using the
128x128 cache exercised normal plotting, configuration/report callbacks, Excel
and final result saving. It does not constitute a full 200x200 interactive solve.
Generated evidence stays in the ignored results directories; no physical golden
values or solver acceptance tolerances were regenerated to make tests pass.

## CoolProp/bicubic timing comparison, 16 November 2026

The expanded runner completed one solve for each of the six combinations. All
converged; the YAML was unchanged and initial designs matched exactly. Each
solver used the same iteration/model-call counts with both backends. Auxiliary
fluids stayed on CoolProp, and the 128x128 bicubic table was loaded from cache.

| Algorithm | CoolProp solve [s] | Bicubic solve [s] | Observed solve-time reduction |
| --- | ---: | ---: | ---: |
| Controlled SLSQP | 9.769 | 6.756 | 30.84% |
| Pounce interior point | 24.307 | 17.877 | 26.45% |
| Pounce SQP | 8.642 | 6.197 | 28.30% |

These are single observed runs with the timing scope described above, not
repeated medians or a general speed guarantee. CoolProp/bicubic efficiencies
agree closely: approximately 41.38199%/41.38198% for SLSQP and SQP, and
41.14815%/41.14839% for interior point. Setup times are retained separately in
the [printed ASCII table](results/solver_comparison/timing_comparison.txt).
The seven focused runner tests pass; the earlier full-suite record remains
the full-suite verification for the dependency change.
