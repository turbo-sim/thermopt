# Recuperated sCO2 cycle examples

Install the checkout with `poetry install`. From the repository root, compare
the SLSQP configurations with:

```bash
poetry run python examples/sCO2_recuperated/compare_slsqp.py
# Save the comparison without opening a window:
poetry run python examples/sCO2_recuperated/compare_slsqp.py --no-show
```

The script runs six independent optimizations from the same physical initial
guess in `case_sCO2_recuperated.yaml`:

| Run | Problem scale | Initial step fraction | Plot style |
| --- | --- | --- | --- |
| SciPy SLSQP | 1 | Native behavior | Black solid |
| SciPy SLSQP | YAML value, currently 10 | Native behavior | Black dashed |
| Custom SLSQP | 1 | 1.0 (full initial step) | Plasma |
| Custom SLSQP | 1 | 0.3 | Plasma |
| Custom SLSQP | 1 | 0.2 | Plasma |
| Custom SLSQP | 1 | 0.1 | Plasma |

The custom fractions multiply the native proposed first step. The cap then
grows by a factor of 1.5 toward the full step; it is not a permanent step cap.
Custom domain backtracking remains enabled. Each run has its own cycle problem,
fluid states, and solver, so one run cannot change another's initial guess or
scaling. No cycle-plotting or file-saving callbacks run during optimization.

The YAML supplies the native stopping tolerance, iteration limit, and numerical
derivative settings. Custom final acceptance uses feasibility tolerance `1e-6`
and optimality tolerance `1e-4`. `--max-iterations N` overrides the iteration
limit for every run; `--verbose` prints the per-iteration solver reports.

The three panels in one row have full axes boxes and show the minimized
objective $f=-\eta_{\mathrm{cycle}}$, maximum normalized constraint violation
(including bounds), and first-order optimality $\|\nabla_x L\|_\infty$.
SciPy observations
are native major-iteration callbacks; custom observations are accepted steps.
The initial and distinct final observations are retained. Histories are not
replaced by a best-so-far envelope. The violation panel uses a symmetric log
scale with a linear region below `1e-8`, retaining exact zeros. The optimality
panel uses the same scale. Dotted lines mark `1e-6` feasibility and `1e-4`
optimality. A stopped run retains its recorded history and is labeled as such.

Optimality is recomputed after each solve at every recorded point, using a
separate cycle problem, bounded numerical derivatives, and the shared KKT
assessment. It is measured in each run's scaled solver coordinates, so changing
the problem scale also changes the gradient's units. These additional model
calls are reported separately as `optimality_model_calls` in `summary.json`;
the original optimization call counts and trajectories are retained. Unavailable
optimality values appear as gaps. `history.csv` includes the optimality column.

Outputs are saved beside the script in `results_slsqp_comparison/`:

- `slsqp_comparison.png`, `.svg`, and `.pdf`: the three-panel comparison.
- One directory per run, containing `history.csv` and `summary.json`.

Use `--output PATH` to select another output directory or `--config PATH` to
read another recuperated-cycle configuration. Paths work from any working
directory after installing ThermOpt. The existing `run_optimization.py`
remains the single-run example with interactive cycle plotting.
