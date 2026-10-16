# Propane and butane comparison at a 140 degree source

This example compares simple propane cycles constrained to dry expansion with
simple n-butane cycles whose expander inlet state is free. Neither cycle has a
recuperator. Both fluids may optimize to either side of their critical pressure.
The results and discussion are in [report.md](report.md).

The conditions follow `examples/rankine_wet2dry`: water source at 140°C and 10 bar,
20°C cooling water warmed to approximately 25°C, 5 K minimum exchanger approaches
at 100 locations, 1% exchanger pressure losses, at least 1 K pump-inlet subcooling,
and 1 MW net cycle output before the external-water pumps.

The sweep varies allowed source exploitation from 0% to 100% in steps of 5
percentage points and working-fluid pump efficiency from 20% to 100% in steps of
10 percentage points. Exploitation sets a **minimum allowed source outlet**:
`Tmin = 140 − 1.2 × exploitation_pct` in °C. The actual outlet remains optimized.
The 0% endpoint has no available heat and is saved as unavailable; the first
finite endpoint is 5%, or a 134°C minimum outlet. The source and cooling-water
pumps retain 70% efficiency. Propane's expander efficiency is 90%, while butane's
is 90%, 80%, 70%, or 60%. This gives 900 finite cases and 45 zero-heat records.

The six variables are source and cooling-water outlet temperatures, pump inlet
pressure and enthalpy, and expander inlet pressure and enthalpy. The high-pressure
bound is three times each fluid's critical pressure; it spans subcritical and
supercritical conditions. Butane has no imposed quality constraint. Propane has
nonnegative outlet superheat plus 17 vapor-enthalpy path constraints; accepted
results also pass a 201-point dry-path check. The path assumes cumulative constant
isentropic efficiency: `h(p) = hin − ηt × [hin − h(p,sin)]`.

Run from the repository root with an environment containing the project dependencies:

```console
python examples/comparison_transcritical_two_phase/run_sensitivity.py
python examples/comparison_transcritical_two_phase/build_report.py
```

Configuration is in the globals of `run_sensitivity.py` and `case_template.yaml`.
The default performs the full grid. Set `OPTIMIZE = False` to regenerate figures
from `results/sweep_results.csv`, or run `plot_results.py` directly. `WORKERS`
controls independent fluid/expander-efficiency processes; use `WORKERS = 1` if
the execution environment disallows multiprocessing. Independent sweeps may also
be launched by calling `run_sweep((family, turbine_efficiency))`.

Each finite case gets two starts, with a third only for recovery or disagreement
above 0.2%. The best converged feasible result is saved. Per-case checkpoints
resume completed work and reject changed model inputs. Use a new `RESULTS`
directory when changing the formulation. Attempt logs, initial guesses, final
vectors, and selected native configurations are kept for review.

`fast_evaluator.py` calculates the same simple-cycle equations through low-level
CoolProp during optimization. Every completed solver candidate is reevaluated
with native ThermoOpt before acceptance; the native model supplies reported
efficiencies and plots. Saturation departure fields are recomputed from returned
states to avoid a wet-state derivative side effect in some jaxprop versions.
The shared library and the earlier example are unchanged.

All comparison plots use net system power after all pumps divided by the source
heat available between **140°C and 20°C** at 10 bar. This fixed denominator makes
efficiencies comparable across the exploitation axis. Native system efficiency
with the varying source reference remains the objective and is retained in the
CSV. At any given grid point the objectives are equivalent. Relative gain is
`100 × (butane efficiency / propane efficiency − 1)`; positive values favor butane.

`results/figures/` contains PNG and PDF versions of:

- Two 2×2 T–s/T–Q overlays, one at 5% and one at 100% allowed exploitation,
  with propane and butane at 90% expander efficiency and colors for all nine
  working-fluid pump efficiencies.
- Two system-efficiency contour plots, propane and butane at 90%.
- Four relative-gain maps, one for each butane expander efficiency.
- Paired pressure, back-work ratio, and actual-exploitation maps.

T–Q retains the native cumulative heat-duty coordinate and 1 MW cycle sizing.
T–s external-stream traces use associated working-fluid entropy coordinates.
The overlays use magma restricted to 0.25–0.75, a labeled pump-efficiency colorbar,
and line styles distinguishing working fluid, source, and cooling water.
Failed points and zero-heat cases are masked in the contour plots.
Gain maps mark propane's critical-pressure boundary and butane's saturated-vapor
inlet boundary. The report also evaluates the subset that directly compares
liquid or wet-admission butane against transcritical propane.

`build_report.py` generates the report and diagnostic tables from saved results;
it does not optimize. `check_results.py` performs a short consistency check of the
saved grid and a few native solutions. Generated `results/` contents are ignored
by Git; keep that directory with the scripts and report when sharing the study.
