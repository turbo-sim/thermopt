# Rankine wet-to-dry study, v4: finer heat exchangers

This version reoptimizes simple n-butane, simple cyclopentane and
recuperated Novec 649 cycles with **100 elements in the evaporator and
condenser**. Novec's recuperator retains 25 elements. Each case starts
from its best matching saved coarse-grid solution and receives **one
optimization only**, without perturbation, multistart or recovery runs.
Results and thermodynamic interpretation are in [report.md](report.md).

## Conditions

| Parameter | Value |
|---|---|
| Geothermal source water | 140 °C, 10 bar |
| Minimum source outlet | 140, 135, …, 20 °C |
| Cooling water | 20 °C inlet, 24.99–25.00 °C outlet |
| Cycle sizing | 1 MW net cycle output, before external-water pumps |
| Superheated baseline | 2 K inlet superheat, dry outlet, expander efficiency 90% |
| Partial evaporation | Optimized inlet state, expander efficiencies 90%, 80%, 70%, 60% |
| Evaporator / condenser | 100 sampled locations each, minimum sampled approach 5 K |
| Novec recuperator | 25 sampled locations, minimum approach 5 K |
| Pressure losses | 1% per exchanger side |
| Pumps | 70% efficiency; working-fluid inlet subcooling at least 1 K |
| Solver | Native SLSQP, tolerance `1e-6`, maximum 250 iterations |

There are 120 finite-heat cases per fluid. The 140 °C outlet-limit
endpoint permits no source cooling, so its five combinations per fluid
are recorded as unavailable. The six-case comparison selects the
70 °C / 90%-expander-efficiency results from the sweeps.

## Initialization and single-run behavior

`prepare_coarse_seeds.py` compares accepted selected results from v3 and
v2 for the same case. It requires matching configuration and bounds,
reevaluates the highest-efficiency saved candidate on its original grid,
and snapshots it in `seed_results/sensitivity_<fluid>/`. This selection
does not run an optimizer. Each manifest records the source version,
coarse efficiency and configuration hash.

The new case uses that exact design-variable vector and configuration,
changing only `heater.num_elements` and `cooler.num_elements` to 100.
The shared YAML supplies the solver options. No random perturbation is
applied and no neighboring fine-grid result supplies the initial guess.

Each finite-heat case has one `attempt_started.json` marker and one row
in `attempts.csv`. A completed checkpoint is reused on subsequent script
runs, including an unsuccessful solve. An interrupted attempt without a
completed checkpoint stops execution rather than silently retrying.
Failed runs are retained and masked in figures; no replacement seed or
second optimization is attempted. The comparison runner never performs
duplicate optimizations.

## Scripts and running

All configuration is in globals; there are no CLI arguments.

For a fresh study, set `OPTIMIZE = True` in the sensitivity runner, then
run these commands from the repository root:

```console
python examples/rankine_wet2dry_v4/run_sensitivity_efficiency_and_exploitation.py
python examples/rankine_wet2dry_v4/run_fluid_comparison.py
python examples/rankine_wet2dry_v4/validate_saved_results.py
python examples/rankine_wet2dry_v4/build_report.py
```

The delivered runner defaults to `OPTIMIZE = False`, which regenerates
plots from saved results. `run_family(family)` can run independent fluid
studies in separate Python processes.

- `prepare_coarse_seeds.py`: validates and snapshots the previously
  selected coarse-grid seeds. Called automatically before optimization.
- `run_sensitivity_efficiency_and_exploitation.py`: performs one solve
  per case or regenerates sensitivity plots from saved states.
- `run_fluid_comparison.py`: exports six matching saved solutions with
  native ThermoOpt plots, workbooks and solver records, and one 2 × 2
  comparison figure per fluid (superheated/partial-evaporation columns;
  T–s/T–Q rows).
- `validate_saved_results.py`: checks balances, constraints, exact seed
  transfer and one optimization call per case. It also evaluates nine
  selected states per fluid and their coarse-grid seeds at 501 exchanger
  locations, without optimizing.
- `build_report.py`: rebuilds the numbered Markdown report from saved
  comparison, sensitivity and validation tables.

## Figures and efficiency definition

Efficiency contours and baseline curves use net system power after all
pumps divided by the source heat available between **140 and 20 °C**:

```text
net system power /
  [source mass flow × (h_water(140 °C, 10 bar) − h_water(20 °C, 10 bar))]
```

The denominator uses the same water enthalpy reference in every case.
ThermoOpt's native system efficiency is retained in the tables and remains
the optimization objective; at a fixed reinjection limit the objectives
differ by a positive constant. Relative gains compare with the same
fluid's 90%-efficient superheated baseline at that limit.

The native T–s/T–Q overlays show every other temperature with thin lines
and magma restricted to 0.25–0.75. Efficiency, exploitation and gain maps
use labeled `RdYlGn` contours. Quality is shown in percent using light blue
for 100% vapor and darker blue for more liquid; superheated vapor is
displayed as 100%.

`results/` contains the fine-grid designs, attempt records, tables and
PNG/PDF figures. `seed_results/` preserves the selected coarse-grid
inputs and provenance. Keep both directories when sharing or archiving
the study; their generated contents are ignored by Git.
