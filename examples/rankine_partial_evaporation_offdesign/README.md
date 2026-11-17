# Partial-evaporation ORC off design (Novec 649)

Design, sizing and off-design operation on fixed hardware of a simple partial-evaporation organic Rankine cycle. A pressurized hot-water source (10 kg/s, 140 °C) heats Novec 649, and cooling water at 20 °C is the sink. All values are example assumptions, not plant data.

Run from this folder. In IPython:

```python
%run run_offdesign.py          # design and sizing, about 1.5 min
expander = f(0.9, 5.0, 4.0)    # one off-design point, 20-90 s
expander["state_in"].p, expander["inlet_quality"], expander["state_out"].p
```

`f(mass_flow_ratio, delta_T_source, delta_T_sink)` changes the source flow (ratio to design) and the source and sink inlet temperatures (K relative to design), solves for maximum net power from the design point, and returns the expander component. `python run_offdesign.py` solves two example points.

Run time of `f` (Windows laptop, one core): about 20 s at design conditions (`f()`), 40 s for `f(0.9, 5, 4)` and 90 s for `f(0.7, -15, 8)`. It grows with the distance from design, because every call starts from the design point.

## What the script does

1. **Design** with the existing `simple` topology ([`case_novec649_design.yaml`](case_novec649_design.yaml)). The design constraints are:
   - 5 K pinch in both exchangers;
   - 1 K pump-inlet subcooling;
   - two-phase turbine inlet, $0.05\le Q_3\le0.95$.

   The objective is the system efficiency, which maximizes net power per kg of source water.
2. **Sizing** with `cycle_power_simple_offdesign.size_hardware`. The intensive design is scaled to 10 kg/s of source water, and the design point is re-evaluated to fix:
   - evaporator and condenser $UA$;
   - the effective throat area of the turbine nozzle (choked at design);
   - the turbine design spouting velocity (the turbine efficiency is $\eta_d\,x(2-x)$ with the velocity ratio relative to design $x=c_{0,d}/c_0$);
   - a variable-speed pump map anchored at 3000 rpm.
3. **Variables** for off design: the design point is the start, with bounds around it.
4. **Off-design points** with `f`, solved directly from the design point with no continuation. The optimizer chooses pump speed and sink flow, with sink flow capped at its design value, for maximum net power. To fix a control instead, put it in `fixed_parameters.controls` and remove it from the variables.

The slow test (`THERMOPT_RUN_SLOW=1 pytest tests/test_offdesign.py`) also checks consistency: at design conditions with the design controls fixed, the off-design model reproduces the design net power within $10^{-6}$.

## Design and off design in ThermoOpt

The off-design cycle is a normal ThermoOpt topology, `power_simple_offdesign`, with the usual `evaluate_cycle(variables, parameters, constraints, objective_function)`. Compared with the design file, mainly the YAML changes ([`case_novec649_offdesign.yaml`](case_novec649_offdesign.yaml)):

| | Design (`simple`) | Off design (`power_simple_offdesign`) |
| --- | --- | --- |
| Variables | pump- and turbine-inlet $p,h$; source and sink outlet $T$ | pump-inlet $p$; turbine-inlet $p,h$; source outlet $T$; pump speed and sink flow (if optimized) |
| Parameters | efficiencies, pressure drops, boundary conditions | sized `hardware`, boundary conditions, fixed `controls` |
| Equalities | subcooling | $\ln(UA_{req}/UA_{av})$ for both exchangers, pump head, source energy balance |
| Inequalities | pinch, inlet quality | inlet quality, pump flow coefficient, turbine velocity ratio relative to design, evaporator specific duty |
| Objective | system efficiency | net power relative to design |

A control listed under `design_variables` is optimized; one listed under `fixed_parameters.controls` is fixed.

## Robust formulation

A start far from the solution should still evaluate, so the solver can move back instead of failing at the first step.
- **Bounded source outlet.** The source outlet temperature is an unknown bounded by the sink and source inlet temperatures, and its energy balance is an equality. Every evaluated water state is valid, even when the start puts far more duty on the source than it can deliver.
- **Crossing extension.** Exchanger node approaches below 0.2 K are replaced by $\Delta T^*=\delta\,e^{(\Delta T-\delta)/\delta}$ with $\delta=0.2$ K, which matches value and slope. The UA match is written as a log ratio, so crossed profiles give a large but continuous residual (about 5 per kelvin of crossing), up to a crossing of about 140 K where UA overflows. Above 0.2 K nothing changes, so a solved point is exact when the reported `min_approach` of both exchangers is at least 0.2 K.
- **Choked or unchoked nozzle.** The turbine is one equivalent nozzle: $\dot m=C_dA_t\,G(\max(p_4,p^*))$, the equilibrium choked flux below the critical pressure $p^*$ and the isentropic flux at the back pressure above it. Value and slope are continuous at $p^*$.

In the case study behind this example (two fluids, sweeps, fold analysis), an earlier version of this formulation (crossing threshold 0.1 K) solved all 105 combinations of source flow (70–110 %), source temperature (±15 K) and sink temperature (±8 K) directly from the design point.

## Expected results

Measured on Windows, Python 3.13, jaxprop 0.5.8, pysolver_view 0.7.4:

| Run | $W_{net}$ [kW] | $p_3$ [bar] | $Q_3$ | $p_4$ [bar] | Pump speed [rad/s] | Sink flow [kg/s] |
| --- | --- | --- | --- | --- | --- | --- |
| Design | 559.48 | 8.841 | 0.050 | 0.452 | 314.16 | 315.52 |
| `f()`: design conditions, optimal controls | 559.52 | 8.845 | 0.056 | 0.451 | 313.95 | 315.52 (at cap) |
| `f(0.9, 5, 4)`: 9 kg/s, 145 °C, 24 °C | 509.28 | 8.737 | 0.164 | 0.516 | 307.77 | 315.52 (at cap) |

`f()` differs slightly from the design point because the design limits (5 K pinch, $Q_3\ge0.05$) do not apply to fixed hardware; the net-power optimum is flat, so $Q_3$ moves for a 0.007 % power gain. With the design controls fixed, the design point is reproduced (slow test). The optimizer keeps the sink flow at its cap, the design value: more cooling water always lowers the condensing pressure, and the sink-pump power grows only linearly with flow because the cold-side pressure drop is a fixed fraction.

## Files

- `components_offdesign.py` (in `thermopt/components/`): pump map, turbine map, nozzle flow (with `find_choked_state` adapted from `flashing_nozzle_design`), UA integration.
- `cycle_power_simple_offdesign.py` (in `thermopt/cycles/`): cycle and sizing.
- Tests: `tests/test_offdesign.py`. The end-to-end run of this example, with the consistency check, is included with `THERMOPT_RUN_SLOW=1`.
