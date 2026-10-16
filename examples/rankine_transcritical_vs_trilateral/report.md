# Propane and butane cycle comparison at a 140 degree source

Simple propane cycles with dry expansion and simple n-butane cycles with a free expander inlet state are optimized to compare the benefit of source matching against working-fluid pump losses and prescribed expander efficiency. Positive relative gains favor butane at the same allowed source cooling and pump efficiency.

With maximum allowed exploitation and a 70%-efficient pump, butane at 90% expander efficiency reaches **10.29%** fixed-reference system efficiency versus **8.84%** for propane, a **+16.37%** relative gain. Butane operates at **23.34 bar** with a **13.05%** back-work ratio, versus **50.12 bar** and **19.93%** for propane.

For propane, reducing pump efficiency from 100% to 20% at maximum allowed exploitation reduces optimized system efficiency by **45.73%** at 90% expander efficiency.

For butane, reducing pump efficiency from 100% to 20% at maximum allowed exploitation reduces optimized system efficiency by **31.01%** at 90% expander efficiency.

## Boundary conditions and comparison basis

The source is water at 140°C and 10 bar. Cooling water enters at 20°C and exits at 24.99–25.00°C. Both cycles have no recuperator, 1% pressure loss on each exchanger side, at least 5 K approach at 100 uniform exchanger locations, and at least 1 K pump-inlet subcooling. Every design is sized to 1 MW after the working-fluid pump, before the external-water pumps. Those external pumps retain 70% efficiency.

Allowed exploitation is `(140 − minimum source outlet temperature) / 120`, sampled in 5 percentage-point steps. The source outlet is optimized above this minimum, so 100% permission does not require cooling to 20°C. Working-fluid pump efficiency is sampled from 20% to 100% in 10 percentage-point steps. Propane's expander efficiency is 90%; butane's is 90%, 80%, 70%, or 60%.

The free variables are source outlet temperature, cooling-water outlet temperature, working-fluid pump inlet pressure and enthalpy, and expander inlet pressure and enthalpy. High-pressure bounds extend to three times the critical pressure for both fluids and permit subcritical or supercritical heating. No wet inlet is imposed on butane. Propane has nonnegative exit superheat and 17 additional vapor-enthalpy constraints along the saturation-capable part of a cumulative constant-isentropic-efficiency expansion path. Accepted propane designs also pass a 201-point check of that path.

The figures use `η_fixed = net system power / [source mass flow × Δh_water(140→20°C at 10 bar)]`, with **Δh = 504.723 kJ/kg**. Net system power subtracts all three pumps. This measures electricity per unit source flow on a common heat reference. The native system efficiency, divided by heat available down to the allowed source outlet, remains the optimization objective and is exported separately. At a given grid point the denominators differ by a positive constant, so both objectives select the same optimum.

Relative change is `100 × (η_fixed,butane / η_fixed,propane − 1)`. Since the denominator is common, this is also the net power change at fixed source mass flow. Percentage-point differences are exported separately. At 0% exploitation no heat can be released; these cases are recorded as unavailable and receive no optimization.

## Optimization and numerical checks

Each finite point receives a continuation or high-pressure start and an independent vapor-oriented start. A third, perturbed start is used when neither succeeds or accepted starts differ by more than 0.2%. The maximum is three starts per point. The best feasible converged result is retained; all attempts are logged. These limited restarts reduce sensitivity to initialization but do not establish a global optimum.

The optimization evaluator computes the same simple-cycle heat balances and uniform exchanger constraints using low-level CoolProp. Native ThermoOpt reevaluates every completed candidate before acceptance and supplies all reported efficiencies and cycle diagrams. Saturation departures used in constraints are recomputed from returned pressure, enthalpy, and temperature because some jaxprop versions perturb their internal state during wet-state derivative calculations. No repository core code is changed.

| Fluid | Expander efficiency | Accepted | Requested including zero heat | Optimization calls |
| --- | --- | --- | --- | --- |
| butane | 60% | 180 | 189 | 502 |
| butane | 70% | 180 | 189 | 504 |
| butane | 80% | 180 | 189 | 497 |
| butane | 90% | 180 | 189 | 496 |
| propane | 90% | 180 | 189 | 526 |

The accepted solutions have maximum normalized constraint violation **7.76e-07** and maximum energy-balance residual **4.33e-06 W**. Failed cases are masked in every plot and comparison.

## Cycle pressure and phase selected by optimization

| Fluid | Expander efficiency | Pressure regime | Inlet phase | Outlet phase | Cases |
| --- | --- | --- | --- | --- | --- |
| butane | 60% | subcritical | two_phase | two_phase | 75 |
| butane | 60% | subcritical | two_phase | vapor | 42 |
| butane | 60% | subcritical | vapor | vapor | 63 |
| butane | 70% | subcritical | two_phase | two_phase | 80 |
| butane | 70% | subcritical | two_phase | vapor | 42 |
| butane | 70% | subcritical | vapor | vapor | 58 |
| butane | 80% | subcritical | liquid | two_phase | 3 |
| butane | 80% | subcritical | two_phase | two_phase | 78 |
| butane | 80% | subcritical | two_phase | vapor | 52 |
| butane | 80% | subcritical | vapor | vapor | 47 |
| butane | 90% | subcritical | liquid | two_phase | 6 |
| butane | 90% | subcritical | two_phase | two_phase | 78 |
| butane | 90% | subcritical | two_phase | vapor | 47 |
| butane | 90% | subcritical | vapor | vapor | 49 |
| propane | 90% | subcritical | vapor | vapor | 29 |
| propane | 90% | transcritical | supercritical_pressure | vapor | 151 |

At maximum allowed exploitation, propane selects subcritical pressure at pump efficiencies 20%, 30% and transcritical pressure at 40%, 50%, 60%, 70%, 80%, 90%, 100%. This transition is an optimization result within the permitted pressure bounds. Poor pump efficiency can change the preferred pressure regime as well as reducing efficiency.

Liquid-admission butane solutions approach the trilateral configuration. The assumed liquid or wet expander performance is still an input: the optimizer establishes the thermodynamic state and its required work, rather than demonstrating that a practical expander attains the prescribed efficiency.

A subcritical propane result is part of the permitted comparison. Likewise a vapor-inlet butane result is a conventional Rankine-like cycle within the flexible formulation. The assumed butane expander efficiency applies to its entire sweep, including vapor-inlet points. These points should not be described as demonstrated advantages of two-phase admission.

## Cycle diagrams at low and high allowed exploitation

![Cycles at minimum finite exploitation](results/figures/cycles_exploitation_005.png)

![Cycles at maximum allowed exploitation](results/figures/cycles_exploitation_100.png)

Both figures show the 90%-efficient expander slices and all accepted pump efficiencies. Magma colors encode pump efficiency; solid, dashed, and dotted lines indicate working fluid, source, and cooling water. The expansion curves follow the assumed cumulative-efficiency path. External-stream traces in T–s use the associated working-fluid entropy coordinate, following the native plot convention. T–Q uses the native cumulative duty coordinate and the common 1 MW net-cycle sizing. Heat duties therefore change with the flow needed by each design.

## Net system efficiency

![Propane system efficiency](results/figures/system_efficiency_propane_turb90.png)

![Butane system efficiency](results/figures/system_efficiency_butane_turb90.png)

The two efficiency maps share their color levels. Actual exploitation can plateau below the allowed limit because additional cooling becomes less valuable than maintaining a favorable heat-addition profile and expansion work.

## Relative change as butane expander efficiency decreases

![Relative gain at butane expander 90%](results/figures/relative_gain_butane_turb90.png)

At **90% butane expander efficiency**, butane wins at **158 of 180** accepted paired conditions. Relative changes range from **-7.87% to +64.85%**. Of 180 accepted butane designs, 125 have two-phase admission.

![Relative gain at butane expander 80%](results/figures/relative_gain_butane_turb80.png)

At **80% butane expander efficiency**, butane wins at **93 of 180** accepted paired conditions. Relative changes range from **-18.89% to +41.15%**. Of 180 accepted butane designs, 130 have two-phase admission.

![Relative gain at butane expander 70%](results/figures/relative_gain_butane_turb70.png)

At **70% butane expander efficiency**, butane wins at **20 of 180** accepted paired conditions. Relative changes range from **-30.00% to +17.69%**. Of 180 accepted butane designs, 122 have two-phase admission.

![Relative gain at butane expander 60%](results/figures/relative_gain_butane_turb60.png)

At **60% butane expander efficiency**, butane wins at **0 of 180** accepted paired conditions. Relative changes range from **-40.96% to -5.26%**. Of 180 accepted butane designs, 117 have two-phase admission.

At maximum allowed exploitation, relative efficiency changes are:

| Pump efficiency | Butane expander 90% | Butane expander 80% | Butane expander 70% | Butane expander 60% |
| --- | --- | --- | --- | --- |
| 20% | +44.12% | +21.78% | +0.00% | -21.07% |
| 30% | +33.20% | +14.12% | -4.65% | -23.00% |
| 40% | +25.68% | +8.49% | -8.50% | -25.21% |
| 50% | +21.10% | +5.04% | -10.88% | -26.60% |
| 60% | +18.09% | +2.76% | -12.44% | -27.50% |
| 70% | +16.37% | +1.50% | -13.23% | -27.86% |
| 80% | +15.02% | +0.53% | -13.86% | -28.16% |
| 90% | +13.95% | -0.13% | -14.28% | -28.31% |
| 100% | +13.38% | -0.66% | -14.56% | -28.39% |

The zero contour gives an interpolated break-even boundary. `sampled_break_even.csv` records the first sampled exploitation with positive gain for each pump and butane-expander efficiency. The 5% exploitation and 10% pump grids limit the precision of any boundary.

The blue dashed line marks propane's critical-pressure boundary; the black dotted line marks the butane saturated-vapor inlet boundary. Both are interpolated from optimized grid values. Gains elsewhere can compare subcritical propane or vapor-admission butane. The subset that directly compares liquid or wet-admission butane with transcritical propane is:

| Butane expander efficiency | Paired conditions | Butane wins | Minimum gain [%] | Maximum gain [%] |
| --- | --- | --- | --- | --- |
| 90% | 110 | 88 | -7.87 | +31.13 |
| 80% | 110 | 33 | -18.89 | +13.83 |
| 70% | 102 | 0 | -30.00 | -4.97 |
| 60% | 100 | 0 | -40.96 | -21.17 |

In this direct architecture comparison, the 70% and 60% butane expanders do not beat the 90% propane expander at any sampled condition. Positive gains in their full fluid maps must be read with the phase and pressure boundaries. In particular, the 70% butane expander is essentially tied at 20% pump efficiency and maximum allowed exploitation, where propane selects a subcritical cycle; that near-zero difference is smaller than the observed local-solution irregularities.

Linear interpolation between neighboring prescribed butane expander efficiencies gives the following approximate break-even values at maximum allowed exploitation. Each fluid is reoptimized at the sampled efficiencies, so these estimates describe a boundary between optimized designs rather than a fixed-design expander calculation.

| Pump efficiency | Approximate butane expander break even [%] |
| --- | --- |
| 20% | 70.00 |
| 30% | 72.48 |
| 40% | 75.00 |
| 50% | 76.83 |
| 60% | 78.18 |
| 70% | 78.98 |
| 80% | 79.63 |
| 90% | 80.09 |
| 100% | 80.47 |

## Pressure levels and back work explain the tradeoff

| Fluid | Pump efficiency | Efficiency [%] | High pressure [bar] | Back work [%] | Source outlet [°C] | Inlet state | Outlet state |
| --- | --- | --- | --- | --- | --- | --- | --- |
| butane | 20% | 7.43 | 14.73 | 29.11 | 38.12 | 0.40 vapor fraction | 0.72 vapor fraction |
| propane | 20% | 5.16 | 30.13 | 34.34 | 59.60 | vapor | vapor |
| butane | 50% | 9.71 | 21.32 | 16.58 | 35.81 | 0.18 vapor fraction | 0.70 vapor fraction |
| propane | 50% | 8.01 | 46.93 | 26.58 | 50.47 | supercritical pressure | vapor |
| butane | 70% | 10.29 | 23.34 | 13.05 | 34.89 | 0.10 vapor fraction | 0.69 vapor fraction |
| propane | 70% | 8.84 | 50.12 | 19.93 | 51.89 | supercritical pressure | vapor |
| butane | 100% | 10.77 | 26.87 | 10.54 | 34.34 | liquid | 0.69 vapor fraction |
| propane | 100% | 9.50 | 50.43 | 14.02 | 51.45 | supercritical pressure | vapor |

![Optimized pressures](results/figures/optimized_pressure.png)

![Backwork ratio](results/figures/backwork_ratio.png)

![Actual exploitation](results/figures/actual_exploitation.png)

The comparison includes changes in both heat recovery and conversion efficiency. Back-work ratio is the working-fluid pump power divided by expander power. Absolute pressure, pressure ratio, inlet density, and specific expansion work all contribute; lower absolute pressure alone does not guarantee a smaller back-work ratio. The optimizer can respond to poor pump efficiency by changing pressure and inlet state, so a fixed-design pump-loss estimate does not predict the optimized transition by itself.

At that 70%-pump, 90%-expander point, butane cools the source to **34.89°C** while propane leaves it at **51.89°C**. Recovered heat per unit source flow changes by **+19.13%**, and net power per unit absorbed heat changes by **-2.32%**. These factors multiply to give the **+16.37%** system-power gain. Thus the benefit combines source recovery and conversion, and can persist even when conversion efficiency by itself is slightly lower.

## Local optimization and interpretation limits

There are **16** neighboring exploitation pairs where fixed-reference efficiency falls by more than 0.01% despite allowing more source cooling; the largest decrease is **0.750%**. These are local-solution irregularities because relaxing the outlet-temperature limit retains every earlier design's feasibility. They are preserved in `local_solution_irregularities.csv` and should not be interpreted as a thermodynamic penalty for allowing deeper cooling.

This is a thermodynamic comparison of two fluids and prescribed component efficiencies. It does not isolate fluid choice from cycle architecture or predict the efficiency, cost, leakage, volume ratio, or lifetime of a real wet expander. The triangular heat-addition profile can improve recovery, but the required expander performance must be achievable at the optimized inlet and outlet qualities. Propane's dry-path requirement is tied to the stated expansion-path model. Exchanger approaches are sampled constraints, so small gains near zero remain sensitive to discretization and local optimization.
