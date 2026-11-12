# Partial evaporation in geothermal Rankine cycles

Updated EMPOWER study: n-butane, n-pentane, cyclopentane and Novec 649. The boundary conditions follow Nishith's 20 October 2026 request and Cristina's 19 October clarification. The standard-cycle reference and the 70 °C comparison use Cristina's 88.6% turbine isentropic efficiency. Partial-evaporation sensitivities retain 90%, 80%, 70% and 60% efficiency.

At equal 88.6% turbine efficiency and a 70 °C source-outlet limit, partial evaporation changes net power per geothermal flow by n-Butane +4.50%, n-Pentane +7.50%, Cyclopentane +13.69%, Novec649 +12.07% relative to each fluid's own standard cycle. The highest output among these eight saved designs is 42.181 kW per kg/s, for novec partial evaporation. Deeper permitted source cooling can increase the advantage by improving heat recovery as well as conversion. These boundary-condition changes include the cooling-water inlet and turbine reference, so differences from the previous report cannot be attributed to pinch alone.

## 1. Boundary conditions and interpretation

| Parameter | Setting |
|---|---|
| Source water | 140 °C; 10 bar |
| Cooling water | 15 °C inlet; exactly 25 °C return |
| External heater / condenser minimum approach | 4.5 K; 100 locations each |
| Novec internal recuperator | 5 K minimum approach; 25 locations; optimized effectiveness |
| Source minimum outlet | 70 °C comparison; 140 to 20 °C in 5 K sensitivity steps |
| Standard-cycle turbine efficiency | 88.6% throughout |
| Partial-evaporation turbine efficiency | 90%, 80%, 70%, 60%; additional 88.6% cases at 70 °C |
| Standard inlet / outlet | 2 K inlet superheat; dry expander outlet |
| Partial evaporation | Optimized inlet pressure and enthalpy; wet inlet and dry outlet not imposed |
| Architectures | Simple n-butane, n-pentane and cyclopentane; recuperated Novec 649 |
| Pumps / pressure drops | 70% isentropic pump efficiencies; 1% per active exchanger side |
| Working-fluid pump inlet | At least 1 K subcooling |
| Sizing | 1 MW expander output minus working-fluid pump; external pumps deducted subsequently |
| Optimization | SLSQP, tolerance 1e-6, maximum 250 iterations per attempt |

The 4.5 K instruction applies to the external hot exchangers and condensation side for both cycle types. The existing 5 K internal-recuperator assumption is retained because the exchange does not explicitly revise it. The original 2 K superheat, pump efficiencies, pressure losses and architecture assumptions are also retained. The spreadsheets do not fully specify these quantities, generator losses or auxiliary-power accounting, so this study aligns the stated boundary conditions rather than claiming an exact reproduction of EXERGY electrical outputs.

A 70 °C source limit is a lower bound, not a fixed outlet. The optimizer can leave the source warmer. Permission to cool the source to 20 °C likewise does not mean the optimum reaches 20 °C.

### Optimization variables

| Variable | Role |
|---|---|
| Source outlet temperature | Optimized above the case's minimum outlet |
| Cooling-water outlet temperature | Retained as a variable and constrained to exactly 25 °C |
| Working-fluid pump inlet pressure and enthalpy | Low-pressure state and at least 1 K subcooling |
| Expander inlet pressure and enthalpy | High-pressure state and admission phase |
| Recuperator effectiveness | Optimized between 0 and 1 for Novec only |

Mass flows follow from the optimized states, heat balances and 1 MW cycle-sizing target. State bounds are broad and property-based. The standard-cycle superheat and dry-outlet constraints restrict the same inlet-state variables that remain flexible in the partial-evaporation formulation.

## 2. Normalization and computation

Fixed-reference efficiency is net system power divided by source mass flow and the water enthalpy drop from 140 to 15 °C at 10 bar. This reference heat is 525.639 kJ/kg. It is common to all designs, and efficiency therefore tracks net power per unit geothermal flow directly. Nominal exploitation is 100 × (140 − minimum outlet) / 125; actual exploitation replaces the minimum by the achieved outlet. A 20 °C limit corresponds to 96% nominal exploitation on this reference.

ThermoOpt's native system-efficiency objective uses available source heat down to each case's minimum outlet. At a fixed limit, maximizing that objective is equivalent to maximizing fixed-reference efficiency or net power per unit source flow. Reported powers at 194.44 kg/s scale each design to the source flow in Cristina's workbook. They remain model net powers before any separately specified generator losses.

Every finite-heat case is optimized through the public thermopt.ThermodynamicCycleOptimization entry point using its native ThermodynamicCycleProblem, cycle equations and component models. No custom fitness function or alternate cycle evaluator is used. The previously converged designs supply initial guesses only. Each attempt preserves its complete case.yaml, initial and final configurations, convergence history, solver log and optimization report. Positive physical flows and work, constraint satisfaction and energy balance are checked before acceptance. Native source-code hashes are recorded so the implementation can be verified. A numerical solution is not proof of global optimality or real two-phase expander performance.

| Fluid | Accepted finite-heat cases | Failed finite-heat cases | Zero-heat unavailable |
|---|---|---|---|
| n-Butane | 121 | 0 | 5 |
| n-Pentane | 121 | 0 | 5 |
| Cyclopentane | 121 | 0 | 5 |
| Novec649 | 121 | 0 | 5 |

There are 484 finite-heat designs: 24 temperatures × five designs × four fluids, plus four partial-evaporation designs at 88.6% for the 70 °C comparison. The 20 requested combinations at a 140 °C minimum outlet cannot release geothermal heat and are unavailable. Each fluid's sensitivity gain uses its standard cycle at 88.6% and the matching source-outlet limit.

## 3. Eight-case comparison at 70 °C and 88.6% turbine efficiency

| Case | Inlet vapor % | Outlet vapor % | Source out °C | Net kW/(kg/s) | Gain vs standard n-butane % | Gain vs same-fluid standard % |
|---|---|---|---|---|---|---|
| butane superheated | 100.00 | 100.00 | 70.00 | 36.780 | +0.00 | +0.00 |
| butane partial evaporation | 58.00 | 86.45 | 70.00 | 38.436 | +4.50 | +4.50 |
| pentane superheated | 100.00 | 100.00 | 70.00 | 36.128 | -1.77 | +0.00 |
| pentane partial evaporation | 23.11 | 71.39 | 70.00 | 38.838 | +5.60 | +7.50 |
| cyclopentane superheated | 100.00 | 100.00 | 70.01 | 34.436 | -6.37 | +0.00 |
| cyclopentane partial evaporation | 16.15 | 51.89 | 70.00 | 39.150 | +6.45 | +13.69 |
| novec superheated | 100.00 | 100.00 | 70.00 | 37.639 | +2.34 | +0.00 |
| novec partial evaporation | 61.15 | 100.00 | 70.00 | 42.181 | +14.69 | +12.07 |

The common cross-fluid reference is standard n-butane, matching Cristina's spreadsheet. The same-fluid column isolates the change in inlet-state formulation at equal turbine efficiency. Superheated vapor is displayed as 100% quality; raw quality extensions remain in the data. Partial evaporation may produce wet-to-wet, wet-to-dry or dry-to-dry endpoints.

The comparison exporter writes optimal_solution.xlsx with native cycle states and energy analysis, and design_variable_bounds.csv with the saved optimum, lower and upper bounds, fractional distance to the nearest bound and active-bound flags. Each case has its own directory under output/results/fluid_comparison/. The HTML comparison page offers the generated downloads when available. A bound is flagged when the optimum lies within 0.01% of its variable range.

![n-Butane: simple cycle: standard and partial evaporation; T-s above T-Q; turbine efficiency 88.6%.](output/results/fluid_comparison/comparison_grid_butane.png)

[Embeddable HTML](output/results/fluid_comparison/comparison_grid_butane.html) · [SVG](output/results/fluid_comparison/comparison_grid_butane.svg) · [PDF](output/results/fluid_comparison/comparison_grid_butane.pdf)

![n-Pentane: simple cycle: standard and partial evaporation; T-s above T-Q; turbine efficiency 88.6%.](output/results/fluid_comparison/comparison_grid_pentane.png)

[Embeddable HTML](output/results/fluid_comparison/comparison_grid_pentane.html) · [SVG](output/results/fluid_comparison/comparison_grid_pentane.svg) · [PDF](output/results/fluid_comparison/comparison_grid_pentane.pdf)

![Cyclopentane: simple cycle: standard and partial evaporation; T-s above T-Q; turbine efficiency 88.6%.](output/results/fluid_comparison/comparison_grid_cyclopentane.png)

[Embeddable HTML](output/results/fluid_comparison/comparison_grid_cyclopentane.html) · [SVG](output/results/fluid_comparison/comparison_grid_cyclopentane.svg) · [PDF](output/results/fluid_comparison/comparison_grid_cyclopentane.pdf)

![Novec 649: recuperated cycle: standard and partial evaporation; T-s above T-Q; turbine efficiency 88.6%.](output/results/fluid_comparison/comparison_grid_novec.png)

[Embeddable HTML](output/results/fluid_comparison/comparison_grid_novec.html) · [SVG](output/results/fluid_comparison/comparison_grid_novec.svg) · [PDF](output/results/fluid_comparison/comparison_grid_novec.pdf)

## 4. Comparison with Cristina's revised workbook

The original and revised workbooks report the same performance values. The revision adds condenser-inlet pressures and clarifies turbine-inlet pressure. The overlapping cases below are compared at 194.44 kg/s source flow. Pressures and inlet qualities in that workbook are optimized outputs, so they are validation targets rather than imposed inputs. The standard cyclopentane condenser-inlet pressure is blank in the workbook.

| Case | Workbook net kWe | Model net kW | Turbine p: workbook / model bar | Condenser p: workbook / model bar | Inlet vapor: workbook / model % |
|---|---|---|---|---|---|
| butane superheated | 6795.8 | 7151.4 | 13.55 / 13.46 | 2.77 / 2.72 | 100.0 / 100.0 |
| pentane superheated | 6396.0 | 7024.8 | 4.34 / 4.64 | 0.78 / 0.77 | 100.0 / 100.0 |
| pentane partial evaporation | 7061.0 | 7551.8 | 8.33 / 9.33 | 0.81 / 0.81 | 31.0 / 23.1 |
| cyclopentane superheated | 5917.0 | 6695.7 | 2.81 / 2.63 | unavailable / 0.50 | 100.0 / 100.0 |
| cyclopentane partial evaporation | 6988.0 | 7612.3 | 6.09 / 6.61 | 0.51 / 0.50 | 19.0 / 16.2 |

Source: summary results_rev1.xlsx, Summary_results!C5:I17; assumptions clarified in RE: EMPOWER Project - ORC Conference, 19–20 October 2026. Net electrical and model mechanical powers may have different generator and auxiliary-loss definitions. The workbook's comparison row has a stray bar unit label; its values are dimensionless relative gains.

## 5. n-Butane: simple cycle

At a 70 °C source limit, standard-cycle fixed-reference efficiency is 7.00%. Partial evaporation at 90% turbine efficiency reaches 7.44%, a 6.29% gain over the 88.6%-efficient standard reference. This sensitivity comparison includes the difference in turbine efficiency; the preceding eight-case comparison holds turbine efficiency equal.

At equal 88.6% turbine efficiency and a 70 °C source limit, standard admission is 13.46 bar and 95.66 °C. Partial-evaporation admission is 17.98 bar and 108.61 °C. Flexible admission changes pressure level and the temperature of heat addition along with the vapor fraction.

At the 20 °C limit, standard and 90%-efficient partial cycles produce 37.37 and 52.92 kW per kg/s of source flow. Their achieved source outlets are 63.06 and 34.49 °C. The partial-expansion inlet/outlet vapor fractions are 0.1071 / 0.6962 before display clipping.

At the 20 °C limit, recovered source heat per kilogram changes by +36.77% and net power per unit absorbed heat changes by +3.53%. Their ratios multiply to give the +41.59% net-power gain at equal geothermal flow. This comparison uses a 90%-efficient partial turbine and the 88.6%-efficient standard turbine.

| Minimum source outlet °C | 90% turbine: gain % | 80% turbine: gain % | 70% turbine: gain % | 60% turbine: gain % |
|---|---|---|---|---|
| 80 | 2.05 | -9.71 | -21.54 | -33.52 |
| 70 | 6.29 | -6.49 | -19.27 | -32.05 |
| 60 | 17.57 | 2.78 | -12.01 | -26.80 |
| 50 | 30.35 | 13.88 | -2.58 | -19.05 |
| 40 | 40.39 | 22.50 | 4.60 | -13.29 |
| 20 | 41.59 | 23.53 | 5.47 | -12.51 |

At 80% turbine efficiency, the highest sampled source-outlet limit with a positive gain is 60 °C. The largest saved gain is 23.53%.

At 70% turbine efficiency, the highest sampled source-outlet limit with a positive gain is 45 °C. The largest saved gain is 5.47%.

At 60% turbine efficiency, no accepted sampled case beats the 88.6%-efficient standard reference.

Expansion endpoint counts across the 96 partial-evaporation sensitivity cases: wet to wet: 45, dry to dry: 40, wet to dry: 11. These endpoint classes do not certify the phase path inside a real expander.

![Standard-cycle efficiency and achieved heat-source exploitation at 88.6% turbine efficiency.](output/results/sensitivity_butane/baseline_performance.png)

[Embeddable HTML](output/results/sensitivity_butane/baseline_performance.html) · [SVG](output/results/sensitivity_butane/baseline_performance.svg) · [PDF](output/results/sensitivity_butane/baseline_performance.pdf)

![Standard-cycle T-s and T-Q evolution at 88.6% turbine efficiency.](output/results/sensitivity_butane/baseline_cycle_evolution.png)

[Embeddable HTML](output/results/sensitivity_butane/baseline_cycle_evolution.html) · [SVG](output/results/sensitivity_butane/baseline_cycle_evolution.svg) · [PDF](output/results/sensitivity_butane/baseline_cycle_evolution.pdf)

![Partial-evaporation T-s and T-Q evolution at 90% turbine efficiency.](output/results/sensitivity_butane/partial_evaporation_cycle_evolution.png)

[Embeddable HTML](output/results/sensitivity_butane/partial_evaporation_cycle_evolution.html) · [SVG](output/results/sensitivity_butane/partial_evaporation_cycle_evolution.svg) · [PDF](output/results/sensitivity_butane/partial_evaporation_cycle_evolution.pdf)

![Fixed-reference net system efficiency for the partial-evaporation sensitivity grid.](output/results/sensitivity_butane/system_efficiency_contours.png)

[Embeddable HTML](output/results/sensitivity_butane/system_efficiency_contours.html) · [SVG](output/results/sensitivity_butane/system_efficiency_contours.svg) · [PDF](output/results/sensitivity_butane/system_efficiency_contours.pdf)

![Achieved heat-source exploitation on the 140-to-15 °C reference.](output/results/sensitivity_butane/actual_exploitation_contours.png)

[Embeddable HTML](output/results/sensitivity_butane/actual_exploitation_contours.html) · [SVG](output/results/sensitivity_butane/actual_exploitation_contours.svg) · [PDF](output/results/sensitivity_butane/actual_exploitation_contours.pdf)

![Expander inlet vapor fraction. Superheated states are displayed as 100%.](output/results/sensitivity_butane/inlet_quality_contours.png)

[Embeddable HTML](output/results/sensitivity_butane/inlet_quality_contours.html) · [SVG](output/results/sensitivity_butane/inlet_quality_contours.svg) · [PDF](output/results/sensitivity_butane/inlet_quality_contours.pdf)

![Expander outlet vapor fraction. Wet admission does not imply dry discharge.](output/results/sensitivity_butane/outlet_quality_contours.png)

[Embeddable HTML](output/results/sensitivity_butane/outlet_quality_contours.html) · [SVG](output/results/sensitivity_butane/outlet_quality_contours.svg) · [PDF](output/results/sensitivity_butane/outlet_quality_contours.pdf)

![Net-power gain against the same fluid's standard cycle at 88.6% turbine efficiency.](output/results/sensitivity_butane/relative_gain_contours.png)

[Embeddable HTML](output/results/sensitivity_butane/relative_gain_contours.html) · [SVG](output/results/sensitivity_butane/relative_gain_contours.svg) · [PDF](output/results/sensitivity_butane/relative_gain_contours.pdf)

## 6. n-Pentane: simple cycle

At a 70 °C source limit, standard-cycle fixed-reference efficiency is 6.87%. Partial evaporation at 90% turbine efficiency reaches 7.51%, a 9.29% gain over the 88.6%-efficient standard reference. This sensitivity comparison includes the difference in turbine efficiency; the preceding eight-case comparison holds turbine efficiency equal.

At equal 88.6% turbine efficiency and a 70 °C source limit, standard admission is 4.64 bar and 91.42 °C. Partial-evaporation admission is 9.33 bar and 121.39 °C. Flexible admission changes pressure level and the temperature of heat addition along with the vapor fraction.

At the 20 °C limit, standard and 90%-efficient partial cycles produce 36.36 and 55.05 kW per kg/s of source flow. Their achieved source outlets are 67.12 and 35.10 °C. The partial-expansion inlet/outlet vapor fractions are 0.0302 / 0.6426 before display clipping.

At the 20 °C limit, recovered source heat per kilogram changes by +43.49% and net power per unit absorbed heat changes by +5.52%. Their ratios multiply to give the +51.42% net-power gain at equal geothermal flow. This comparison uses a 90%-efficient partial turbine and the 88.6%-efficient standard turbine.

| Minimum source outlet °C | 90% turbine: gain % | 80% turbine: gain % | 70% turbine: gain % | 60% turbine: gain % |
|---|---|---|---|---|
| 80 | 3.69 | -8.14 | -19.97 | -31.80 |
| 70 | 9.29 | -3.47 | -16.22 | -28.97 |
| 60 | 21.10 | 6.75 | -7.74 | -21.96 |
| 50 | 35.40 | 19.29 | 3.29 | -12.76 |
| 40 | 48.82 | 31.15 | 13.49 | -4.17 |
| 20 | 51.42 | 33.44 | 15.50 | -2.44 |

At 80% turbine efficiency, the highest sampled source-outlet limit with a positive gain is 65 °C. The largest saved gain is 33.69%.

At 70% turbine efficiency, the highest sampled source-outlet limit with a positive gain is 50 °C. The largest saved gain is 15.71%.

At 60% turbine efficiency, no accepted sampled case beats the 88.6%-efficient standard reference.

Expansion endpoint counts across the 96 partial-evaporation sensitivity cases: wet to wet: 54, dry to dry: 29, wet to dry: 13. These endpoint classes do not certify the phase path inside a real expander.

![Standard-cycle efficiency and achieved heat-source exploitation at 88.6% turbine efficiency.](output/results/sensitivity_pentane/baseline_performance.png)

[Embeddable HTML](output/results/sensitivity_pentane/baseline_performance.html) · [SVG](output/results/sensitivity_pentane/baseline_performance.svg) · [PDF](output/results/sensitivity_pentane/baseline_performance.pdf)

![Standard-cycle T-s and T-Q evolution at 88.6% turbine efficiency.](output/results/sensitivity_pentane/baseline_cycle_evolution.png)

[Embeddable HTML](output/results/sensitivity_pentane/baseline_cycle_evolution.html) · [SVG](output/results/sensitivity_pentane/baseline_cycle_evolution.svg) · [PDF](output/results/sensitivity_pentane/baseline_cycle_evolution.pdf)

![Partial-evaporation T-s and T-Q evolution at 90% turbine efficiency.](output/results/sensitivity_pentane/partial_evaporation_cycle_evolution.png)

[Embeddable HTML](output/results/sensitivity_pentane/partial_evaporation_cycle_evolution.html) · [SVG](output/results/sensitivity_pentane/partial_evaporation_cycle_evolution.svg) · [PDF](output/results/sensitivity_pentane/partial_evaporation_cycle_evolution.pdf)

![Fixed-reference net system efficiency for the partial-evaporation sensitivity grid.](output/results/sensitivity_pentane/system_efficiency_contours.png)

[Embeddable HTML](output/results/sensitivity_pentane/system_efficiency_contours.html) · [SVG](output/results/sensitivity_pentane/system_efficiency_contours.svg) · [PDF](output/results/sensitivity_pentane/system_efficiency_contours.pdf)

![Achieved heat-source exploitation on the 140-to-15 °C reference.](output/results/sensitivity_pentane/actual_exploitation_contours.png)

[Embeddable HTML](output/results/sensitivity_pentane/actual_exploitation_contours.html) · [SVG](output/results/sensitivity_pentane/actual_exploitation_contours.svg) · [PDF](output/results/sensitivity_pentane/actual_exploitation_contours.pdf)

![Expander inlet vapor fraction. Superheated states are displayed as 100%.](output/results/sensitivity_pentane/inlet_quality_contours.png)

[Embeddable HTML](output/results/sensitivity_pentane/inlet_quality_contours.html) · [SVG](output/results/sensitivity_pentane/inlet_quality_contours.svg) · [PDF](output/results/sensitivity_pentane/inlet_quality_contours.pdf)

![Expander outlet vapor fraction. Wet admission does not imply dry discharge.](output/results/sensitivity_pentane/outlet_quality_contours.png)

[Embeddable HTML](output/results/sensitivity_pentane/outlet_quality_contours.html) · [SVG](output/results/sensitivity_pentane/outlet_quality_contours.svg) · [PDF](output/results/sensitivity_pentane/outlet_quality_contours.pdf)

![Net-power gain against the same fluid's standard cycle at 88.6% turbine efficiency.](output/results/sensitivity_pentane/relative_gain_contours.png)

[Embeddable HTML](output/results/sensitivity_pentane/relative_gain_contours.html) · [SVG](output/results/sensitivity_pentane/relative_gain_contours.svg) · [PDF](output/results/sensitivity_pentane/relative_gain_contours.pdf)

## 7. Cyclopentane: simple cycle

At a 70 °C source limit, standard-cycle fixed-reference efficiency is 6.55%. Partial evaporation at 90% turbine efficiency reaches 7.57%, a 15.55% gain over the 88.6%-efficient standard reference. This sensitivity comparison includes the difference in turbine efficiency; the preceding eight-case comparison holds turbine efficiency equal.

At equal 88.6% turbine efficiency and a 70 °C source limit, standard admission is 2.63 bar and 83.60 °C. Partial-evaporation admission is 6.61 bar and 120.69 °C. Flexible admission changes pressure level and the temperature of heat addition along with the vapor fraction.

At the 20 °C limit, standard and 90%-efficient partial cycles produce 34.44 and 55.24 kW per kg/s of source flow. Their achieved source outlets are 70.01 and 35.17 °C. The partial-expansion inlet/outlet vapor fractions are 0.0191 / 0.4652 before display clipping.

At the 20 °C limit, recovered source heat per kilogram changes by +49.27% and net power per unit absorbed heat changes by +7.48%. Their ratios multiply to give the +60.43% net-power gain at equal geothermal flow. This comparison uses a 90%-efficient partial turbine and the 88.6%-efficient standard turbine.

| Minimum source outlet °C | 90% turbine: gain % | 80% turbine: gain % | 70% turbine: gain % | 60% turbine: gain % |
|---|---|---|---|---|
| 80 | 3.15 | -8.46 | -20.07 | -31.67 |
| 70 | 15.55 | 2.24 | -11.08 | -24.39 |
| 60 | 28.70 | 13.68 | -1.47 | -16.37 |
| 50 | 43.94 | 27.08 | 10.32 | -6.48 |
| 40 | 58.03 | 39.55 | 21.08 | 2.62 |
| 20 | 60.43 | 41.67 | 22.93 | 4.21 |

At 80% turbine efficiency, the highest sampled source-outlet limit with a positive gain is 70 °C. The largest saved gain is 41.67%.

At 70% turbine efficiency, the highest sampled source-outlet limit with a positive gain is 55 °C. The largest saved gain is 22.93%.

At 60% turbine efficiency, the highest sampled source-outlet limit with a positive gain is 40 °C. The largest saved gain is 4.21%.

Expansion endpoint counts across the 96 partial-evaporation sensitivity cases: wet to wet: 56, dry to dry: 34, wet to dry: 6. These endpoint classes do not certify the phase path inside a real expander.

![Standard-cycle efficiency and achieved heat-source exploitation at 88.6% turbine efficiency.](output/results/sensitivity_cyclopentane/baseline_performance.png)

[Embeddable HTML](output/results/sensitivity_cyclopentane/baseline_performance.html) · [SVG](output/results/sensitivity_cyclopentane/baseline_performance.svg) · [PDF](output/results/sensitivity_cyclopentane/baseline_performance.pdf)

![Standard-cycle T-s and T-Q evolution at 88.6% turbine efficiency.](output/results/sensitivity_cyclopentane/baseline_cycle_evolution.png)

[Embeddable HTML](output/results/sensitivity_cyclopentane/baseline_cycle_evolution.html) · [SVG](output/results/sensitivity_cyclopentane/baseline_cycle_evolution.svg) · [PDF](output/results/sensitivity_cyclopentane/baseline_cycle_evolution.pdf)

![Partial-evaporation T-s and T-Q evolution at 90% turbine efficiency.](output/results/sensitivity_cyclopentane/partial_evaporation_cycle_evolution.png)

[Embeddable HTML](output/results/sensitivity_cyclopentane/partial_evaporation_cycle_evolution.html) · [SVG](output/results/sensitivity_cyclopentane/partial_evaporation_cycle_evolution.svg) · [PDF](output/results/sensitivity_cyclopentane/partial_evaporation_cycle_evolution.pdf)

![Fixed-reference net system efficiency for the partial-evaporation sensitivity grid.](output/results/sensitivity_cyclopentane/system_efficiency_contours.png)

[Embeddable HTML](output/results/sensitivity_cyclopentane/system_efficiency_contours.html) · [SVG](output/results/sensitivity_cyclopentane/system_efficiency_contours.svg) · [PDF](output/results/sensitivity_cyclopentane/system_efficiency_contours.pdf)

![Achieved heat-source exploitation on the 140-to-15 °C reference.](output/results/sensitivity_cyclopentane/actual_exploitation_contours.png)

[Embeddable HTML](output/results/sensitivity_cyclopentane/actual_exploitation_contours.html) · [SVG](output/results/sensitivity_cyclopentane/actual_exploitation_contours.svg) · [PDF](output/results/sensitivity_cyclopentane/actual_exploitation_contours.pdf)

![Expander inlet vapor fraction. Superheated states are displayed as 100%.](output/results/sensitivity_cyclopentane/inlet_quality_contours.png)

[Embeddable HTML](output/results/sensitivity_cyclopentane/inlet_quality_contours.html) · [SVG](output/results/sensitivity_cyclopentane/inlet_quality_contours.svg) · [PDF](output/results/sensitivity_cyclopentane/inlet_quality_contours.pdf)

![Expander outlet vapor fraction. Wet admission does not imply dry discharge.](output/results/sensitivity_cyclopentane/outlet_quality_contours.png)

[Embeddable HTML](output/results/sensitivity_cyclopentane/outlet_quality_contours.html) · [SVG](output/results/sensitivity_cyclopentane/outlet_quality_contours.svg) · [PDF](output/results/sensitivity_cyclopentane/outlet_quality_contours.pdf)

![Net-power gain against the same fluid's standard cycle at 88.6% turbine efficiency.](output/results/sensitivity_cyclopentane/relative_gain_contours.png)

[Embeddable HTML](output/results/sensitivity_cyclopentane/relative_gain_contours.html) · [SVG](output/results/sensitivity_cyclopentane/relative_gain_contours.svg) · [PDF](output/results/sensitivity_cyclopentane/relative_gain_contours.pdf)

## 8. Novec 649: recuperated cycle

At a 70 °C source limit, standard-cycle fixed-reference efficiency is 7.16%. Partial evaporation at 90% turbine efficiency reaches 8.14%, a 13.71% gain over the 88.6%-efficient standard reference. This sensitivity comparison includes the difference in turbine efficiency; the preceding eight-case comparison holds turbine efficiency equal.

At equal 88.6% turbine efficiency and a 70 °C source limit, standard admission is 3.55 bar and 92.87 °C. Partial-evaporation admission is 5.88 bar and 111.49 °C. Flexible admission changes pressure level and the temperature of heat addition along with the vapor fraction.

At the 20 °C limit, standard and 90%-efficient partial cycles produce 39.60 and 56.12 kW per kg/s of source flow. Their achieved source outlets are 46.96 and 34.37 °C. The partial-expansion inlet/outlet vapor fractions are 0.0660 / 1.0556 before display clipping.

At the 20 °C limit, recovered source heat per kilogram changes by +13.42% and net power per unit absorbed heat changes by +24.97%. Their ratios multiply to give the +41.73% net-power gain at equal geothermal flow. This comparison uses a 90%-efficient partial turbine and the 88.6%-efficient standard turbine.

| Minimum source outlet °C | Standard effectiveness, 88.6% turbine | Partial turbine efficiency % | Partial effectiveness |
|---|---|---|---|
| 70 | 0.621 | 88.6 | 0.655 |
| 20 | 0.034 | 90 | 0.098 |

Internal recuperation reduces external heat demand and can limit how deeply the geothermal stream is cooled. Optimized effectiveness balances these effects, so its selected value changes with the permitted source outlet and inlet-state formulation.

| Minimum source outlet °C | 90% turbine: gain % | 80% turbine: gain % | 70% turbine: gain % | 60% turbine: gain % |
|---|---|---|---|---|
| 80 | 5.95 | -4.98 | -16.19 | -27.69 |
| 70 | 13.71 | 1.94 | -10.08 | -22.53 |
| 60 | 22.82 | 9.98 | -2.71 | -15.97 |
| 50 | 36.44 | 22.25 | 6.87 | -7.78 |
| 40 | 40.11 | 23.09 | 7.63 | -8.86 |
| 20 | 41.73 | 25.44 | 9.31 | -6.95 |

At 80% turbine efficiency, the highest sampled source-outlet limit with a positive gain is 70 °C. The largest saved gain is 25.44%.

At 70% turbine efficiency, the highest sampled source-outlet limit with a positive gain is 55 °C. The largest saved gain is 9.31%.

At 60% turbine efficiency, no accepted sampled case beats the 88.6%-efficient standard reference.

Expansion endpoint counts across the 96 partial-evaporation sensitivity cases: wet to dry: 64, dry to dry: 32. These endpoint classes do not certify the phase path inside a real expander.

![Standard-cycle efficiency and achieved heat-source exploitation at 88.6% turbine efficiency.](output/results/sensitivity_novec/baseline_performance.png)

[Embeddable HTML](output/results/sensitivity_novec/baseline_performance.html) · [SVG](output/results/sensitivity_novec/baseline_performance.svg) · [PDF](output/results/sensitivity_novec/baseline_performance.pdf)

![Standard-cycle T-s and T-Q evolution at 88.6% turbine efficiency.](output/results/sensitivity_novec/baseline_cycle_evolution.png)

[Embeddable HTML](output/results/sensitivity_novec/baseline_cycle_evolution.html) · [SVG](output/results/sensitivity_novec/baseline_cycle_evolution.svg) · [PDF](output/results/sensitivity_novec/baseline_cycle_evolution.pdf)

![Partial-evaporation T-s and T-Q evolution at 90% turbine efficiency.](output/results/sensitivity_novec/partial_evaporation_cycle_evolution.png)

[Embeddable HTML](output/results/sensitivity_novec/partial_evaporation_cycle_evolution.html) · [SVG](output/results/sensitivity_novec/partial_evaporation_cycle_evolution.svg) · [PDF](output/results/sensitivity_novec/partial_evaporation_cycle_evolution.pdf)

![Fixed-reference net system efficiency for the partial-evaporation sensitivity grid.](output/results/sensitivity_novec/system_efficiency_contours.png)

[Embeddable HTML](output/results/sensitivity_novec/system_efficiency_contours.html) · [SVG](output/results/sensitivity_novec/system_efficiency_contours.svg) · [PDF](output/results/sensitivity_novec/system_efficiency_contours.pdf)

![Achieved heat-source exploitation on the 140-to-15 °C reference.](output/results/sensitivity_novec/actual_exploitation_contours.png)

[Embeddable HTML](output/results/sensitivity_novec/actual_exploitation_contours.html) · [SVG](output/results/sensitivity_novec/actual_exploitation_contours.svg) · [PDF](output/results/sensitivity_novec/actual_exploitation_contours.pdf)

![Expander inlet vapor fraction. Superheated states are displayed as 100%.](output/results/sensitivity_novec/inlet_quality_contours.png)

[Embeddable HTML](output/results/sensitivity_novec/inlet_quality_contours.html) · [SVG](output/results/sensitivity_novec/inlet_quality_contours.svg) · [PDF](output/results/sensitivity_novec/inlet_quality_contours.pdf)

![Expander outlet vapor fraction. Wet admission does not imply dry discharge.](output/results/sensitivity_novec/outlet_quality_contours.png)

[Embeddable HTML](output/results/sensitivity_novec/outlet_quality_contours.html) · [SVG](output/results/sensitivity_novec/outlet_quality_contours.svg) · [PDF](output/results/sensitivity_novec/outlet_quality_contours.pdf)

![Net-power gain against the same fluid's standard cycle at 88.6% turbine efficiency.](output/results/sensitivity_novec/relative_gain_contours.png)

[Embeddable HTML](output/results/sensitivity_novec/relative_gain_contours.html) · [SVG](output/results/sensitivity_novec/relative_gain_contours.svg) · [PDF](output/results/sensitivity_novec/relative_gain_contours.pdf)

## 9. Validation and limits

The consolidated run passes comparison with the saved v2 reference across 484 finite-heat sweep cases, 20 unavailable endpoints and 8 comparison cases. The largest specific-power difference is 0.000e+00% relative and 0.000e+00 kW per kg/s absolute. Efficiencies, powers and source mass flow use 0.01% relative tolerance, pressures use 0.1%, temperatures use 0.02 K absolute tolerance and raw qualities use 0.001. Full metric comparisons and tolerances are recorded in output/results/reproducibility_comparison.json and its CSV companion.

| Fluid | Cases with native KKT records | All strict KKT checks passed |
|---|---|---|
| n-Butane | 121 | 67 |
| n-Pentane | 121 | 68 |
| Cyclopentane | 121 | 60 |
| Novec649 | 121 | 89 |

The native first-order KKT diagnostics use a threshold of five times the solver tolerance. Case acceptance follows the native solver-success, physical-flow, energy-balance and constraint checks used in the v2 reference. The additional first-order diagnostics qualify the reported local solutions, and global optimality remains unproven.

| Fluid | Native ThermoOpt optimizations | Dense heater min K | Dense condenser min K | Dense recuperator min K |
|---|---|---|---|---|
| n-Butane | 121 | 4.162 | 4.433 | n.a. |
| n-Pentane | 121 | 4.115 | 4.406 | n.a. |
| Cyclopentane | 121 | 4.052 | 4.407 | n.a. |
| Novec649 | 121 | 4.276 | 4.411 | 5.000 |

The minimum approaches are enforced on the optimization grids. Selected accepted designs at 135, 70 and 20 °C are reevaluated at 501 locations without reoptimization. Approaches below the prescribed limit reveal sampling error between optimization locations. These are subset diagnostics, not a global worst-case bound. Small gains, precise break-even contours and close fluid rankings need grid refinement before design conclusions. Contour lines interpolate between four sampled turbine efficiencies; the plots are not additional optimizations.

The study compares three simple-fluid architectures and one recuperated architecture. Differences involving Novec combine working-fluid choice and recuperation. Real wet-expander efficiencies, internal phase trajectories, equipment sizes and economics are not predicted by a prescribed constant isentropic efficiency.

## 10. Embedding figures

Every report figure is saved as PNG, PDF, SVG and a self-contained HTML file using the same Matplotlib drawing. The HTML uses inline SVG, requires no network or JavaScript, and scales with its container. Embed the HTML in an iframe or the SVG in an image/object element. These are static vector figures; they do not provide Plotly-style hover or data filtering. The standalone HTML report embeds its figures and data. output/figures.html is a local gallery.
