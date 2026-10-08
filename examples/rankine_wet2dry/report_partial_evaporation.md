# Partial evaporation in geothermal Rankine cycles: butane, cyclopentane and Novec 649

Partial evaporation in organic Rankine cycles can improve heat recovery from sources with a temperature glide by reducing the isothermal evaporation duty and modifying the temperature profile of heat addition. However, these benefits must compensate for any reduction in expander efficiency associated with two-phase operation. This study examines how the minimum permissible heat-source outlet temperature and expander efficiency determine when partial evaporation increases net power production. Simple cycles using n-butane and cyclopentane and a recuperated cycle using Novec 649 are optimized for a 140 °C water source and 20 °C cooling-water inlet. A sensitivity analysis comprising 360 optimized cases varies the minimum source outlet temperature from 135 to 20 °C and the expander efficiency from 60% to 90%. Each configuration is compared with its corresponding superheated cycle at 90% expander efficiency, using net power per unit source flow as the comparison basis. At equal expander efficiencies of 90% and a minimum source outlet temperature of 70 °C, partial evaporation increases net power by 6.5%, 15.2%, and 13.1% for n-butane, cyclopentane, and Novec 649, respectively. Relaxing this limit to 20 °C increases the gains to 40.8%, 59.3%, and 42.7%, through combined improvements in heat recovery and heat conversion. At this limit, gains remain positive with an 80%-efficient expander, reaching 22.8%, 40.6%, and 26.2%, respectively; only cyclopentane retains a positive gain at 60% efficiency. The corresponding designs at 90% efficiency require wet-to-wet expansion for n-butane and cyclopentane, whereas Novec 649 achieves wet-to-dry expansion. These results identify deeper source cooling as a major driver of the benefit of partial evaporation and quantify the expander performance needed to retain that benefit. 

- [Partial evaporation in geothermal Rankine cycles: butane, cyclopentane and Novec 649](#partial-evaporation-in-geothermal-rankine-cycles-butane-cyclopentane-and-novec-649)
  - [1. Boundary conditions and optimization formulation](#1-boundary-conditions-and-optimization-formulation)
  - [2. Comparison of fluid and cycle configurations](#2-comparison-of-fluid-and-cycle-configurations)
  - [3. Sensitivity rationale and initialization](#3-sensitivity-rationale-and-initialization)
  - [4. Butane: simple cycle](#4-butane-simple-cycle)
  - [5. Cyclopentane: simple cycle](#5-cyclopentane-simple-cycle)
  - [6. Novec 649: recuperated cycle](#6-novec-649-recuperated-cycle)
  - [7. Cross-fluid interpretation and limits](#7-cross-fluid-interpretation-and-limits)

## 1. Boundary conditions and optimization formulation

### 1.1 Common conditions

| Parameter | Value |
|---|---|
| Source water | 140 °C, 10 bar |
| Source minimum outlet | 70 °C for the six-case comparison; 140–20 °C in 5 K steps for sensitivities |
| Cooling water | 20 °C inlet; 24.99–25.00 °C outlet |
| Net cycle power used to size every design | 1 MW after the working-fluid pump, before external-water pumps |
| Working fluids and architectures | Simple n-butane; simple cyclopentane; recuperated Novec 649 |
| Expander efficiency | 90% for superheated baselines; 90%, 80%, 70%, 60% for partial evaporation |
| Pump efficiencies | 70% |
| Exchanger pressure losses | 1% per active side |
| Minimum exchanger approach | 5 K at 100 evaporator/condenser locations; 25 recuperator locations |
| Working-fluid pump inlet | At least 1 K subcooling |
| Superheated formulation | 2 K inlet superheat and dry expander outlet |
| Partial-evaporation formulation | Optimized inlet state; no fixed quality or dry-outlet constraint |

### 1.2 Optimization variables

| Variable | Role |
|---|---|
| Source outlet temperature | Amount of source cooling, subject to the reinjection limit |
| Cooling-water outlet temperature | Approximately 5 K total warming |
| Working-fluid pump inlet pressure and enthalpy | Low-pressure state and subcooling |
| Expander inlet pressure and enthalpy | High-pressure state and inlet phase |
| Recuperator effectiveness | 0–1, Novec only |

Mass flows follow from the states, heat balances and 1 MW sizing target. Pressure and enthalpy bounds are broad and property-based. Each initial state comes from a saved coarse-grid optimum; the partial-evaporation formulation does not constrain the inlet quality.

### 1.3 Fixed-reference efficiency and comparison basis

The efficiency contours and baseline performance curves use net system power after all pumps divided by the source heat available between **140 and 20 °C**, at the source pressure of 10 bar:

```text
ηfixed = net system power / [source mass flow × (hwater(140 °C, 10 bar) − hwater(20 °C, 10 bar))]
```

The reference enthalpy drop is **504.723 kJ/kg** for every case. This efficiency is directly proportional to net system power per source flow, so deeper permissible cooling is assessed on the same energy basis. The 1 MW cycle sizing does not imply equal geothermal flow.

## 2. Comparison of fluid and cycle configurations

All six designs use a 70 °C reinjection limit and 90% expander efficiency. This table uses **superheated cyclopentane** as the common baseline; the sensitivities below use each fluid's own baseline.

| Case | Inlet vapor [%] | Outlet vapor [%] | Source outlet [°C] | Gain vs cyclopentane baseline [%] | Variable-reference system efficiency [%] | Fixed-reference system efficiency [%] |
|---|---:|---:|---:|---:|---:|---:|
| cyclopentane superheated | 100.00 | 100.00 | 70.11 | +0.00 | 11.5486 | 6.7663 |
| cyclopentane partial evaporation | 15.95 | 51.32 | 70.00 | +15.24 | 13.3084 | 7.7974 |
| butane superheated | 100.00 | 100.00 | 70.48 | +6.17 | 12.2615 | 7.1840 |
| butane partial evaporation | 57.00 | 85.40 | 70.00 | +13.03 | 13.0528 | 7.6476 |
| novec superheated | 100.00 | 100.00 | 70.00 | +8.92 | 12.5792 | 7.3702 |
| novec partial evaporation | 60.80 | 100.00 | 70.00 | +23.20 | 14.2273 | 8.3358 |

Each fluid has a separate 2 × 2 comparison: superheated and partial-evaporation designs occupy the left and right columns, with T–s diagrams above T–Q diagrams.

![Cyclopentane: simple cycle: superheated and partial evaporation, T–s above T–Q](results/fluid_comparison/comparison_grid_cyclopentane.png)

![Butane: simple cycle: superheated and partial evaporation, T–s above T–Q](results/fluid_comparison/comparison_grid_butane.png)

![Novec 649: recuperated cycle: superheated and partial evaporation, T–s above T–Q](results/fluid_comparison/comparison_grid_novec.png)

Relative to each fluid's own superheated design, the partial-evaporation gains at this limit are **+6.45% for Butane**, **+15.24% for Cyclopentane**, **+13.10% for Novec 649**. The highest fixed-reference efficiency is **8.34%**, for **novec partial evaporation**. The optimized source outlets differ, so the comparison includes both heat-recovery and heat-conversion effects.


## 3. Sensitivity rationale and initialization

### 3.1 Physical question and grid

Wet admission can reduce the isothermal evaporation duty and alter source matching, while an efficiency penalty reduces work recovered in expansion. The sensitivity tests where the former benefit outweighs the latter. There are **24 finite-heat temperatures × five designs = 120 optimized points per fluid**, or **360 total**. The 140 °C limit permits no geothermal temperature drop and makes ThermoOpt's available-heat denominator zero. Its 15 requested combinations are recorded as unavailable, not as zero-efficiency cycles.



### 3.2 Plot conventions

Nominal exploitation is `100 × (140 − minimum outlet) / 120`; actual exploitation replaces the minimum with the optimized outlet. Both use 20 °C as the reference. Quality is shown in percent: 0% includes liquid and 100% includes superheated vapor. Dark blue indicates more liquid. Contour lines label the interpolated bands. The T–s/T–Q overlays use every other temperature at 90% expander efficiency, thin lines, and magma restricted to colormap coordinates 0.25–0.75. Their native cumulative heat-duty coordinate uses the 1 MW sizing basis.

## 4. Butane: simple cycle

### 4.1 Superheated baseline

![Butane: simple cycle: baseline net efficiency on a fixed heat reference and actual exploitation](results/sensitivity_butane/baseline_performance.png)

At 90% expander efficiency, the baseline's fixed-reference net efficiency changes from **0.69%** at the 135 °C reinjection limit to **7.31%** at 20 °C. Its corresponding net power per unit source flow changes from **3.48 to 36.88 kW per kg/s**. These two measures have exactly the same trend because their heat reference is fixed.

At the 20 °C limit, the baseline actually leaves the source at **63.41 °C**, or **63.82% actual exploitation**. Permission for further cooling does not require the optimizer to use it: the selected pressure levels, phase-change temperatures and exchanger approaches determine useful recovery.

![Butane: simple cycle: superheated T–s and T–Q evolution at 90% expander efficiency](results/sensitivity_butane/baseline_cycle_evolution.png)

The overlays show how the pressure levels and heat-addition profile change as the reinjection limit is relaxed. 

### 4.2 Evolution with partial evaporation

![Butane: simple cycle: optimized-inlet T–s and T–Q evolution at 90% efficiency](results/sensitivity_butane/partial_evaporation_cycle_evolution.png)

The optimized inlet has **57.00% vapor** at a 70 °C limit and **9.90%** at 20 °C. Values below 100% shorten the evaporation portion of external heating. The balance between sensible heating and evaporation can then match the cooling source differently from the superheated cycle. 

At the 70 °C limit, the partial-evaporation design admits fluid at **17.82 bar and 108.11 °C**, compared with **13.40 bar and 95.45 °C** for the baseline. Inlet-state flexibility therefore changes pressure level and temperature of heat addition, not just the amount of liquid entering the expander.

### 4.3 Net efficiency on the common heat basis

![Butane: simple cycle: fixed-reference net-efficiency contours](results/sensitivity_butane/system_efficiency_contours.png)

At 90% expander efficiency, the plotted efficiency is **0.69%** at a 135 °C limit, **7.65%** at 70 °C and **10.29%** at 20 °C. Unlike efficiency divided by heat available only down to each reinjection limit, this figure directly tracks net electricity per unit geothermal flow across the whole horizontal axis. Any plateau means that permitting further source cooling adds little useful power.

### 4.4 Actual source exploitation

![Butane: simple cycle: actual heat-source exploitation](results/sensitivity_butane/actual_exploitation_contours.png)

At the 20 °C limit, the four partial-evaporation designs actually leave the source at **34.66–35.78 °C**, giving **86.85–87.78%** actual exploitation. The map distinguishes allowed cooling from cooling selected by the optimizer. More recovered heat alone does not establish greater power: its temperature and the expansion efficiency also matter.

### 4.5 Expander inlet and outlet quality

![Butane: simple cycle: expander inlet vapor quality in percent](results/sensitivity_butane/inlet_quality_contours.png)

Across the 96 partial-formulation designs, **64** have two-phase inlets, **32** have dry-vapor inlets and **0** have liquid inlets. The formulation optimizes inlet state rather than imposing a wet inlet. A displayed value of 100% includes superheated vapor.

![Butane: simple cycle: expander outlet vapor quality in percent](results/sensitivity_butane/outlet_quality_contours.png)

| Expansion endpoint phases | Number of designs |
|---|---:|
| wet to wet | 45 |
| dry to dry | 32 |
| wet to dry | 19 |

At 90% efficiency, displayed outlet quality is **85.40%** at the 70 °C limit and **69.11%** at 20 °C. The deep-recovery design remains two-phase at discharge; its advantage therefore requires an expander suited to substantial liquid at both ends.

### 4.6 Relative gain and its physical origin

![Butane: simple cycle: relative net-power gain against the same fluid's 90%-efficient superheated baseline](results/sensitivity_butane/relative_gain_contours.png)

At the **70 °C** limit and 90% expander efficiency, partial evaporation changes fixed-reference efficiency from **7.18% to 7.65%**, a **+6.45% relative gain**. The actual source outlets are **70.48 °C** for the baseline and **70.00 °C** for partial evaporation. The differing source outlets show that heat recovery also contributes; the gain cannot be interpreted as a change in conversion efficiency alone.

At the **20 °C** limit, partial evaporation at 90% efficiency changes recovered source heat per kilogram by **+37.17%** and net system power per unit absorbed heat by **+2.66%**. These factors multiply, giving **+40.82% net-power gain** at equal geothermal flow. This separates the value of additional source cooling from the value of converting the recovered heat.

| Minimum source outlet [°C] | Nominal exploitation [%] | Gain at expander efficiency 90% | 80% | 70% | 60% |
|---:|---:|---:|---:|---:|---:|
| 80 | 50.00 | +1.23% | -10.69% | -22.65% | -34.21% |
| 70 | 58.33 | +6.45% | -6.36% | -19.18% | -31.99% |
| 60 | 66.67 | +17.86% | +2.99% | -11.87% | -26.74% |
| 50 | 75.00 | +30.61% | +14.06% | -2.49% | -19.04% |
| 40 | 83.33 | +40.13% | +22.19% | +4.25% | -13.69% |
| 20 | 100.00 | +40.82% | +22.81% | +4.80% | -13.09% |

### 4.7 Efficiency penalty and break-even conditions

At **80% expander efficiency**, the largest modeled gain is **+22.81%**. The highest sampled reinjection limit with positive gain is **60 °C**, where the gain is **+2.99%**.

At **70% expander efficiency**, the largest modeled gain is **+4.80%**. The highest sampled reinjection limit with positive gain is **45 °C**, where the gain is **+1.52%**.

At **60% expander efficiency**, the largest modeled gain is **-13.08%**. No sampled reinjection limit beats the 90%-efficient superheated baseline.

The zero contour is an approximate break-even boundary between the sampled efficiencies. A positive gain means the inlet-state and heat-recovery benefits offset the prescribed expander penalty. It does not demonstrate that a real expander can attain that efficiency at the computed inlet and outlet qualities. Near-zero gains and abrupt phase transitions are especially sensitive to local optimization and exchanger discretization.

The largest decrease in saved power per source flow relative to an earlier, stricter outlet limit is **0.01%**. Relaxing the limit preserves the earlier design's feasibility, so such a decrease is a local-solution irregularity rather than a thermodynamic penalty for allowing more cooling. It is retained under the single-start protocol.

At equal 90% expander efficiency, the partial formulation falls below the superheated baseline by more than 0.01% at 1 sampled outlet limit (largest loss **0.04%**). The flexible formulation includes the baseline's inlet state, so these negative values identify local-solution limitations, not an inherent disadvantage of allowing partial evaporation.

## 5. Cyclopentane: simple cycle

### 5.1 Superheated baseline

![Cyclopentane: simple cycle: baseline net efficiency on a fixed heat reference and actual exploitation](results/sensitivity_cyclopentane/baseline_performance.png)

At 90% expander efficiency, the baseline's fixed-reference net efficiency changes from **0.78%** at the 135 °C reinjection limit to **6.75%** at 20 °C. Its corresponding net power per unit source flow changes from **3.95 to 34.05 kW per kg/s**. 

At the 20 °C limit, the baseline actually leaves the source at **68.95 °C**, or **59.21% actual exploitation**. Permission for further cooling does not require the optimizer to use it: the selected pressure levels, phase-change temperatures and exchanger approaches determine useful recovery.

![Cyclopentane: simple cycle: superheated T–s and T–Q evolution at 90% expander efficiency](results/sensitivity_cyclopentane/baseline_cycle_evolution.png)

The overlays show how the pressure levels and heat-addition profile change as the reinjection limit is relaxed. Coincident low-limit curves indicate that a similar design is selected despite allowing more source cooling.

### 5.2 Evolution with partial evaporation

![Cyclopentane: simple cycle: optimized-inlet T–s and T–Q evolution at 90% efficiency](results/sensitivity_cyclopentane/partial_evaporation_cycle_evolution.png)

The optimized inlet has **15.95% vapor** at a 70 °C limit and **2.05%** at 20 °C. Values below 100% shorten the evaporation portion of external heating. The balance between sensible heating and evaporation can then match the cooling source differently from the superheated cycle. The expansion outlet must also be examined; wet admission does not imply dry discharge.

At the 70 °C limit, the partial-evaporation design admits fluid at **6.54 bar and 120.19 °C**, compared with **2.60 bar and 83.17 °C** for the baseline. Inlet-state flexibility therefore changes pressure level and temperature of heat addition, not just the amount of liquid entering the expander.

### 5.3 Net efficiency on the common heat basis

![Cyclopentane: simple cycle: fixed-reference net-efficiency contours](results/sensitivity_cyclopentane/system_efficiency_contours.png)

At 90% expander efficiency, the plotted efficiency is **0.78%** at a 135 °C limit, **7.80%** at 70 °C and **10.74%** at 20 °C. Unlike efficiency divided by heat available only down to each reinjection limit, this figure directly tracks net electricity per unit geothermal flow across the whole horizontal axis. Any plateau means that permitting further source cooling adds little useful power.

### 5.4 Actual source exploitation

![Cyclopentane: simple cycle: actual heat-source exploitation](results/sensitivity_cyclopentane/actual_exploitation_contours.png)

At the 20 °C limit, the four partial-evaporation designs actually leave the source at **35.46–36.62 °C**, giving **86.15–87.12%** actual exploitation. The map distinguishes allowed cooling from cooling selected by the optimizer. More recovered heat alone does not establish greater power: its temperature and the expansion efficiency also matter.

### 5.5 Expander inlet and outlet quality

![Cyclopentane: simple cycle: expander inlet vapor quality in percent](results/sensitivity_cyclopentane/inlet_quality_contours.png)

Across the 96 partial-formulation designs, **64** have two-phase inlets, **32** have dry-vapor inlets and **0** have liquid inlets. The formulation optimizes inlet state rather than imposing a wet inlet. A displayed value of 100% includes superheated vapor.

![Cyclopentane: simple cycle: expander outlet vapor quality in percent](results/sensitivity_cyclopentane/outlet_quality_contours.png)

| Expansion endpoint phases | Number of designs |
|---|---:|
| wet to wet | 56 |
| dry to dry | 32 |
| wet to dry | 8 |

At 90% efficiency, displayed outlet quality is **51.32%** at the 70 °C limit and **46.14%** at 20 °C. The deep-recovery design remains two-phase at discharge; its advantage therefore requires an expander suited to substantial liquid at both ends.

### 5.6 Relative gain and its physical origin

![Cyclopentane: simple cycle: relative net-power gain against the same fluid's 90%-efficient superheated baseline](results/sensitivity_cyclopentane/relative_gain_contours.png)

At the **70 °C** limit and 90% expander efficiency, partial evaporation changes fixed-reference efficiency from **6.77% to 7.80%**, a **+15.24% relative gain**. The actual source outlets are **70.11 °C** for the baseline and **70.00 °C** for partial evaporation. The differing source outlets show that heat recovery also contributes; the gain cannot be interpreted as a change in conversion efficiency alone.

At the **20 °C** limit, partial evaporation at 90% efficiency changes recovered source heat per kilogram by **+45.04%** and net system power per unit absorbed heat by **+9.81%**. These factors multiply, giving **+59.27% net-power gain** at equal geothermal flow. This separates the value of additional source cooling from the value of converting the recovered heat.

| Minimum source outlet [°C] | Nominal exploitation [%] | Gain at expander efficiency 90% | 80% | 70% | 60% |
|---:|---:|---:|---:|---:|---:|
| 80 | 50.00 | +2.28% | -9.23% | -20.75% | -32.26% |
| 70 | 58.33 | +15.24% | +1.95% | -11.34% | -24.62% |
| 60 | 66.67 | +28.90% | +13.83% | -1.36% | -16.30% |
| 50 | 75.00 | +44.09% | +27.24% | +10.44% | -6.45% |
| 40 | 83.33 | +57.76% | +39.28% | +20.83% | +2.36% |
| 20 | 100.00 | +59.27% | +40.63% | +22.02% | +3.41% |

### 5.7 Efficiency penalty and break-even conditions

At **80% expander efficiency**, the largest modeled gain is **+40.64%**. The highest sampled reinjection limit with positive gain is **70 °C**, where the gain is **+1.95%**.

At **70% expander efficiency**, the largest modeled gain is **+22.02%**. The highest sampled reinjection limit with positive gain is **55 °C**, where the gain is **+4.62%**.

At **60% expander efficiency**, the largest modeled gain is **+3.42%**. The highest sampled reinjection limit with positive gain is **40 °C**, where the gain is **+2.36%**.

The zero contour is an approximate break-even boundary between the sampled efficiencies. A positive gain means the inlet-state and heat-recovery benefits offset the prescribed expander penalty. It does not demonstrate that a real expander can attain that efficiency at the computed inlet and outlet qualities. Near-zero gains and abrupt phase transitions are especially sensitive to local optimization and exchanger discretization.

The largest decrease in saved power per source flow relative to an earlier, stricter outlet limit is **0.33%**. Relaxing the limit preserves the earlier design's feasibility, so such a decrease is a local-solution irregularity rather than a thermodynamic penalty for allowing more cooling. It is retained under the single-start protocol.

At equal 90% expander efficiency, the partial formulation falls below the superheated baseline by more than 0.01% at 1 sampled outlet limit (largest loss **0.06%**). The flexible formulation includes the baseline's inlet state, so these negative values identify local-solution limitations, not an inherent disadvantage of allowing partial evaporation.

## 6. Novec 649: recuperated cycle

### 6.1 Superheated baseline

![Novec 649: recuperated cycle: baseline net efficiency on a fixed heat reference and actual exploitation](results/sensitivity_novec/baseline_performance.png)

At 90% expander efficiency, the baseline's fixed-reference net efficiency changes from **0.80%** at the 135 °C reinjection limit to **7.57%** at 20 °C. Its corresponding net power per unit source flow changes from **4.05 to 38.19 kW per kg/s**.

At the 20 °C limit, the baseline actually leaves the source at **47.50 °C**, or **77.09% actual exploitation**. Permission for further cooling does not require the optimizer to use it: the selected pressure levels, phase-change temperatures and exchanger approaches determine useful recovery.

![Novec 649: recuperated cycle: superheated T–s and T–Q evolution at 90% expander efficiency](results/sensitivity_novec/baseline_cycle_evolution.png)

The overlays show how the pressure levels and heat-addition profile change as the reinjection limit is relaxed. Coincident low-limit curves indicate that a similar design is selected despite allowing more source cooling.

### 6.2 Evolution with partial evaporation

![Novec 649: recuperated cycle: optimized-inlet T–s and T–Q evolution at 90% efficiency](results/sensitivity_novec/partial_evaporation_cycle_evolution.png)

The optimized inlet has **60.80% vapor** at a 70 °C limit and **5.58%** at 20 °C. Values below 100% shorten the evaporation portion of external heating. The balance between sensible heating and evaporation can then match the cooling source differently from the superheated cycle. The expansion outlet must also be examined; wet admission does not imply dry discharge.

At the 70 °C limit, the partial-evaporation design admits fluid at **5.81 bar and 110.99 °C**, compared with **3.51 bar and 92.37 °C** for the baseline. Inlet-state flexibility therefore changes pressure level and temperature of heat addition, not just the amount of liquid entering the expander.

For Novec, recuperator effectiveness is **0.617/0.648** for baseline/partial evaporation at 70 °C and **0.033/0.097** at 20 °C. Internal preheating reduces external heat demand, but can also prevent deep cooling of the geothermal stream. The optimizer balances these effects; a recuperated architecture need not use high effectiveness at every condition.

### 6.3 Net efficiency on the common heat basis

![Novec 649: recuperated cycle: fixed-reference net-efficiency contours](results/sensitivity_novec/system_efficiency_contours.png)

At 90% expander efficiency, the plotted efficiency is **0.80%** at a 135 °C limit, **8.34%** at 70 °C and **10.79%** at 20 °C. Unlike efficiency divided by heat available only down to each reinjection limit, this figure directly tracks net electricity per unit geothermal flow across the whole horizontal axis. Any plateau means that permitting further source cooling adds little useful power.

### 6.4 Actual source exploitation

![Novec 649: recuperated cycle: actual heat-source exploitation](results/sensitivity_novec/actual_exploitation_contours.png)

At the 20 °C limit, the four partial-evaporation designs actually leave the source at **31.54–34.04 °C**, giving **88.30–90.38%** actual exploitation. The map distinguishes allowed cooling from cooling selected by the optimizer. More recovered heat alone does not establish greater power: its temperature and the expansion efficiency also matter.

### 6.5 Expander inlet and outlet quality

![Novec 649: recuperated cycle: expander inlet vapor quality in percent](results/sensitivity_novec/inlet_quality_contours.png)

Across the 96 partial-formulation designs, **64** have two-phase inlets, **32** have dry-vapor inlets and **0** have liquid inlets. The formulation optimizes inlet state rather than imposing a wet inlet. A displayed value of 100% includes superheated vapor.

![Novec 649: recuperated cycle: expander outlet vapor quality in percent](results/sensitivity_novec/outlet_quality_contours.png)

| Expansion endpoint phases | Number of designs |
|---|---:|
| wet to dry | 64 |
| dry to dry | 32 |

At 90% efficiency, displayed outlet quality is **100.00%** at the 70 °C limit and **100.00%** at 20 °C. The deep-recovery design reaches dry vapor at discharge. If its inlet is wet, this is a wet-to-dry expansion; the dry outlet does not remove the need to accommodate liquid admission.

### 6.6 Relative gain and its physical origin

![Novec 649: recuperated cycle: relative net-power gain against the same fluid's 90%-efficient superheated baseline](results/sensitivity_novec/relative_gain_contours.png)

At the **70 °C** limit and 90% expander efficiency, partial evaporation changes fixed-reference efficiency from **7.37% to 8.34%**, a **+13.10% relative gain**. The actual source outlets are **70.00 °C** for the baseline and **70.00 °C** for partial evaporation. Because both reach essentially the same outlet, this gain comes from converting the recovered heat more effectively, not from a larger recovered heat quantity.

At the **20 °C** limit, partial evaporation at 90% efficiency changes recovered source heat per kilogram by **+15.66%** and net system power per unit absorbed heat by **+23.35%**. These factors multiply, giving **+42.66% net-power gain** at equal geothermal flow. This separates the value of additional source cooling from the value of converting the recovered heat.

| Minimum source outlet [°C] | Nominal exploitation [%] | Gain at expander efficiency 90% | 80% | 70% | 60% |
|---:|---:|---:|---:|---:|---:|
| 80 | 50.00 | +4.46% | -6.35% | -17.43% | -28.79% |
| 70 | 58.33 | +13.10% | +1.37% | -10.61% | -23.01% |
| 60 | 66.67 | +23.11% | +10.22% | -2.51% | -15.79% |
| 50 | 75.00 | +37.11% | +22.56% | +7.26% | -7.52% |
| 40 | 83.33 | +42.14% | +25.06% | +8.93% | -7.77% |
| 20 | 100.00 | +42.66% | +26.22% | +9.74% | -6.80% |

### 6.7 Efficiency penalty and break-even conditions

At **80% expander efficiency**, the largest modeled gain is **+26.22%**. The highest sampled reinjection limit with positive gain is **70 °C**, where the gain is **+1.37%**.

At **70% expander efficiency**, the largest modeled gain is **+9.75%**. The highest sampled reinjection limit with positive gain is **55 °C**, where the gain is **+2.78%**.

At **60% expander efficiency**, the largest modeled gain is **-6.80%**. No sampled reinjection limit beats the 90%-efficient superheated baseline.

The zero contour is an approximate break-even boundary between the sampled efficiencies. A positive gain means the inlet-state and heat-recovery benefits offset the prescribed expander penalty. It does not demonstrate that a real expander can attain that efficiency at the computed inlet and outlet qualities. Near-zero gains and abrupt phase transitions are especially sensitive to local optimization and exchanger discretization.

The largest decrease in saved power per source flow relative to an earlier, stricter outlet limit is **0.07%**. Relaxing the limit preserves the earlier design's feasibility, so such a decrease is a local-solution irregularity rather than a thermodynamic penalty for allowing more cooling. It is retained under the single-start protocol.

## 7. Cross-fluid interpretation and limits

At a 20 °C reinjection limit and 90% expander efficiency:

| Configuration | Baseline power per source flow [kW per kg/s] | Partial power per source flow [kW per kg/s] | Gain [%] | Partial fixed-reference efficiency [%] |
|---|---:|---:|---:|---:|
| Butane: simple cycle | 36.88 | 51.93 | +40.82 | 10.29 |
| Cyclopentane: simple cycle | 34.05 | 54.22 | +59.27 | 10.74 |
| Novec 649: recuperated cycle | 38.19 | 54.48 | +42.66 | 10.79 |

Relative gain and absolute electricity production answer different questions. A larger gain can reflect a weaker baseline rather than a better final cycle. The fixed-reference efficiency columns permit direct comparison of power per geothermal flow across both fluids and reinjection limits. 

### 7.1 Remaining sampling and optimization limits

The **5 K approach is imposed at 100 evaporator/condenser locations and 25 recuperator locations, not continuously**. Nine selected designs per fluid and their coarse-grid seeds were reevaluated at 501 locations without reoptimization:

| Fluid | Coarse seed: lowest checked heater approach [K] | Fine result: lowest checked heater approach [K] | Fine result: lowest checked condenser approach [K] | Fine result: lowest checked recuperator approach [K] |
|---|---:|---:|---:|---:|
| Butane | 3.451 | 4.658 | 4.968 | — |
| Cyclopentane | 2.983 | 4.553 | 4.970 | — |
| Novec 649 | 4.118 | 4.780 | 4.955 | 5.000 |

These subset diagnostics do not locate the worst approach over all designs. The reported results retain the requested 100-element evaporator and condenser model. Approaches below 5 K mean that a finer discretization and reoptimization are needed before treating precise break-even efficiencies or small gains as design conclusions. Close fluid rankings should be read with the same caution.

