
## Assignment: Optimization of conventional vs. partial-evaporation ORC

This assignment is part of the *Training42Phase* school and focuses on optimizing Organic Rankine Cycle (ORC) power systems.

You will compare a **conventional ORC** (with superheating before the expansion) to a **partial-evaporation ORC** (two-phase expansion) under identical heat source and sink conditions. The goal is to analyze how expander efficiency affects system performance and to determine the conditions under which partial-evaporation becomes competitive.


### 🎯 Objective

Optimize two types of ORCs for a fixed net power output of 1 MW:
- A **conventional ORC**, with superheated vapor at the expander inlet (5 °C superheating)
- A **partial-evaporation ORC**, with two-phase flow at the expander inlet (vapor quality of 30 %)

You will:
1. Compare system efficiencies under fixed expander isentropic efficiency
2. Perform sensitivity analysis on two-phase expander performance
3. Identify the break-even efficiency between configurations



### ⚙️ Starting point

You are given a baseline configuration file: `ORC_single_phase.yaml`, which includes:
- Working fluid: butane
- Heat source: pressurized water at 140 °C, cooled to a minimum temperature of 20 °C
- Heat sink: water at 20 °C, heated to a maximum temperature of 30 °C
- Fixed net power: 1 MW
- Component efficiencies:
  - Expander: 90 %
  - Pumps: 80 %
  - Pressure drop in heat exchangers: 1 %
- Constraints:
  - Minimum **temperature difference** (i.e., pinch point) in heat exchangers : 5 °C
  - Minimum **subcooling** at pump inlet: 1 °C
  - Minimum **superheating** at expander inlet: 5 °C


### Task 1. Baseline optimization (conventional ORC)

Optimize a **conventional ORC** with superheated vapor at the expander inlet.

🔍 **Goal:** Determine the maximum system efficiency

💡 **Hints**:
- System efficiency = thermal efficiency × heat source utilization
- Use the provided `run_optimization.py` script to run the case.
- After optimization, system efficiency is stored in:
  ```python
  problem.cycle_data["energy_analysis"]["system_efficiency"]
  ```

---

### Task 2. Partial-evaporation ORC optimization

Now switch to a **partial-evaporation** configuration (two-phase expansion).

🔍 **Goal:** Determine the maximum system efficiency assuming the two-phase expander has an isentropic efficiency of 90 %, identical to the conventional ORC case.
🛠️ **Required changes**:

💡 **Hints**:
- Remove the superheating constraint:
  ```yaml
  # - variable: $components.expander.state_in.superheating
  #   type: ">"
  #   value: 5.0
  ```
- Add a vapor quality constraint at the expander inlet:
  ```yaml
  - variable: $components.expander.state_in.Q
    type: "="
    value: 0.3
  ```


### Task 3. Sensitivity analysis

🔍 **Goal:** Analyze how the two-phase expander isentropic efficiency affects the system efficiency.  
Generate a plot of system efficiency vs. expander efficiency, covering the range from 100 % to 10 %  in steps of 5 %.

🛠️ **What to do**:
- Expand your Python script to loop over expander efficiencies (e.g., 1.00 → 0.10)
- For each efficiency:
  - Update the expander efficiency in the configuration dictionary
  - Run the optimization
  - Store the resulting system efficiency for plotting

💡 **Hints**:
- Use `np.arange()` or `np.linspace()` to define the range of expander efficiencies.  
  ⚠️ *Remember that `np.arange(start, stop, step)` excludes the stop value and requires a **negative step** if you're looping from 1.00 down to 0.10.*

- When creating the `ThermodynamicCycleOptimization` object, you can specify a custom output directory to separate results for each case:
  ```python
  cycle = ThermodynamicCycleOptimization(CONFIG_FILE, out_dir=f"results_two_phase_eff_{value:0.2f}")
  ```

- You can update specific parameters in the configuration using the provided `set_config_value()` method:
  ```python
  cycle.set_config_value("problem_formulation.fixed_parameters.expander.efficiency", 0.9)
  ```

- The optimization method accepts an optional initial guess via the x0 argument:
  - If `x0=None`, the optimizer starts from the initial guess provided in the YAML configuration file.
  - You can speed up convergence by reusing the solution from the previous run as the initial guess (`x0`) for the next case:
  ```python
  cycle.run_optimization(x0=x0)
  x0 = cycle.solver.x_final
  ```

- To avoid memory issues when running multiple optimizations in a loop, make sure to close the plot figure after each run:
  ```python
  plt.close(cycle.problem.figure)
  ```
  Alternatively, you can disable plotting altogether by setting `callbacks.plot_cycle: False` in the YAML configuration under `solver_options`:
  ```yaml
  solver_options:
    callbacks:
      plot_cycle: False
  ```  
  This is recommended if you don't need to inspect the plots during automated runs.

- Don’t forget to store the resulting **system efficiency** for each case, so you can later plot efficiency vs. expander performance.


### Task 4. Break-even analysis

From your plot, determine the minimum two-phase expander efficiency at which the partial-evaporation ORC becomes as efficient as the conventional ORC (90 % expander efficiency). Report both the break-even efficiency value and your interpretation of the result.








