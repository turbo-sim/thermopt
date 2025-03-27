## Assignment: Optimization of transcritical CO₂ heat pumps with and without work recovery

This assignment is part of the *Training42Phase* school and focuses on optimizing **transcritical CO₂ heat pump systems**. You will analyze the performance of a CO₂-based cycle in two configurations:
- **Baseline system** with expansion through a Joule-Thomson (JT) valve
- **Modified system** with a two-phase expander recovering work from the throttling process

The goal is to understand how expander efficiency impacts heat pump performance and to quantify the trade-off between complexity (i.e., adding work recovery) and system-level benefits.


### 🎯 Objective

Optimize a **recuperated transcritical CO₂ heat pump** for a fixed heating duty of 35 MW:
- A **baseline case** with an isentropic efficiency of 0 % (i.e., a JT valve)
- A **work recovery case** with a 100 % efficient two-phase expander

You will:
1. Compare COPs of the two configurations
2. Estimate the spouting velocity across the two-phase expander
3. Perform sensitivity analysis on expander efficiency


### ⚙️ Starting point

You are given a baseline configuration file: `transcritical_heat_pump.yaml`, which includes:
- **Working fluid**: carbon dioxide (CO₂)
- **Heat source**: nitrogen, from 283.15 K to 278.15 K
- **Heat sink**: high-pressure water, from 313.15 K to 398.15 K
- **Heating duty**: 35 MW
- **Components**:
  - Expander (JT valve by default): 0 % isentropic efficiency
  - Compressor: 85 % isentropic efficiency
  - Pumps: 80 % isentropic efficiency
- **Constraints**:
  - Minimum temperature difference in heat exchangers: 10 °C
  - Minimum superheating at compressor inlet: 1 °C


### Task 1. Baseline optimization (JT valve configuration)

Optimize the system with a Joule-Thomson isenthalpic valve.

🔍 **Goal:** Determine the maximum achievable COP in the baseline configuration.

💡 **Hints**:
- An isenthalpic process is the same as a process with 0 % isentropic efficiency
- After solving, extract the system COP from:
  ```python
  cycle.problem.cycle_data["energy_analysis"]["COP_heat_pump"]
  ```


### Task 2. Optimization with isentropic expander

Re-run the optimization with a 100 % efficient (i.e., isentropic) two-phase expander to replace the JT valve.

🔍 **Goal:** Determine the improvement in COP when full work recovery is available.

💡 **Hints**:
- Change the expander efficiency in the config:
  ```python
  cycle.set_config_value("problem_formulation.fixed_parameters.expander.efficiency", 1.0)
  ```

- Keep all other parameters and constraints the same to ensure a fair comparison.


### Task 3. Spouting velocity analysis

Compute the **spouting velocity** associated with the enthalpy drop across the two-phase expander from Task 2.

🔍 **Goal:**
- Estimate the theoretical velocity achieved by the fluid during isentropic expansion (spouting velocity)
- Reflect on the physical relevance of this velocity in the context of two-phase turbine design
- Compare the result to erosion threshold velocities for typical turbine materials (you can find information about this in Elliot 1982).
- Is the spouting velocity above or below those thresholds? What does this imply?
- Consider whether the two-phase jets actually impact the rotating turbine blades at spouting velocity. Justify your answer

💡 **Hints**:
- The **spouting velocity** represents the ideal velocity reached if the fluid expands **isentropically** from the inlet state to the outlet pressure:
  $$
  v = \sqrt{2 \cdot (h_{\text{in}} - h_{\text{out}})}
  $$
  (with $h$ in J/kg and $v$ in m/s)

- You can extract the necessary enthalpy values from the simulation results:
  ```python
  h_in = cycle.problem.cycle_data["components"]["expander"]["state_in"]["h"]
  h_out = cycle.problem.cycle_data["components"]["expander"]["state_out"]["h"]
  ```



### Task 4. Sensitivity analysis – Expander efficiency

Study how the **isentropic efficiency of the two-phase expander** affects system performance.

🔍 **Goal:** Plot:
- **COP** of the heat pump vs. expander efficiency
- **Work recovery ratio** (expander power / compressor power) vs. expander efficiency

🛠️ **Steps**:
- Loop over expander efficiencies from 1.00 to 0.10 in steps of 0.05
- For each case:
  - Update the expander efficiency in the config
  - Run the optimization
  - Extract and store:
    - `COP_heat_pump`
    - `compressor_power`
    - `expander_power`

💡 **Hints**:
- Use `np.arange(1.00, 0.10, -0.05)` to define the efficiency range
- Use the final solution of the previous case as the initial guess (`x0`) to improve convergence
- Avoid memory issues by closing each figure after the run:
  ```python
  plt.close(cycle.problem.figure)
  ```
