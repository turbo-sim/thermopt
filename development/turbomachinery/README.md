# Functional radial turbomachinery models

The module implements design-point correlations from [Parisi et al., Energy Conversion and Management 342 (2025), 120085](https://doi.org/10.1016/j.enconman.2025.120085). Its functions accept an explicit jaxprop fluid and return their results. The former machine classes, mutable machine state, `set_fluid()` and temperature-iteration branch have been removed.

## Use

```python
import jaxprop as props
from thermopt.components.turbomachinery_nondimensional import solve_compressor, solve_turbine

fluid = props.Fluid("Air", backend="HEOS")
inlet = fluid.get_state(props.PT_INPUTS, 1e5, 300.0)
result = solve_compressor(
    fluid, inlet, p_out=2e5, mass_flow=2.0,
    design={"flow_coefficient": 0.07, "backsweep": 40.0},
)
print(result["state_out"].T)
print(result["performance"])
print(result["residual"], result["solver"])
```

`solve_turbine()` takes the same arguments, with an outlet pressure below the inlet pressure. `compressor_defaults()` and `turbine_defaults()` return fresh design dictionaries. Partial overrides are supported. Each result contains inlet, outlet and isentropic-outlet FluidStates, performance, completed design, efficiency residual and solver diagnostics.

Performance contains efficiencies as fractions, positive consumed/generated power [W], diameter [m], peripheral speed [m/s], angular speed [rad/s] and estimated mechanical stress [Pa]. Compressor `correlation_efficiency` is polytropic; turbine `correlation_efficiency` is overall isentropic. The turbine's specific speed uses **rotor-exit static density**, not overall total-outlet density.

Every function in the model has a docstring explaining its inputs, physical role, units, assumptions and applicable paper equations. The model stays single phase: enthalpy iteration does not establish validity for wet machinery. Published study ranges are documented; the solver checks basic input domains but does not silently clamp extrapolated designs.

## Run the demonstrations

From the thermopt root:

```powershell
.\.venv\Scripts\python.exe development/turbomachinery/demo_compressor.py
.\.venv\Scripts\python.exe development/turbomachinery/demo_turbine.py
```

Each writes `trends.csv`, `metadata.json`, `efficiency_trends.png` and `efficiency_trends.svg` under its own ignored `results/` directory. Use `--points 7` for smaller sweeps or `--output PATH` to choose a destination. Scripts work independently of the calling directory and do not write files on import.

- Compressor: flow coefficient versus polytropic efficiency, comparing pressure ratios, mass flows, vaned/vaneless diffusers, backsweep and clearance.
- Turbine: rotor-exit specific speed versus overall isentropic efficiency, comparing mass flows, scalloping, diffuser loss and clearance; also efficiency versus rpm.

These are design-point sweeps, not measured off-design maps or reproductions of the paper's raw validation datasets. A failed solve stops the script instead of dropping the point.

## Reproduce the paper figures

```powershell
.\.venv\Scripts\python.exe development/turbomachinery/paper_figures.py
```

This reproduces component Figures 1-7 with production correlation functions and independently extracted published curves/mean-line markers. PNG/SVG figures, an overview, calculated CSV and comparison reports are written under `results/paper_figures/`. See the [figure-by-figure method, findings and provenance](paper_reference/README.md). The compressor velocity-ratio correlation follows Figure 4 and its mean-line data; the production docstring documents the erroneous coefficients in printed Equation 23. Figure 3 is explicitly labeled as a perfect-gas, Mach-only reconstruction.

## Numerical implementation and JAX

`compressor_performance()` evaluates the arithmetic at a trial outlet without property calls. `turbine_performance()` additionally performs one rotor-static h,s flash through the supplied fluid. Both return the same eight performance quantities as the solve functions. Input states use the standard jaxprop names and can also be numerical mappings.

The `xp` keyword selects NumPy (default) or `jax.numpy`. This keeps host evaluations free of eager JAX dispatch overhead while allowing the same correlation equations to be compiled and differentiated:

```python
import jax
import jax.numpy as jnp
from thermopt.components.turbomachinery_nondimensional import compressor_performance

# Fixed trial thermodynamic states; this is a partial derivative of the
# correlation, not the total derivative of a solved design point.
outlet = result["state_out"]
isentropic = result["state_out_isentropic"]

def efficiency(flow_coefficient):
    """Evaluate the trial-state polytropic correlation at a chosen flow coefficient."""
    return compressor_performance(
        inlet, outlet, isentropic, 2.0,
        {"flow_coefficient": flow_coefficient}, xp=jnp,
    )["correlation_efficiency"]

value, derivative = jax.jit(jax.value_and_grad(efficiency))(0.07)
```

For traced turbine calculations, the fluid must itself support h,s under JAX; tests exercise the complete calculation with `FluidPerfectGas`. HEOS remains a host property backend. `solve_compressor()` and `solve_turbine()` still use SciPy Brent iteration and are **not differentiable/JIT-compatible**. Full solved-design derivatives require a differentiable root solver or implicit differentiation of the efficiency residual, plus an appropriate property model. No claim is made that a partial correlation derivative equals that total derivative. Gradients at piecewise correction thresholds need the usual care.

Outlet iteration always uses enthalpy. The logarithmic mean in Eq. (9) uses `log1p` and a continuous series limit at equal temperatures. The turbine entropy allocation is algebraically simplified without changing the model. Inlet states and isentropic flashes are reused; cycle integration preserves generalized-quality annotations using the solved enthalpy, without converting through temperature.

A bounded comparison made before the Equation 23 correction used the saved original classes and the same Air/HEOS fluid, included the inlet flash in both paths, and took the median of five groups of 150 warm solves. At the demonstration defaults, compressor time changed from 8.85 to 7.66 ms and turbine time from 7.61 to 8.04 ms. These shared-machine samples show a modest compressor improvement and no turbine speed gain from restructuring; they are not timings of the later model correction. Native property evaluation and a compiled outlet solve remain the opportunities for larger gains.

## Paper review and corrected defect

| Comparison | Finding |
|---|---|
| Section 3.2, Table 4 / Eqs. 12–27 | Compressor baseline and corrections retained, with the velocity-ratio correction below. |
| Section 3.2, Eq. 23 / Figure 4 | Use $w/U=0.38+3\phi-7.8\phi^2$, supported by Figure 4 and its mean-line markers. The printed $0.42+5\phi-14\phi^2$ expression is a paper error, documented in the production docstring. |
| Section 3.3, Eq. 29 / Figure 6 / Table 6 | **0.50 retained** to match the plot and author code. Table 6's **0.30** is a typo, documented beside the coefficient. |
| Turbine Eqs. 30–41, 46–48 | Scalloping, clearance, diffuser, Reynolds and rotor-state calculations retained. |
| Appendix B | Fixed missing shrouded-compressor stress keys; interpolate 0.67–0.75 between flow coefficients 0.05–0.12, clipping beyond them. |
| Section 4 | Retained scalar efficiency matching; outlet enthalpy replaces the temperature/enthalpy branch. |

The [author implementation](https://github.com/simone-parisi/NonDimensionalTurbo) is an additional reference, not an independent validation dataset.

## Tests and reference preservation

```powershell
.\.venv\Scripts\python.exe -m pytest tests/test_components.py tests/test_machine_property_access.py tests/test_turbomachinery_correlations.py tests/test_turbomachinery_functional.py
.\.venv\Scripts\python.exe -m pytest
```

Coverage includes independent compressor equations at reference size/speed, the turbine Figure 6 correlation away from default speed, all 85 frozen demo designs, four additional high-pressure Air/CO₂ designs, physical scale/loss trends, shrouded stress interpolation/endpoints, equal-temperature limits, solver failures, explicit-fluid access, and JAX JIT/gradient/batching checks against NumPy and finite differences. Cycle evaluations and complete optimization regressions exercise the integration.

`tests/data/turbomachinery_trends.json` was established before restructuring. `tests/data/turbomachinery_real_fluid.json` records four additional designs, which had also been checked using temperature and enthalpy roots. The functional refactor preserved those results. The later authorized correction to Equation 23 intentionally changes compressor Reynolds corrections, so affected compressor and cycle reference values are explicitly refreshed and documented in the test suite. Ordinary tests never update references. The recorder requires an explicit reviewed overwrite:

```powershell
.\.venv\Scripts\python.exe development/turbomachinery/record_trends.py --overwrite
```

Physical validation against raw experimental/mean-line data remains outstanding; those datasets were unavailable in the public example repository. Passing equation and regression tests is not experimental certification.
