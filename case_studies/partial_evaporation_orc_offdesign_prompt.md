# Implementation prompt: design and off-design partial-evaporation ORC

Implement and validate a reproducible design and off-design model of a partial-evaporation organic Rankine cycle (ORC). Complete the implementation and demonstrate it with a runnable example; do not stop after proposing a plan.

## 1. Workspace and scope

Locate these repositories on this computer. The paths below are examples from the original computer, not paths to hard-code:

- `C:/Users/roagr/git_projects/thermopt`
- `C:/Users/roagr/git_projects/flashing_nozzle_design`
- Optionally, a compatible sibling `pysolver_view` checkout if the installed package lacks the required solver capabilities.

Put **all new code, configuration, documentation, tests, environment instructions, and generated results** inside `thermopt/case_studies/partial_evaporation_orc_offdesign/`. Do not edit `thermopt/thermopt/`, existing examples, repository dependency manifests or lockfiles, or the nozzle and solver repositories. Do not monkey-patch their functions or global registries. Use imports, composition, and case-local adapters or subclasses where needed. Preserve unrelated working-tree changes. Do not commit or push unless requested.

Read applicable `AGENTS.md` files. Use `$...$` for inline Markdown math and `$$...$$` for display math.

Inspect the current code before choosing APIs. Useful reference locations, verified in the original checkouts, are:

| Purpose | Reference |
| --- | --- |
| Native partial-evaporation design configuration | `thermopt/examples/rankine_organic/case_butane_PEORC.yaml` |
| Native design optimizer and problem interface | `thermopt/thermopt/optimize_cycle.py` |
| Simple-cycle evaluator | `thermopt/thermopt/cycles/cycle_power_simple.py`, delegating to `cycle_power_recuperated.py` with recuperation disabled |
| Component and HX discretization functions | `thermopt/thermopt/components/basic_components.py` |
| Controlled SLSQP integration | `thermopt/examples/sCO2_recuperated/compare_slsqp.py` |
| Choked-flow calculation | `flashing_nozzle_design/src/flashing_nozzle_design/functions.py`, function `find_choked_state` |

The original repositories declared incompatible JaxProp requirements: ThermoOpt required `^0.5.7`, while the nozzle repository pinned `0.5.6`. Resolve imports and verify the needed nozzle function in a compatible, reproducible case environment. Do not blindly install both editable projects with dependency resolution or silently downgrade ThermoOpt's dependencies. A documented case-local source-path adapter is acceptable after compatibility checks. Record actual import locations, package versions, repository commits, and dirty-tree status. Verify solver capabilities rather than relying on a version label alone.

## 2. Configuration and cycle definition

Provide one clearly organized master YAML containing all independent model inputs: fluids and property backends; nominal and off-design boundary conditions; design and operating bounds; component parameters; closure policy; objectives; solver options; initial guesses; sweep definitions; output settings; and validation tolerances. Use SI units internally, explain units in comments, validate inputs, and reject misspelled or unused keys. Derived quantities belong in a saved design artifact rather than being duplicated as independent inputs.

Use a simple, steady-state loop:

1. Condenser outlet / working-fluid pump inlet: state 1.
2. Pump outlet / evaporator inlet: state 2.
3. Evaporator outlet / turbine inlet: state 3.
4. Turbine outlet / condenser inlet: state 4.

Use ThermoOpt property and component functions wherever appropriate. The nozzle repository supplies choking physics; do not substitute its complete test-rig cycle, which contains additional assumptions and can include an upstream throttle. Use consistent equilibrium two-phase properties throughout the baseline model and document their limitations for flashing flow.

The prescribed operating boundary conditions are heat-source mass flow $\dot m_h$, heat-source inlet temperature $T_{h,in}$, and heat-sink supply temperature $T_{k,in}$. Fluid identities and supply pressures are fixed configuration inputs. Define precisely where temperatures are prescribed relative to auxiliary pumps, follow the same convention at design and off design, and include pump heating consistently.

The adjustable controls are heat-sink mass flow $\dot m_k$ and working-fluid pump angular speed $\omega_p$. Turbine blade speed and nozzle throat area remain fixed after sizing. Working-fluid mass flow, cycle pressures, outlet temperatures, and turbine inlet quality must follow from the coupled model and operating policy.

Use configurable pressure-loss assumptions consistently at design and off design. The baseline may retain ThermoOpt's fractional pressure-loss convention; document that approximation. Include source/sink pumping power consistently when pressure losses require it. Impose finite control bounds so unlimited cooling flow cannot be an artificial optimum.

Provide a runnable butane example informed by the existing PEORC configuration. Clearly label any fluid, nominal boundary value, efficiency, or operating bound you select as an example assumption rather than user-specified data.

## 3. Optimize and size the design point

Perform a real design optimization using ThermoOpt's functionality, with user-configurable design turbine and pump efficiencies. Optimize cycle states and heat recovery subject to configured pinch, pressure, temperature, subcooling, and phase constraints. For the baseline partial-evaporation mode, require a genuinely two-phase turbine inlet using configurable quality bounds strictly between zero and one; do not silently return a conventional superheated ORC. Do not fix inlet quality to an arbitrary value.

Support objectives `net_power` and `thermal_efficiency`, with `net_power` as the default at the specified heat-source flow. Define net power as turbine power minus working-fluid pump power and included auxiliaries. Define thermal efficiency as net power divided by actual evaporator heat input. Report the heat-recovery fraction separately. ThermoOpt's `system_efficiency` uses a maximum available source-heat denominator; label that metric separately if reported.

**Account for ThermoOpt's native size convention.** The inspected native cycle computes working-fluid mass flow from a prescribed `net_power`, then derives source and sink flows. That is a reference size, not the off-design closure. For design at fixed source flow, a suitable native power objective is net system power per source mass flow. Optimize the intensive design, then scale all extensive flows, duties, and powers to the configured design source flow before sizing hardware. Verify the assumptions permitting this scaling. Do not accidentally optimize a fixed reference power or carry the reference-power equation into off design.

Use the same closure policy in design and off-design recovery. Save a machine-readable design artifact containing states, flows, duties, efficiencies, HX conductances, throat sizing, fixed blade speed, pump reference quantities, resolved configuration, provenance, and units. Subsequent off-design runs must load this hardware without redesigning it.

## 4. Heat exchangers: design UA and off-design scaling

Model counter-current evaporator and condenser using ThermoOpt's discretized property profiles, with **30 physical elements by default**, configurable in YAML. In the inspected implementation, `heat_exchanger(..., num_steps=...)` creates that many **nodes**, and the cycle's `num_elements` is forwarded as `num_steps`. Check the local version; use 31 nodes for 30 intervals when that behavior is unchanged. Explain the mapping in the case adapter.

Compute design conductance by integrating along each exchanger. For element $i$:

$$
\Delta T_{lm,i}=\frac{\Delta T_i-\Delta T_{i+1}}{\ln(\Delta T_i/\Delta T_{i+1})},
\qquad
(UA)_{req}=\sum_i\frac{\Delta\dot Q_i}{\Delta T_{lm,i}}.
$$

Use the continuous equal-temperature-difference limit, correct counter-current alignment, and positive heat-duty increments. Check element duties against both stream energy balances. Reject temperature crosses and nonpositive approaches; do not hide them with absolute values or artificial temperature floors. A single whole-exchanger LMTD is insufficient when phase change and temperature-dependent properties occur.

Reuse saturation-node support if appropriate, recognizing that it moves nodes and can make duty spacing nonuniform. Always compute $\Delta\dot Q_i$ from the actual aligned enthalpy/duty coordinates after node adjustment. Check discretization accuracy around phase boundaries.

At off design impose:

$$
(UA)_{evap,req}=(UA)_{evap,d}
\left(\frac{\dot m_h}{\dot m_{h,d}}\right)^{n_e},
$$

$$
(UA)_{cond,req}=(UA)_{cond,d}
\left(\frac{\dot m_k}{\dot m_{k,d}}\right)^{n_c}.
$$

Both exponents default to **0.8** and are user-configurable. The scaling flows are explicitly the source stream for the evaporator and sink stream for the condenser, reflecting dominant external-side thermal resistance. Do not use working-fluid flow for these scaling laws.

Area is fixed. Since a thermal design point determines $UA$, not $U$ and $A$ independently, use fixed design $UA$ plus these multipliers as the baseline. If physical areas are requested, require configured design $U$ values or areas and infer the complementary quantity. Do not claim an independently identified area from $UA$ alone. With uniform effective $U$ per exchanger, close the model using the total required conductance; uniform-duty integration cells are not fixed physical area cells, and their individual design $UA_i$ values should not be frozen as phase fronts move.

## 5. Turbine and fixed choked throat

Use `flashing_nozzle_design.functions.find_choked_state(fluid, h0, s0, T_min)` or a thin case-local adapter. Inspect the implementation: in the checked version it returns **`(state_star, G_star, v_star)`**, despite an inconsistent return-value docstring. It maximizes equilibrium isentropic mass flux. Validate that the computed maximum is interior to an adequate property-domain search interval rather than an artificial temperature or velocity bound.

Treat the cycle inlet as a stagnation state under the negligible-port-kinetic-energy approximation. Size total effective throat area from the final scaled design:

$$
A_t=\frac{\dot m_{wf,d}}{G_d^*},
\qquad
\dot m_{wf}=A_tG^*(h_3,s_3)
\quad\text{at off design}.
$$

Keep $A_t$ fixed. If an optional discharge coefficient is introduced, configure it explicitly and apply it consistently in both equations. Do not apply turbine efficiency to the ideal choking mass flux.

Assume choked operation, but verify compatibility with downstream pressure and flag infeasible points. For a simplified converging-throat representation, use the corresponding critical-pressure criterion. Document that a real converging-diverging nozzle or turbine can require a different backpressure criterion; do not claim geometry-resolved choking from this reduced model.

For the turbine efficiency map:

$$
\Delta h_s=h_3-h(p_4,s_3),
\qquad c_0=\sqrt{2\Delta h_s},
\qquad U=0.5\sqrt{2\Delta h_{s,d}},
$$

$$
\nu=\frac{U}{c_0},
\qquad \eta_t=4\eta_{t,d}\nu(1-\nu),
\qquad h_4=h_3-\eta_t\Delta h_s.
$$

Here $U$ is constant after sizing. This puts the design at the impulse-turbine peak $\nu_d=0.5$, reproduces the configured design efficiency, and gives zero efficiency at $\nu=0$ and $\nu=1$. Require $\Delta h_s>0$ and an admissible velocity-ratio range. Do not clip invalid efficiencies into apparently feasible results. Treat this as a preliminary performance correlation, not a calibrated wet-turbine map.

## 6. Variable-speed pump with consistent endpoint behavior

Use a centrifugal-pump similarity model based on **inlet volumetric flow**, not uncorrected mass flow. Let

$$
Q_v=\frac{\dot m_{wf}}{\rho_1},
\qquad r=\frac{\omega_p}{\omega_{p,d}},
\qquad \phi=\frac{Q_v}{rQ_{v,d}},
\qquad k=\frac{H_{shutoff,d}}{H_d}>1,
$$

where $k$ defaults to **1.20**, and $\omega_{p,d}$ is configurable. Use

$$
H(Q_v,r)=H_dr^2\left[k-(k-1)\phi^2\right],
\qquad \phi_0=\sqrt{\frac{k}{k-1}}.
$$

This gives the desired zero slope at zero flow, design head at $\phi=1$, and shutoff head $kH_d$ at design speed. Use the affinity laws $Q_v\sim\omega_p$, $H\sim\omega_p^2$, and shaft power $\sim\rho\omega_p^3$ at corresponding operating points, with invariant efficiency versus $\phi$ as the baseline approximation.

**Explain and resolve the incompatible original parabola requirements.** A symmetric efficiency parabola with its peak at $\phi=1$ and one zero at zero flow has its other zero at $\phi=2$. The head curve above reaches zero there only if $k=4/3$. With $k=1.20$, its zero-head point is $\phi_0=\sqrt{6}\simeq2.449$, while head is still $0.4H_d$ at $\phi=2$ and design speed. Combining those two parabolas would make calculated shaft power diverge. Do not silently impose all these inconsistent assumptions.

Use the following proposed engineering surrogate as the default, retaining the requested head curve and relaxing efficiency symmetry:

$$
\eta_p(\phi)=\eta_{p,d}
\frac{\phi\left[k-(k-1)\phi^2\right]}
{(k-\tfrac12)+(\tfrac32-k)\phi^2},
\qquad 0<\phi<\phi_0.
$$

It has its maximum $\eta_{p,d}$ at $\phi=1$, vanishes at zero flow and the zero-head endpoint, and yields finite shaft-power limits. It is an analytically constructed surrogate, **not** a published empirical pump correlation or manufacturer map. Derive and verify these properties. If strict symmetric parabolas are selected as an alternative, require $k=4/3$ and reject inconsistent configurations.

For compatibility with ThermoOpt's isentropic compression model, define pump head consistently as specific isentropic work divided by $g$:

$$
H_{req}=\frac{h(p_2,s_1)-h_1}{g},
\qquad H_{req}=H(Q_v,r),
\qquad h_2=h_1+\frac{gH_{req}}{\eta_p}.
$$

For a nearly incompressible liquid, $H_{req}\simeq(p_2-p_1)/(\rho_1g)$, recovering the requested pressure-rise curve. State the efficiency convention explicitly; do not mix pressure-head hydraulic efficiency with isentropic efficiency without accounting for the approximation.

Compute pump power and torque consistently. Under the head convention above, the cancellation gives a numerically useful map expression:

$$
P_p=\dot m_{wf}(h_2-h_1)
=\frac{\rho_1gQ_{v,d}H_dr^3}{\eta_{p,d}}
\left[(k-\tfrac12)+(\tfrac32-k)\phi^2\right],
\qquad \tau_p=\frac{P_p}{\omega_p}.
$$

At $k=1.20$, the bracket approaches 0.7 at shutoff and 2.5 at the zero-head endpoint, so power and torque remain finite at fixed nonzero speed. Do not evaluate the cycle exactly at zero flow or zero efficiency: finite shaft power does not remove the singular specific heating per unit mass near shutoff. Configure practical positive flow, speed, efficiency, and liquid-inlet/subcooling limits. Endpoint coincidence alone does not prove sensible torque; verify the ratio and its limits. Do not claim cavitation prediction without an NPSH model and its required data.

## 7. Equation closure and operating optimization

Before implementing the solver, write an explicit unknown/equation table and check the numerical rank of the scaled equality Jacobian. Distinguish physical equalities, definitions eliminated by component evaluations, operating policies, inequalities, and adjustable controls. Do not add an overall energy equality if it is already implied by component balances.

The boundary conditions alone do not prescribe liquid inventory or condenser-outlet subcooling. Make this modeling choice explicit in YAML. Use **fixed pump-inlet subcooling** as the baseline closure, with a configurable example value of 1 K, applied consistently at design and off design. Identify this as an idealized operating policy, not a derived hardware law or a rigorous inventory model. Provide an optional minimum-subcooling-only mode only if its additional degree of freedom is handled explicitly.

For illustration, with pressure relations and component enthalpy relations evaluated internally, a candidate unknown vector is

$$
x=[p_3,p_1,h_3,h_1,\dot m_{wf},T_{h,out},T_{k,out},\dot m_k,\omega_p].
$$

The baseline has seven candidate equalities: two HX energy balances, two HX conductance matches, one choking mass-flow condition, one pump-head condition, and one fixed-subcooling condition. If independent, this leaves **two operating degrees of freedom**. Fixing both controls can then give a square system. Removing fixed subcooling without another closure introduces another degree of freedom. Verify this count for the actual implementation and account for any additional physics or algebraic eliminations; do not present the illustrative count as proof of full rank.

Support:

- **Fixed-control simulation:** prescribe $\dot m_k$ and $\omega_p$, eliminate fixed variables, and solve the remaining physical equalities as a feasibility problem.
- **Optimal operation:** keep those controls free within bounds and maximize the configured net-power or thermal-efficiency objective subject to the same equations and physical inequalities.

Do not overconstrain the model by simultaneously prescribing source flow, design net power, working-fluid flow, outlet temperatures, and both controls. If closure is inconsistent or rank deficient, diagnose the particular equations or missing assumptions and correct the case formulation explicitly.

## 8. Required solver and numerical handling

Use `pysolver_view` with **`library="pysolver"`, `method="slsqp"`** for both design optimization and off-design feasibility/optimization. Calling SciPy SLSQP through its standard backend does not provide the required first-trial control. A scalar maximization inside the imported choking routine is acceptable as a component calculation; the coupled cycle solves must use the requested backend.

Expose the controlled line-search parameters in YAML, following the checked ThermoOpt comparison example:

```yaml
solver_options:
  library: pysolver
  method: slsqp
  derivative_method: 2-point
  derivative_abs_step: 1.0e-6
  tolerance: 1.0e-8
  max_iterations: 300
  problem_scale: 1.0
  extra_options:
    first_step_fraction: 0.2
    max_step_fraction: 1.0
    step_growth_factor: 1.5
    domain_backtracking: true
    feasibility_tol: 1.0e-6
    optimality_tol: 1.0e-4
```

These are starting numerical settings, not guaranteed convergence parameters. Check support in the local package and record the effective options. The inspected custom solver bridge requires SciPy **1.17.1**; verify the requirement in the version actually used. Preserve integer types for iteration limits when ThermoOpt's YAML conversion produces NumPy scalars. Distinguish `first_step_fraction` from the finite-difference derivative step; prove that it reaches the backend rather than merely appearing in the YAML.

Use well-scaled variables, residuals, and objectives. Document manual scaling versus `problem_scale` to avoid scaling twice. Verify inequality sign conventions: the inspected public `OptimizationProblem` interface uses inequalities $c\leq0$, while the internal SLSQP driver uses nonnegative slacks and the package adapter converts them.

For square systems, minimize a constant zero objective with the physical residuals as equality constraints, using a case-local `OptimizationProblem` or a compatible root-to-optimization adapter. Independently verify the original residuals after solving. Do not accept nonzero residuals just because a least-squares objective has become stationary. Eliminate fixed variables rather than passing equal lower and upper bounds to the custom backend; remove redundant equalities.

Handle recognized thermodynamic-domain failures through `pysolver_view.DomainError` and supported backtracking. Do not catch programming errors indiscriminately, fabricate thermodynamic states, or return fake successful penalties. Note that domain recovery for line-search trials does not guarantee valid finite-difference probes. Use physically informed initial guesses, design-point recovery, and continuation/warm starts for sweeps. Report failed or infeasible conditions honestly.

## 9. Deliverables and validation

Keep the implementation modular and reasonably small. Provide a master YAML, configuration loader/validator, design-sizing code, component adapters, coupled off-design problem, CLI runner, focused validation, and a README. Save generated artifacts under a case-local `results/` directory. Avoid a notebook-only workflow or hard-coded computer-specific paths.

Document exact environment setup and commands for:

1. Optimizing and sizing the design from YAML.
2. Reloading the saved hardware and recovering the design point with fixed controls.
3. Simulating a specified off-design condition with fixed controls.
4. Optimizing controls at a specified boundary condition.
5. Running small configured sweeps of source flow, source temperature, and sink temperature without resizing hardware.

Validate at least:

- Design recovery reproduces states, flows, duties, pump/turbine efficiencies, and net power within configured tolerances.
- Both HX stream balances, cycle energy balance, conductance residuals, pump-head residual, nozzle mass-flow residual, and subcooling closure pass independent post-solve checks.
- A direct check confirms area/design conductances, throat area, blade speed, and pump reference geometry/map do not change across off-design points.
- Pump head, efficiency peak, shared high-flow zero, finite power/torque limits, and affinity scaling behave as derived. Turbine efficiency has the specified peak and zeros.
- The 30-element baseline is checked against finer reevaluation, for example 60 and 120 elements, for UA and minimum approach at representative design/off-design states. Report any material discretization error rather than concealing it.
- At least one cooler/warmer source, lower/higher source flow, and cooler/warmer sink case is attempted. Separate infeasibility from solver failure and retain useful diagnostics.
- Optimization acceptance includes physical feasibility and appropriate first-order optimality checks, with effective backend settings and convergence histories saved. Do not claim a global optimum from one SLSQP run.
- A final diff or source fingerprint confirms protected source files are unchanged.

Export state tables, operating/control variables, heat duties, gross and net powers, efficiencies with explicit denominators, minimum approaches, required/available UA, choking diagnostics, residuals, constraint violations, solver status, and provenance to readable machine formats such as YAML/JSON and CSV. Include clear T-s and T-Q plots, pump maps with operating points, turbine efficiency versus velocity ratio, and small off-design performance plots. Failed points must remain marked as failed rather than being interpolated into reported solutions.

In the final implementation report, explain the adopted closure and verified degrees of freedom, the pump-model correction, dependency resolution, commands executed, validation results, remaining approximations, and files created. Complete the work autonomously within this scope; ask for clarification only when missing information prevents a defensible implementation.

## Background for the pump modeling assumptions

The speed scaling follows the [Hydraulic Institute's pump curves and affinity rules](https://datatool.pumps.org/pump-fundamentals/pump-curves.html) and [KSB's affinity-law definitions](https://www.ksb.com/en-global/centrifugal-pump-lexicon/article/affinity-laws-1115682). These references support the similarity assumptions and distinction between head, efficiency, and shaft-power curves. The rational efficiency function above is a proposed analytical surrogate derived to satisfy this task's constraints; it is not taken from those references. Real pump operating limits and maps require manufacturer or experimental data.
