# Repository instructions

## Architecture direction (user preference, 15 November 2026)

ThermOpt will be ported to JAX in a later, separate sweep. Keep that direction
in mind when defining functions, data structures and module boundaries now.
The current migration remains eager NumPy/SciPy; do not introduce a cycle JAX
port or replace host solvers without a new instruction.

- Prefer small, composable numerical functions with explicit inputs and outputs.
  Avoid new mutable numerical objects, hidden global state, and calculations
  that depend on changing `self` attributes. Existing configured JaxProp model
  objects are valid inputs; this preference is not a ban on their use.
- Keep configuration parsing, fluid construction, reference preparation, solver
  orchestration, logging, plotting and file output outside numerical functions.
  Place validation at setup boundaries whenever possible, with no repeated
  metadata/capability checks added to thermodynamic state calculations.
- Keep numerical result structures predictable. Consider immutable records or
  dataclasses with a fixed schema when they clarify the model, but do not add
  them solely in anticipation of JAX. A frozen dataclass is not automatically
  a JAX PyTree: future transformation support must distinguish dynamic numeric
  leaves from static metadata and provide deliberate registration when needed.
- Treat static metadata equality/hash behavior as a measured performance concern.
  Equivalent JaxProp bicubic model instances can incur deep metadata equality
  during JIT-cache lookup. Counterbalance cache/model identity when benchmarking;
  do not add global sharing of mutable evaluators to disguise that cost.
- Preserve physical equations, failure semantics, useful docstrings, units and
  derivations. Future transformability does not justify weakening correctness
  checks or duplicating JaxProp thermodynamic models in ThermOpt.
- Follow the phase checklist in `z_JAXPROP_060_IMPLEMENTATION_PLAN.md`. Verify
  changes with meaningful tests and complete-consumer benchmarks. Record timing
  flags and numerical/model changes honestly; do not relax gates silently.
- Group generated example/case-study outputs under that example's `results/`
  directory, including figures, histories, table caches and diagnostic logs.
  Use subdirectories there instead of sibling `results_*` directories or files.
- Keep maintained runnable cases under `case_studies/`. Update active entry
  points and documentation when moving them; preserve dated benchmark and
  physical-reference records as historical evidence.
- Standard case YAMLs use PySolverView controlled SLSQP (`library: pysolver`,
  `method: slsqp`) with explicit initial-step controls; see `case_studies/README.md`.
  Keep intentionally named alternative-solver configurations and comparisons explicit.

## Markdown mathematics

Use `$content$` for inline LaTeX and `$$content$$` for displayed mathematics.
Do not use `\( ... \)` or `\[ ... \]` delimiters.
