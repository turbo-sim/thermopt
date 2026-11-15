# State access and capability inventory (Phase 3)

ThermOpt consumes JaxProp 0.6 evaluator objects and immutable `FluidState` values.
Numerical cycle/component functions remain eager NumPy/SciPy. The future JAX
direction in [AGENTS.md](../AGENTS.md) favors explicit numerical inputs/outputs;
this phase adds setup functions and metadata, without a new numerical class or
dataclass hierarchy.

## Consumers and requirements

| Consumer | Required canonical fields or operations | Optional/model-specific behavior |
|---|---|---|
| Heat-transfer profiles, exchanger duties and cycle energy accounting | `enthalpy`, `pressure`, `temperature`; HP/PT routes where used | Saturation-node snapping explicitly requires saturation support. Plain profiles do not enable source diagnostics. |
| Isentropic machinery | HP/PS, `enthalpy`, `entropy`, `pressure` | Source diagnostics are enabled only when the supplied model supports them. |
| Polytropic machinery | HP, `density`, plus the isentropic reference state | Scalar SciPy RHS calls and original solver tolerances are preserved. No capability checks enter the RHS. |
| Radial compressor | Above fields plus `viscosity`, `speed_of_sound`, `temperature` | Numerical mappings use canonical keys. Design validation stays outside correlation arithmetic. |
| Radial turbine | `viscosity`, `density`, HS rotor-static flash and ordinary machinery fields | All flashes use the supplied model; PG/JAX transformation regressions remain active. |
| Constraints and objective | Canonical state fields or JaxProp's public aliases in configured expressions | Known component-state references in variables, targets and objective are checked at setup. Unsupported constraints remain errors. Arbitrary expression semantics still belong to the existing resolver. |
| Phase-diagram backgrounds and initialization | Explicit source reference/phase operations | These queries remain separate from a table's evaluation rectangle. Additional lifecycle/reference work is Phase 5. |
| Cycle/process/pinch plots | `temperature`, `entropy`, `enthalpy`, `pressure`, `density` and configured axes | Plot aliases normalize at the plotting boundary. Canonical axes preserve auxiliary streams and existing display labels. |
| Excel exporter | Canonical property keys, existing friendly columns/units | Unsupported, disabled or locally unavailable optional values remain blank; no extra flashes or fabricated zeros. |
| Generic state serialization | Public `to_data()` / `from_data()` | Preserves provenance, enabled fields, status and masks. `to_dict()` alone is a property projection and is not a lossless state archive. |
| Current Rankine study evaluators | Public source states and configured state constraints | Public aliases remain valid input notation. Full standalone/example modernization and cache identities remain Phase 7. |

`cycles/state_requirements.py` describes component-to-stream ownership for every
public topology and alias. Tests compare these declarations to all 14 actual
cycle cases, including PTES's historical storage-role naming. Required input
pairs/fields and component-state expression requirements are checked only in
`update_problem`, after model construction and before reference preparation.
Local thermodynamic finiteness cannot be inferred from capabilities at setup.
Repeated fitness calls do not invoke this validator.

## Availability and quality

- **Unsupported by a model:** absent from `fluid.capabilities.properties`.
  Configuration requiring it raises `UnsupportedPropertyError` with the fluid,
  provider, canonical field and expression context.
- **Supported but disabled on a returned path:** absent from that state's
  enabled fields. Exchanger profiles do not request source supersaturation
  diagnostics; a constraint requiring them is rejected explicitly, without
  silently adding diagnostic flashes to every node.
- **Enabled but locally undefined:** an optional value may be NaN while core
  properties and status remain valid. Native guarded access raises
  `PropertyUnavailableError` under the normal raising policy. Expression
  resolution preserves this type and indices and adds expression context.
- **Failed/masked state:** status and NaNs survive serialization. Normal cycle
  calls retain the raising policy; this phase introduces no automatic clipping,
  masked optimization trials or per-property validity wrapper.

Physical quality is meaningful inside an actual two-phase equilibrium region.
Source generalized quality outside that region is an entropy indicator and may
exceed one; it is not a vapor mass fraction. Trimmed source quality outside the
dome is a 0/1 entropy label. The PG model does not support quality. A table's
declared support/local values belong to that table/model; support alone never
establishes that the current state is two phase. Tests retain source-specific
generalized behavior, including an enabled but undefined volume-quality field.

## Public state operations and remaining scope

Component paths now consistently contain a native batched `FluidState`, including
nondimensional machines; endpoints remain scalar. Public addition/stacking,
`at_index` and `flipped` preserve identity and status. No state private storage,
dataclass mutation or `states_to_dict` helper is used in production.

The heat-transfer blending example is the one active caller of the obsolete
collection helper found in this scan. It now uses `FluidState.stack`, native
construction/constants and an explicit smoothing-function import. A headless
test evaluates all 300 points and writes both figures to an isolated directory.
Its constant 500/5000 heat-transfer coefficients remain illustrative assumptions,
not a validated heat-transfer correlation.

Historical archives and isolated JAX prototypes retain their provenance. Phase 7
still owns the comprehensive standalone/example execution inventory. Phase 4
owns broader component/domain checks; Phase 6 owns additional output round trips
and any explicit opt-in exchanger diagnostic design. None is inferred complete
from this state-contract audit.
