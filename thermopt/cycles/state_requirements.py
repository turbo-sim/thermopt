"""Resolve cycle state requirements once, outside numerical evaluations.

Routes describe the actual stream supplied to each component, including PTES's
historical hot/cold storage naming. They are setup metadata, not surrogate fluid
objects or a second cycle implementation. Tests compare them with evaluated
states so a future topology change cannot silently select the wrong backend.
"""

import re

import jaxprop as jp

_STANDARD = frozenset(
    {
        "simple",
        "power_simple",
        "recuperated",
        "power_recuperated",
        "refrigeration_simple",
        "refrigeration_recuperated",
    }
)
_SPLIT = frozenset({"recompression", "split_compression", "power_split_compression"})
_STORAGE = frozenset({"PTES_recuperated", "PTES_recuperated_turbo"})
_DIAGNOSTICS = frozenset(
    {
        "subcooling",
        "superheating",
        "pressure_saturation",
        "temperature_saturation",
        "supersaturation_degree",
        "supersaturation_ratio",
    }
)
_STATE_REFERENCE = re.compile(
    r"\$components\.(\w+(?:\.(?:hot_side|cold_side))?)\."
    r"(?:state_in|state_out|states)\.(\w+)"
)


def cycle_state_routes(topology):
    """Return component paths mapped to (fluid role, diagnostic-enabled path).

    Machinery requests optional source diagnostics at component entry. Ordinary
    exchanger profiles do not enable them; neither route changes a model's own
    supported properties. A fresh dictionary is returned for setup callers.
    """
    working = "working_fluid"
    if topology in _STORAGE:
        heating, cooling = "hot_storage_fluid", "cold_storage_fluid"
        machines = {
            f"{kind}_{mode}": working
            for mode in ("charge", "discharge")
            for kind in ("compressor", "expander")
        }
        exchangers = {
            f"{name}_{mode}": roles
            for mode in ("charge", "discharge")
            for name, roles in (
                ("recuperator", (working, working)),
                ("heater", (heating, working)),
                ("cooler", (working, cooling)),
            )
        }
    elif topology in _STANDARD or topology in _SPLIT:
        heating, cooling = "heating_fluid", "cooling_fluid"
        names = (
            ("main_compressor", "split_compressor", "expander")
            if topology in _SPLIT
            else ("compressor", "expander")
        )
        machines = {name: working for name in names}
        machines.update(heat_source_pump=heating, heat_sink_pump=cooling)
        recuperators = (
            ("recuperator_lowT", "recuperator_highT")
            if topology in _SPLIT
            else ("recuperator",)
        )
        exchangers = {name: (working, working) for name in recuperators}
        exchangers.update(heater=(heating, working), cooler=(working, cooling))
    else:
        raise ValueError(f"Invalid cycle topology: {topology!r}")
    routes = {name: (role, True) for name, role in machines.items()}
    for name, (hot, cold) in exchangers.items():
        routes[f"{name}.hot_side"] = (hot, False)
        routes[f"{name}.cold_side"] = (cold, False)
    return routes


def _require_properties(fluid, properties, *, context):
    """Check declared support only; local state availability remains numerical."""
    for name in sorted(set(properties) - fluid.capabilities.properties):
        raise jp.UnsupportedPropertyError(
            fluid_name=fluid.name,
            provider=fluid.provider,
            model=fluid.model,
            request=name,
            reason=f"ThermOpt {context} requires this property.",
            supported=sorted(fluid.capabilities.properties),
        )


def validate_cycle_state_requirements(
    topology, parameters, fluids, constraints, objective
):
    """Validate model support and expression fields during problem configuration.

    Inspect constraint variables AND targets and the objective. Do not enable
    disabled diagnostics, substitute zeros or drop unsupported constraints.
    Capabilities cannot establish whether an enabled optional property is finite
    at a particular state: native state access must retain that typed failure.
    No parsing, model validation or new wrapper enters repeated fitness calls.
    """
    routes = cycle_state_routes(topology)
    for role, fluid in fluids.items():
        _require_properties(
            fluid, {"pressure", "temperature", "enthalpy"}, context=role
        )
        required_pairs = {jp.HmassP_INPUTS: "HmassP_INPUTS", jp.PT_INPUTS: "PT_INPUTS"}
        machines = [name for name, route in routes.items() if route == (role, True)]
        if machines:
            _require_properties(fluid, {"entropy", "density"}, context=role)
            required_pairs[jp.PSmass_INPUTS] = "PSmass_INPUTS"
        for name in machines:
            if parameters[name].get("efficiency_type") == "non-dimensional":
                properties = {"viscosity"}
                if name.startswith("expander"):
                    required_pairs[jp.HmassSmass_INPUTS] = "HmassSmass_INPUTS"
                else:
                    properties.add("speed_of_sound")
                _require_properties(
                    fluid, properties, context=f"{role} radial machinery"
                )
        for pair, name in required_pairs.items():
            if pair not in fluid.capabilities.input_pairs:
                raise jp.UnsupportedInputPairError(
                    fluid_name=fluid.name,
                    provider=fluid.provider,
                    model=fluid.model,
                    request=name,
                    reason=f"ThermOpt {role} requires this state input pair.",
                    supported=tuple(sorted(fluid.capabilities.input_pairs)),
                )

    expressions = [objective["variable"]]
    for constraint in constraints:
        expressions.extend((constraint["variable"], constraint["value"]))
    for expression in expressions:
        if not isinstance(expression, str):
            continue
        for path, field in _STATE_REFERENCE.findall(expression):
            if path not in routes:
                raise jp.InvalidConfigurationError(
                    request=f"components.{path}",
                    reason=f"No state stream with this component path exists in {topology}.",
                    supported=sorted(routes),
                )
            if field in ("status", "valid", "clipped"):
                continue  # Public status fields are not thermodynamic properties.
            canonical = jp.ALIAS_TO_CANONICAL.get(field)
            if canonical is None:
                raise jp.InvalidConfigurationError(
                    request=field,
                    reason=f"Unknown fluid-state property in {expression!r}.",
                    supported=jp.PROPERTIES_CANONICAL,
                )
            role, diagnostics = routes[path]
            fluid = fluids[role]
            _require_properties(
                fluid, {canonical}, context=f"expression {expression!r} ({role})"
            )
            if canonical in _DIAGNOSTICS and (
                not diagnostics or "supersaturation" not in fluid.capabilities.options
            ):
                raise jp.UnsupportedPropertyError(
                    fluid_name=fluid.name,
                    provider=fluid.provider,
                    model=fluid.model,
                    request=canonical,
                    reason=f"Supported by the model but not enabled on component path {path!r}; expression {expression!r} cannot use it.",
                    supported=sorted(fluid.capabilities.properties - _DIAGNOSTICS),
                )
