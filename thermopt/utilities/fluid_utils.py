"""Shared property models for study scripts, owned entirely by jaxprop."""

from collections.abc import Mapping

import jaxprop
import numpy as np

__all__ = [
    "fluid_from_configuration",
    "get_cycle_fluid",
    "process_state_options",
]


def _scalar_settings(value):
    """Copy static configuration mappings into native Python scalar metadata.

    YAML processing uses NumPy scalars, while read-only mappings are convenient
    for sharing configuration safely. Normalize both once at setup; numerical
    arrays and configured model objects retain their original meaning.
    """
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Mapping):
        return {key: _scalar_settings(item) for key, item in value.items()}
    return value


def fluid_from_configuration(configuration, *, identifier=None):
    """Construct a native evaluator using the JaxProp 0.6 vocabulary directly.

    The configuration is copied, never mutated. JaxProp validates its own
    selectors/options and reports unsupported constructor fields. In particular,
    error policies belong on ``get_state`` calls, not in model construction.
    No ThermOpt compatibility translation or property-call wrapper is applied.
    """
    settings = _scalar_settings(dict(configuration))
    if identifier is not None:
        settings["identifier"] = identifier
    return jaxprop.Fluid(**settings)


def get_cycle_fluid(parameters, role, fluids=None):
    """Consume a stream configuration and reuse its problem-owned evaluator.

    Standalone cycle calls can omit ``fluids`` and construct their own models.
    Optimization problems supply preconfigured objects so table loading and
    metadata validation occur during setup rather than each cycle evaluation.
    """
    settings = parameters.pop(role)
    return (
        fluids[role]
        if fluids is not None
        else fluid_from_configuration(settings, identifier=role)
    )


def process_state_options(fluid, *, generalize_quality=False):
    """Select optional source diagnostics once at component entry.

    Work/energy equations need ordinary complete states. Preserve source-EOS
    reporting where supported, without asking bicubic/perfect-gas models for
    saturation diagnostics or inventing replacement properties for them.
    """
    supported = fluid.capabilities.options
    options = {}
    if "supersaturation" in supported:
        options["supersaturation"] = True
    if generalize_quality and "generalize_quality" in supported:
        options["generalize_quality"] = True
    return options
