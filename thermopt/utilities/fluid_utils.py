"""Shared property models for study scripts, owned entirely by jaxprop."""

from functools import lru_cache

import jaxprop
import numpy as np

__all__ = [
    "get_property_fluid",
    "fluid_from_configuration",
    "get_cycle_fluid",
    "process_state_options",
]


def fluid_from_configuration(configuration, *, identifier=None):
    """Construct an actual JaxProp evaluator from a cycle's YAML settings.

    Old case files spell the default raising policy as ``exceptions: True``.
    Accept that setup spelling only; ThermOpt needs finite states and does not
    silently turn a disabled exception flag into masked or clipped calculations.
    New cases select ``method`` and pass ``table_options`` directly to JaxProp.
    """

    def scalar_settings(value):
        # YAML processing represents numbers as NumPy scalars. Model setup
        # options are static Python metadata, not numerical array parameters.
        # Normalize them once to avoid Equinox's static-array warnings and
        # keep every repeated property call free of this conversion work.
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, dict):
            return {key: scalar_settings(item) for key, item in value.items()}
        return value

    settings = scalar_settings(dict(configuration))
    exceptions = settings.pop("exceptions", True)
    if exceptions is not True:
        raise ValueError(
            "ThermOpt requires raising fluid errors; remove exceptions: False. "
            "JaxProp's mask/clip policies must be requested explicitly on state calls."
        )
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


@lru_cache(maxsize=32)
def get_property_fluid(name, backend="HEOS"):
    """Reuse a jaxprop model within a process; obtain states with get_state."""
    return jaxprop.Fluid(name=name, backend=backend)
