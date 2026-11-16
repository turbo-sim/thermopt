"""Test-only control of an external backend's process-global energy convention.

JaxProp has no public setter for this external setting. Import only that setter;
all thermodynamic evaluation and reference inspection in the tests use JaxProp.
The property-boundary audit permits this exact import in this file only.
"""

from CoolProp.CoolProp import set_reference_state


def set_energy_reference(name, convention):
    """Change external source configuration; callers must restore it finally."""
    set_reference_state(name, convention)
