"""Study-file boundaries and source identity, never numerical-loop helpers."""

import hashlib
from importlib.metadata import version
from pathlib import Path

import thermopt as th
from jaxprop.bicubic.cache_metadata import metadata_fingerprint, source_signature


def read_formulation(filename):
    """Read a formulation export or select it from a complete manager snapshot.

    This only selects the documented YAML payload. Fluid recipes still pass
    unchanged to ThermOpt/JaxProp; obsolete constructor options are not adapted.
    """
    data = th.read_configuration_file(filename)
    return data.get("problem_formulation", data)


def calculation_identity(config, study_sources):
    """Fingerprint study/core code, dependencies and active EOS references.

    Call at setup/cache boundaries. This deliberately does not memoize global
    source/reference state or add provenance probes to property calculations.
    Old saved results remain historical and are rejected after a real change.
    """
    package = Path(th.__file__).resolve().parent
    core = {
        str(path.relative_to(package)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(package.rglob("*.py"))
    }
    sources = {
        Path(path).name: hashlib.sha256(Path(path).read_bytes()).hexdigest()
        for path in study_sources
    }
    fluids = {}
    for role, settings in config["problem_formulation"]["fixed_parameters"].items():
        if role.endswith("fluid") and isinstance(settings, dict) and "name" in settings:
            fluids[role] = metadata_fingerprint(
                source_signature(settings["name"], settings.get("backend", "HEOS"))
            )
    return dict(
        schema=1,
        core_sources=core,
        study_sources=sources,
        fluid_sources=fluids,
        dependencies={
            name: version(name)
            for name in ("jaxprop", "pysolver-view", "CoolProp", "numpy", "scipy")
        },
    )
