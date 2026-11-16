"""Prepare reusable native models and reference data outside numerical calls.

Reuse is local to a configured problem. Source signatures use JaxProp's own
table/serialization provenance, including active energy references and native
definitions. Probing and comparing that metadata is setup work only; fitness
and component state calculations never call these functions.
"""

from copy import deepcopy

import jaxprop as jp
from jaxprop.bicubic.cache_metadata import metadata_fingerprint, source_signature

from .. import utilities as utils
from ..utilities.fluid_utils import _scalar_settings

STREAM_ROLES = (
    "working_fluid",
    "heating_fluid",
    "cooling_fluid",
    "hot_storage_fluid",
    "cold_storage_fluid",
)


def prepare_cycle_models(
    parameters, previous_models=None, previous_settings=None, previous_sources=None
):
    """Return problem-owned evaluators, construction settings and source keys.

    A real configuration or source change reconstructs only affected roles.
    Backends without deterministic provenance and composite models without one
    source backend are conservatively reconstructed on reload. No native model
    is shared globally or reconstructed during cycle evaluation.
    """
    previous_models = {} if previous_models is None else previous_models
    previous_settings = {} if previous_settings is None else previous_settings
    previous_sources = {} if previous_sources is None else previous_sources
    probed = {}

    def identity(model):
        backend = getattr(model, "backend", None)
        if backend is None:
            return None
        pair = (model.name, backend)
        if pair not in probed:
            signature = source_signature(*pair)
            probed[pair] = (
                (
                    metadata_fingerprint(signature)
                    if signature["deterministic_cache_identity"]
                    else None
                ),
                signature["source_reference_id"],
            )
        return probed[pair]

    models, settings, sources = {}, {}, {}
    for role in STREAM_ROLES:
        if role not in parameters:
            continue
        current = _scalar_settings(dict(parameters[role]))
        current["identifier"] = role
        previous = previous_models.get(role)
        reusable = False
        if (
            previous is not None
            and getattr(previous, "backend", None) is not None
            and current == previous_settings.get(role)
        ):
            key = identity(previous)
            reusable = (
                key is not None
                and key[0] is not None
                and key == previous_sources.get(role)
            )
        model = previous if reusable else utils.fluid_from_configuration(current)
        key = identity(model)
        if key is not None and model.references.reference_id != key[1]:
            raise jp.InvalidConfigurationError(
                fluid_name=model.name,
                provider=model.provider,
                model=model.model,
                request="source_reference",
                reason="The source convention changed while preparing cycle models; retry setup under a stable source configuration.",
            )
        models[role], settings[role], sources[role] = model, current, key
    return models, settings, sources


def prepare_cycle_references(
    parameters,
    working_model,
    source_key,
    previous_reference=None,
    previous_key=None,
    previous_points=None,
):
    """Return an independently owned source evaluator and immutable point data.

    Critical/triple records remain the source's own records, including honest
    pseudo-pure lower-saturation endpoint interpretation. Unsupported references
    raise JaxProp's informative error rather than inventing replacement values.
    Source initialization can lie outside a numerical table's declared domain.
    """
    # Resolve availability before trying to construct a single source for a
    # component mixture that has no unique source critical/triple point.
    for name in ("critical_point", "triple_point_liquid", "triple_point_vapor"):
        if not working_model.references[name].available:
            getattr(working_model, name)  # Raise the native contextual error.
    special = parameters["special_points"]
    key = (
        source_key,
        float(special["ambient_temperature"]),
        float(special["maximum_temperature"]),
    )
    source_reusable = (
        source_key is not None
        and source_key[0] is not None
        and previous_reference is not None
        and previous_key is not None
        and source_key == previous_key[0]
    )
    reference = (
        previous_reference
        if source_reusable
        else utils.fluid_from_configuration(
            dict(
                name=working_model.name,
                method="coolprop",
                backend=working_model.backend,
                identifier="working_fluid_reference",
            )
        )
    )
    if reference.references.reference_id != working_model.references.reference_id:
        raise jp.InvalidConfigurationError(
            fluid_name=working_model.name,
            provider=working_model.provider,
            model=working_model.model,
            request="source_reference",
            reason="Cycle reference preparation and the numerical model use different source conventions.",
        )
    if source_reusable and key == previous_key and previous_points is not None:
        points = previous_points
    else:
        points = calculate_reference_points(parameters, reference)
    params = deepcopy(parameters)
    params["working_fluid"] = dict(points)
    return reference, key, points, params


def calculate_reference_points(parameters, reference):
    """Evaluate source-EOS cycle reference points with explicit inputs.

    Below the critical temperature use saturated liquid at the sink temperature;
    above it retain the established constant-critical-density construction.
    Evaluate dilute gas at 1.01 times the liquid lower-endpoint pressure. These
    equations intentionally preserve the existing scientific initialization.
    """
    try:
        critical = reference.critical_point
        liquid, vapor = reference.triple_point_liquid, reference.triple_point_vapor
        sink = parameters["special_points"]["ambient_temperature"]
        source = parameters["special_points"]["maximum_temperature"]
        if sink < critical.temperature:
            state_sat = reference.get_state(jp.QT_INPUTS, 0.0, sink)
        else:
            state_sat = reference.get_state(jp.DmassT_INPUTS, critical.density, sink)
        state_dilute = reference.get_state(jp.PT_INPUTS, 1.01 * liquid.pressure, source)
    except jp.JaxPropError as error:
        error.add_note(
            "While preparing working_fluid source-EOS initialization points for the cycle."
        )
        raise
    return dict(
        critical_point=critical,
        triple_point_liquid=liquid,
        triple_point_vapor=vapor,
        liquid_at_ambient_temperature=state_sat,
        gas_at_maximum_temperature=state_dilute,
    )
