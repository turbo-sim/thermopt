"""The new constructor vocabulary and unchanged complete case inputs."""

from copy import deepcopy
import hashlib
import json
from pathlib import Path
from types import MappingProxyType
import warnings

import jaxprop as jp
import numpy as np
import pytest
import yaml

from thermopt.utilities.fluid_utils import fluid_from_configuration

ROOT = Path(__file__).resolve().parents[1]
BASELINES = json.loads(
    (ROOT / "tests/data/fluid_configuration_phase2.json").read_text(encoding="utf-8")
)
ROLES = {
    "working_fluid",
    "heating_fluid",
    "cooling_fluid",
    "hot_storage_fluid",
    "cold_storage_fluid",
}


def test_all_maintained_cycle_inputs_are_in_the_reviewed_inventory():
    found = set()
    archived = {"initializers", "thesis_results", "figure_preparation", "archive"}
    for directory in (
        "assignments",
        "docs/dev",
        "examples",
        "tests/data/cases",
        "case_studies",
    ):
        for path in (ROOT / directory).rglob("*.yaml"):
            parts = path.relative_to(ROOT).parts
            if any(
                part in archived or part.startswith(("output", "results"))
                for part in parts
            ):
                continue
            data = yaml.safe_load(path.read_text(encoding="utf-8-sig"))
            if not isinstance(data, dict):
                continue
            formulation = data.get("problem_formulation", data)
            if isinstance(formulation, dict) and ROLES.intersection(
                formulation.get("fixed_parameters", {})
            ):
                found.add(path.relative_to(ROOT).as_posix())
    assert found == set(BASELINES["cases"])


def physical_values(value):
    """Compare physics, allowing only the new explicit source selector."""
    if isinstance(value, dict):
        result = {k: physical_values(v) for k, v in value.items()}
        if (
            "name" in result
            and "backend" in result
            and result.get("method") == "coolprop"
        ):
            result.pop("method")
        return result
    if isinstance(value, list):
        return [physical_values(v) for v in value]
    return value


@pytest.mark.parametrize("filename", sorted(BASELINES["cases"]))
def test_maintained_case_schema_and_physical_values_are_preserved(filename):
    data = yaml.safe_load((ROOT / filename).read_text(encoding="utf-8"))
    formulation = data.get("problem_formulation", data)
    streams = {k: v for k, v in formulation["fixed_parameters"].items() if k in ROLES}
    assert set(streams) == set(BASELINES["cases"][filename]["streams"])
    for settings in streams.values():
        assert settings["method"] == "coolprop"
        assert settings["backend"] == "HEOS"
        assert "exceptions" not in settings
        assert set(settings) <= {
            "name",
            "method",
            "backend",
            "identifier",
            "table_options",
            "model_options",
        }
    digest = hashlib.sha256(
        json.dumps(physical_values(data), sort_keys=True, allow_nan=False).encode()
    ).hexdigest()
    assert digest == BASELINES["cases"][filename]["physical_sha256"]


@pytest.mark.parametrize(
    "method,expected",
    [
        ("coolprop", jp.FluidCoolProp),
        ("perfect_gas", jp.FluidPerfectGas),
        ("bicubic", jp.FluidBicubic),
    ],
)
def test_configuration_returns_native_evaluator_without_mutating_input(
    method, expected, tmp_path
):
    settings = dict(name="CO2", method=method, backend="HEOS", identifier="configured")
    if method == "perfect_gas":
        settings["model_options"] = dict(T_ref=np.float64(600.0), p_ref=np.float64(1e5))
    if method == "bicubic":
        settings["table_options"] = dict(
            h_min=np.float64(600e3),
            h_max=850e3,
            p_min=5e6,
            p_max=12e6,
            N_h=np.int64(4),
            N_p=4,
            N_p_sat=4,
            coarse_step=1,
            table_dir=str(tmp_path),
        )
    original = deepcopy(settings)
    with warnings.catch_warnings(record=True) as messages:
        warnings.simplefilter("always")
        model = fluid_from_configuration(settings, identifier="working_fluid")
    assert type(model) is expected
    assert model.identifier == "working_fluid"
    assert settings == original
    assert not any("array is being set as static" in str(w.message) for w in messages)
    state = model.get_state(jp.HmassP_INPUTS, 700e3, 8e6)
    assert state.shape == () and int(state.status) == 0
    assert np.isfinite(float(state.temperature))


def test_readonly_nested_configuration_normalizes_static_scalars_once():
    options = MappingProxyType(dict(T_ref=np.float64(600.0), p_ref=np.float64(1e5)))
    settings = MappingProxyType(
        dict(
            name="Air", method="perfect_gas", model_options=options, identifier="inlet"
        )
    )
    model = fluid_from_configuration(settings)
    assert model.identifier == "inlet"
    assert type(model.configuration.model_options["T_ref"]) is float
    assert isinstance(options["T_ref"], np.float64)


@pytest.mark.parametrize("legacy_value", [True, False])
def test_legacy_constructor_policy_is_rejected_not_translated(legacy_value):
    with pytest.raises(jp.InvalidConfigurationError, match="error policy") as caught:
        fluid_from_configuration(
            dict(name="Air", method="coolprop", exceptions=legacy_value)
        )
    assert caught.value.request == "exceptions"


@pytest.mark.parametrize(
    "extra",
    [
        dict(method="HEOS"),
        dict(provider="coolprop"),
        dict(unknown=True),
        dict(errors="mask"),
    ],
)
def test_invalid_constructor_fields_are_not_dropped(extra):
    with pytest.raises(jp.InvalidConfigurationError):
        fluid_from_configuration(dict(name="Air", method="coolprop") | extra)


def test_old_positional_backend_is_not_the_new_method_selector():
    with pytest.raises(jp.InvalidConfigurationError, match="method"):
        jp.Fluid("Air", "HEOS")


def test_error_policy_belongs_to_the_state_call():
    model = fluid_from_configuration(
        dict(name="Air", method="coolprop", backend="HEOS")
    )
    with pytest.raises(jp.InvalidStateError):
        model.get_state(jp.PT_INPUTS, -1.0, 300.0)
    state = model.get_state(jp.PT_INPUTS, np.array([1e5, -1.0]), 300.0, errors="mask")
    assert state.status[0] == 0 and state.status[1] != 0
    assert np.isnan(state.temperature[1])
