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
# The dated baseline retains its original paths and physical hashes. Remap only
# lookup locations after the examples-to-case_studies directory migration.
BASELINES["cases"] = {
    (
        "case_studies/" + name[len("examples/") :]
        if name.startswith("examples/")
        else name
    ): record
    for name, record in BASELINES["cases"].items()
}
SOLVER_POLICY = json.loads(
    (ROOT / "tests/data/case_solver_defaults_20261116.json").read_text(encoding="utf-8")
)
SOLVER_DEFAULTS = SOLVER_POLICY["cases"]
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
        "case_studies",
        "tests/data/cases",
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
    recompression = (
        filename == "case_studies/sCO2_recompression/case_sCO2_recompression.yaml"
    )
    for role, settings in streams.items():
        expected_method = (
            "bicubic" if recompression and role == "working_fluid" else "coolprop"
        )
        assert settings["method"] == expected_method
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
    compared = deepcopy(data)
    if filename in SOLVER_DEFAULTS:
        # Review only the requested solver policy; compare all other setup and
        # physical values against their original, unchanged historical record.
        options = compared["solver_options"]
        previous = SOLVER_DEFAULTS[filename]["previous_solver_options"]
        # Actual launcher checks exposed order-sensitive finite differences.
        # Review explicit numerical fixes separately from immutable physics.
        for key, expected in (
            SOLVER_POLICY.get("reviewed_numerical_overrides", {})
            .get(filename, {})
            .items()
        ):
            assert options[key] == expected
            options[key] = previous[key]
        assert options["library"] == "pysolver" and options["method"] == "slsqp"
        assert options["update_on"] == "iteration"
        extra = options["extra_options"]
        assert set(extra) == {
            "first_step_fraction",
            "max_step_fraction",
            "step_growth_factor",
            "feasibility_tol",
            "optimality_tol",
        }
        assert (
            extra["first_step_fraction"],
            extra["max_step_fraction"],
            extra["step_growth_factor"],
        ) == (0.2, 1.0, 1.5)
        assert float(extra["feasibility_tol"]) == float(previous["tolerance"])
        assert float(extra["optimality_tol"]) == 1e-4
        unchanged_keys = set(options) - {"library", "update_on", "extra_options"}
        assert {key: options[key] for key in unchanged_keys} == {
            key: value
            for key, value in previous.items()
            if key not in {"library", "update_on", "extra_options"}
        }
        compared["solver_options"] = deepcopy(previous)
    if recompression:
        # Explicit user-requested evaluator/optimizer migration, 16 November 2026.
        # Keep the original physical hash: only these reviewed setup changes
        # are normalized; all design bounds, fluids and cycle physics stay checked.
        fluid = compared["problem_formulation"]["fixed_parameters"]["working_fluid"]
        assert fluid.pop("table_options") == dict(
            h_min="150e3",
            h_max="1.3e6",
            p_min="8e6",
            p_max="42e6",
            N_h=200,
            N_p=200,
            N_p_sat=8,
            gradient_method="central",
            coarse_step=1,
            table_dir="./results/property_tables",
        )
        fluid["method"] = "coolprop"
        options = compared["solver_options"]
        assert options["library"] == "pysolver"
        assert options.pop("extra_options") == dict(
            first_step_fraction=0.2,
            max_step_fraction=1.0,
            step_growth_factor=1.5,
            feasibility_tol="1e-6",
            optimality_tol="1e-4",
        )
        options["library"] = "scipy"
    digest = hashlib.sha256(
        json.dumps(physical_values(compared), sort_keys=True, allow_nan=False).encode()
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
