"""Contracts for configuration, residual signs, and report serialization.

These checks use explicit expected values so later property/backend refactors
cannot silently change the optimizer's numerical conventions.
"""

from copy import deepcopy

import jax.numpy as jnp
import jaxprop as jp
import numpy as np
import pytest
import yaml

from thermopt.config import convert_configuration_options, read_configuration_file
from thermopt.utilities.file_utils import convert_numpy_to_python
from thermopt.utilities.optimization_utils import (
    evaluate_constraints,
    evaluate_objective_function,
    render_and_evaluate,
    render_nested_value,
)


def test_yaml_converts_units_expressions_without_changing_model_selectors(tmp_path):
    path = tmp_path / "configuration.yaml"
    path.write_text(
        "temperature: 25 + 273.15\n"
        "pressure: 2e5\n"
        "backend: HEOS\n"
        "enabled: true\n"
        "steps: [5, 10, 20]\n"
        "initial_enthalpy: 1.2*$working_fluid.critical_point.h\n"
        "nested:\n  pressure_drop: 1 / 100\n",
        encoding="utf-8",
    )
    result = read_configuration_file(path)
    assert result["temperature"] == pytest.approx(298.15)
    assert result["pressure"] == pytest.approx(200_000)
    assert result["nested"]["pressure_drop"] == pytest.approx(0.01)
    np.testing.assert_array_equal(result["steps"], [5, 10, 20])
    assert result["backend"] == "HEOS"
    assert result["enabled"] is True
    assert result["initial_enthalpy"] == "1.2*$working_fluid.critical_point.h"


def test_configuration_conversion_does_not_mutate_the_input():
    original = {"temperature": "300 + 20", "elements": [5, 10], "flags": []}
    before = deepcopy(original)
    result = convert_configuration_options(original)
    assert original == before
    assert result["temperature"] == 320
    assert isinstance(result["elements"], np.ndarray)
    assert result["flags"] == []


def test_configuration_rejects_mixed_numeric_and_text_array():
    with pytest.raises(ValueError, match="elements.*different types"):
        convert_configuration_options({"elements": [5, "HEOS"]})


@pytest.mark.parametrize("content", [None, "working_fluid: [unterminated\n"])
def test_unreadable_or_malformed_yaml_has_file_context(tmp_path, content):
    path = tmp_path / "bad_configuration.yaml"
    if content is not None:
        path.write_text(content, encoding="utf-8")
    with pytest.raises(Exception, match="Error parsing configuration file") as caught:
        read_configuration_file(path)
    assert str(path) in str(caught.value)


def test_expressions_resolve_fluid_state_aliases_and_vector_constraints():
    data = {
        "component": {
            "state": jp.FluidState(enthalpy=np.float64(400_000), temperature=np.float64(320)),
            "temperature_difference": np.array([10.0, 20.0, 30.0]),
        }
    }
    assert render_and_evaluate("$component.state.h / 1000", data) == 400
    assert render_nested_value("component.state.T", data) == 320
    np.testing.assert_allclose(
        render_and_evaluate("$component.temperature_difference - 5", data),
        [5, 15, 25],
    )


def test_expressions_accept_jax_arrays_without_changing_numeric_results():
    result = render_and_evaluate("$pinch / 10", {"pinch": jnp.array([20.0, 30.0])})
    np.testing.assert_allclose(result, [2, 3])


@pytest.mark.parametrize(
    ("expression", "error", "match"),
    [
        ("$missing.temperature", TypeError, "missing.temperature"),
        ("$label", TypeError, "not numeric"),
        ("$pressure +", SyntaxError, "Syntax error"),
    ],
)
def test_expression_errors_identify_the_invalid_expression(expression, error, match):
    with pytest.raises(error, match=match):
        render_and_evaluate(expression, {"pressure": 100_000.0, "label": "HEOS"})


@pytest.mark.parametrize(("kind", "sign"), [("minimize", 1), ("maximize", -1)])
def test_net_work_objective_has_expected_optimization_sign(kind, sign):
    data = {"work": {"turbine": 12.0, "compressor": 4.0}}
    result = evaluate_objective_function(
        data,
        {"variable": "$work.turbine - $work.compressor", "type": kind},
    )
    assert float(result) == pytest.approx(sign * 8)


@pytest.mark.parametrize(
    ("definition", "match"),
    [
        ({"variable": "$values", "type": "minimize"}, "must be scalar"),
        ({"variable": "$power", "type": "unknown"}, "Unknown objective"),
    ],
)
def test_invalid_objective_is_rejected(definition, match):
    with pytest.raises(ValueError, match=match):
        evaluate_objective_function({"values": np.array([1, 2]), "power": 10}, definition)


@pytest.mark.parametrize(("sense", "sign"), [("<", 1), (">", -1)])
@pytest.mark.parametrize(("normalize", "divisor"), [(False, 1), (True, 20), (100.0, 100)])
def test_exchanger_constraints_keep_node_order_sign_and_normalization(sense, sign, normalize, divisor):
    # The solver convention is c_ineq <= 0, regardless of the YAML sense.
    equalities, inequalities, report = evaluate_constraints(
        {"pinch": np.array([10.0, 20.0, 40.0]), "limit": 20.0},
        [{"variable": "$pinch", "type": sense, "value": "$limit", "normalize": normalize}],
    )
    assert equalities == []
    np.testing.assert_allclose(inequalities, sign * np.array([-10, 0, 20]) / divisor)
    assert [row["name"] for row in report] == ["$pinch[0]", "$pinch[1]", "$pinch[2]"]
    np.testing.assert_allclose([row["value"] for row in report], [10, 20, 40])
    expected_satisfied = [True, True, False] if sense == "<" else [False, True, True]
    assert [bool(row["satisfied"]) for row in report] == expected_satisfied


def test_energy_balance_equalities_and_inequalities_are_returned_separately():
    data = {"energy_balance": 3.0, "power": 90.0}
    equalities, inequalities, report = evaluate_constraints(
        data,
        [
            {"variable": "$power", "type": ">", "value": 100.0},
            {"variable": "$energy_balance", "type": "=", "value": 0.0, "normalize": 100.0},
        ],
    )
    np.testing.assert_allclose(equalities, [0.03])
    np.testing.assert_allclose(inequalities, [10.0])
    assert [row["name"] for row in report] == ["$power", "$energy_balance"]
    assert not any(bool(row["satisfied"]) for row in report)


@pytest.mark.parametrize(
    ("definition", "match"),
    [
        ({"variable": "$balance", "type": "=", "value": 0.0, "normalize": True}, "normalization factor is 0"),
        ({"variable": "$balance", "type": "=", "value": 5.0, "normalize": 0.0}, "normalization factor is 0"),
        ({"variable": "$balance", "type": "invalid", "value": 0.0}, "Unknown constraint"),
    ],
)
def test_invalid_constraint_definitions_fail_before_solver_use(definition, match):
    with pytest.raises(ValueError, match=match):
        evaluate_constraints({"balance": 3.0}, [definition])


def test_numeric_configuration_can_be_exported_and_read_back_as_yaml(tmp_path):
    data = {
        "temperature": np.array(298.15),
        "pressures": np.array([100_000.0, 200_000.0]),
        "enabled": np.bool_(True),
        "iterations": np.int64(20),
        "backend": np.str_("HEOS"),
        "optional": None,
        "nested": {"coefficient": np.float64(1 / 3)},
    }
    serializable = convert_numpy_to_python(data, precision=8)
    path = tmp_path / "roundtrip.yaml"
    path.write_text(yaml.safe_dump(serializable), encoding="utf-8")
    result = read_configuration_file(path)
    assert result["temperature"] == pytest.approx(298.15)
    np.testing.assert_allclose(result["pressures"], data["pressures"])
    assert result["enabled"] is True
    assert result["iterations"] == 20
    assert result["backend"] == "HEOS"
    assert result["optional"] is None
    assert result["nested"]["coefficient"] == pytest.approx(1 / 3, abs=1e-8)
