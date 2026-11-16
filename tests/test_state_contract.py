"""ThermOpt consumers preserve native state shape, identity and availability."""

from copy import deepcopy
from types import MappingProxyType
from types import SimpleNamespace
from dataclasses import replace
from pathlib import Path
import runpy

import jaxprop as jp
import numpy as np
import pytest

import thermopt as th
from thermopt.components import (
    compression_process,
    expansion_process,
    heat_transfer_process,
)
from thermopt.config import object_to_dict
from thermopt.utilities import (
    fluid_from_configuration,
    render_and_evaluate,
    render_nested_value,
)
from thermopt.cycles.state_requirements import (
    cycle_state_routes,
    validate_cycle_state_requirements,
)
from thermopt.optimize_cycle import CYCLE_TOPOLOGIES
from regression_support import EXAMPLE_CASES, load_case, make_problem


@pytest.fixture(scope="module", params=["coolprop", "perfect_gas", "bicubic"])
def native_model(request, test_workspace):
    settings = dict(name="Air", method=request.param, identifier="working_fluid")
    if request.param == "bicubic":
        settings["table_options"] = dict(
            h_min=200e3,
            h_max=850e3,
            p_min=40e3,
            p_max=350e3,
            N_h=32,
            N_p=32,
            N_p_sat=4,
            coarse_step=1,
            table_dir=str(test_workspace / "state_contract_tables"),
        )
    return fluid_from_configuration(settings)


@pytest.mark.parametrize("steps", [1, 4])
def test_heat_transfer_keeps_scalar_endpoints_and_batched_native_path(
    native_model, steps
):
    inlet = native_model.get_state(jp.PT_INPUTS, 1e5, 300.0)
    result = heat_transfer_process(
        native_model,
        float(inlet.enthalpy),
        1e5,
        float(inlet.enthalpy) + (20e3 if steps > 1 else 0),
        1e5,
        num_steps=steps,
    )
    path = result["states"]
    assert isinstance(path, jp.FluidState) and path.shape == (steps,)
    for name, index in (("state_in", 0), ("state_out", -1)):
        state = result[name]
        assert state.shape == ()
        assert state.identifier == "working_fluid"
        assert state.provider == path.provider == native_model.provider
        assert state.reference_id == path.reference_id == native_model.reference_id
        np.testing.assert_array_equal(state.enthalpy, path.enthalpy[index])
        assert int(state.status) == 0
    before = np.array(path.enthalpy)
    with pytest.raises((AttributeError, TypeError)):
        path.enthalpy = before
    np.testing.assert_array_equal(path.enthalpy, before)


@pytest.mark.parametrize(
    "process,p_in,temperature,p_out",
    [
        (compression_process, 1e5, 300.0, 2e5),
        (expansion_process, 2e5, 600.0, 1e5),
    ],
)
@pytest.mark.parametrize(
    "kind,steps", [("isentropic", 2), ("polytropic", 4), ("non-dimensional", 2)]
)
def test_every_machine_path_is_a_native_batched_state(
    native_model, process, p_in, temperature, p_out, kind, steps
):
    inlet = native_model.get_state(jp.PT_INPUTS, p_in, temperature)
    result = process(
        native_model,
        float(inlet.enthalpy),
        p_in,
        p_out,
        efficiency=0.85,
        efficiency_type=kind,
        mass_flow=2.0,
        num_steps=steps,
    )
    path = result["states"]
    assert isinstance(path, jp.FluidState)
    assert path.shape == (steps,)
    for name, index in (("state_in", 0), ("state_out", -1)):
        state = result[name]
        assert state.shape == ()
        assert path.provider == state.provider == native_model.provider
        assert path.identifier == state.identifier == native_model.identifier
        assert path.reference_id == state.reference_id == native_model.reference_id
        np.testing.assert_array_equal(path.enthalpy[index], state.enthalpy)
    np.testing.assert_array_equal(path.status, 0)


def test_readonly_expression_context_is_supported():
    assert (
        render_nested_value(
            "cycle.pressure",
            MappingProxyType(dict(cycle=MappingProxyType(dict(pressure=1e5)))),
        )
        == 1e5
    )


def test_expression_reads_each_guarded_attribute_once():
    class Counted:
        reads = 0

        @property
        def temperature(self):
            self.reads += 1
            return 320.0

    state = Counted()
    assert render_and_evaluate("$state.temperature + 1", {"state": state}) == 321.0
    assert state.reads == 1


@pytest.mark.parametrize("enabled", [False, True])
def test_constraint_expression_preserves_typed_optional_property_failures(enabled):
    model = jp.Fluid("CO2")
    state = model.get_state(jp.PT_INPUTS, 1e5, 320.0, supersaturation=enabled)
    field = "surface_tension" if enabled else "superheating"
    error_type = jp.PropertyUnavailableError if enabled else jp.UnsupportedPropertyError
    with pytest.raises(error_type) as caught:
        render_and_evaluate(f"$state.{field}", {"state": state})
    assert caught.value.request == field
    assert caught.value.provider == model.provider
    assert any("$state." in note for note in caught.value.__notes__)


def test_state_serialization_uses_lossless_public_data(native_model):
    state = native_model.get_state(
        jp.PT_INPUTS, np.array([1e5, -1.0]), 300.0, errors="mask"
    )
    serialized = object_to_dict(dict(state=state))["state"]
    assert "_data" not in serialized and "_clipped" not in serialized
    restored = jp.FluidState.from_data(serialized)
    assert restored.provider == state.provider
    assert restored.identifier == state.identifier
    assert restored.reference_id == state.reference_id
    assert restored.available_properties == state.available_properties
    np.testing.assert_array_equal(restored.status, state.status)
    np.testing.assert_array_equal(restored.temperature, state.temperature)


@pytest.mark.parametrize("field", ["quality_mass", "Q", "superheating"])
@pytest.mark.parametrize("location", ["variable", "target", "objective"])
def test_unsupported_state_requirements_fail_during_problem_setup(
    field, location, tmp_path
):
    config = deepcopy(load_case("sco2_recuperated")["problem_formulation"])
    config["constraints"] = list(config["constraints"])
    config["fixed_parameters"]["working_fluid"]["method"] = "perfect_gas"
    expression = f"$components.expander.state_out.{field}"
    if location == "objective":
        config["objective_function"]["variable"] = expression
    else:
        config["constraints"].append(
            dict(
                variable=(
                    expression
                    if location == "variable"
                    else "$energy_analysis.cycle_efficiency"
                ),
                type=">",
                value=expression if location == "target" else 0.0,
            )
        )
    with pytest.raises(jp.UnsupportedPropertyError) as caught:
        th.ThermodynamicCycleProblem(config, out_dir=str(tmp_path))
    assert caught.value.request == ("quality_mass" if field == "Q" else field)
    assert "expander" in str(caught.value)


def test_supported_but_disabled_heat_exchanger_diagnostic_is_explicit(tmp_path):
    config = deepcopy(load_case("sco2_recuperated")["problem_formulation"])
    config["constraints"] = list(config["constraints"])
    config["constraints"].append(
        dict(
            variable="$components.heater.cold_side.state_out.superheating",
            type=">",
            value=0.0,
        )
    )
    with pytest.raises(jp.UnsupportedPropertyError, match="not enabled"):
        th.ThermodynamicCycleProblem(config, out_dir=str(tmp_path))


def test_unknown_state_path_is_rejected_at_setup(tmp_path):
    config = deepcopy(load_case("sco2_recuperated")["problem_formulation"])
    config["constraints"] = list(config["constraints"])
    config["constraints"].append(
        dict(variable="$components.missing.state_out.temperature", type=">", value=0.0)
    )
    with pytest.raises(jp.InvalidConfigurationError, match="missing"):
        th.ThermodynamicCycleProblem(config, out_dir=str(tmp_path))


@pytest.mark.parametrize("case", sorted(EXAMPLE_CASES))
def test_declared_routes_match_actual_component_models(case, tmp_path):
    problem = make_problem(case, tmp_path)
    problem.fitness(problem.x0)
    for path, (role, diagnostic_path) in cycle_state_routes(
        problem.cycle_topology
    ).items():
        stream = problem.cycle_data["components"]
        for part in path.split("."):
            stream = stream[part]
        model = problem.fluids[role]
        assert isinstance(stream["states"], jp.FluidState)
        for state in (stream["states"], stream["state_in"], stream["state_out"]):
            assert state.identifier == role
            assert state.provider == model.provider
            assert state.reference_id == model.reference_id
            expected = (
                diagnostic_path and "supersaturation" in model.capabilities.options
            )
            assert state.supports_property("superheating") == expected


@pytest.mark.parametrize("topology", sorted(CYCLE_TOPOLOGIES))
def test_all_public_topology_aliases_have_state_routes(topology):
    assert cycle_state_routes(topology)


def test_capability_checks_are_setup_only(tmp_path, monkeypatch):
    import thermopt.optimize_cycle as module

    problem = make_problem("sco2_recuperated", tmp_path)

    def unexpected(*args):
        raise AssertionError(
            "repeated fitness must not parse or validate model requirements"
        )

    monkeypatch.setattr(module, "validate_cycle_state_requirements", unexpected)
    problem.fitness(problem.x0)
    problem.fitness(problem.x0)


@pytest.mark.parametrize(
    "missing", ["viscosity", "speed_of_sound", "HmassSmass_INPUTS"]
)
def test_radial_machine_requirements_fail_before_numerical_calls(missing):
    model = jp.Fluid("CO2")
    config = load_case("ptes_turbo")["problem_formulation"]
    roles = ("working_fluid", "hot_storage_fluid", "cold_storage_fluid")
    caps = model.capabilities
    if missing == "HmassSmass_INPUTS":
        caps = replace(caps, input_pairs=caps.input_pairs - {jp.HmassSmass_INPUTS})
    else:
        caps = replace(caps, properties=caps.properties - {missing})
    restricted = SimpleNamespace(
        name=model.name, provider=model.provider, model=model.model, capabilities=caps
    )
    fluids = {role: restricted if role == "working_fluid" else model for role in roles}
    kind = (
        jp.UnsupportedInputPairError
        if missing == "HmassSmass_INPUTS"
        else jp.UnsupportedPropertyError
    )
    with pytest.raises(kind) as caught:
        validate_cycle_state_requirements(
            "PTES_recuperated_turbo",
            config["fixed_parameters"],
            fluids,
            [],
            config["objective_function"],
        )
    assert caught.value.request == missing


@pytest.mark.parametrize(
    "x_alias,x_canonical",
    [("h", "enthalpy"), ("s", "entropy"), ("heat_flow", "heat_flow")],
)
def test_auxiliary_plot_streams_support_canonical_axis_names(
    x_alias, x_canonical, tmp_path
):
    problem = make_problem("sco2_recuperated", tmp_path)
    problem.fitness(problem.x0)
    expected = problem._get_process_data("heater_hot_side", x_alias, "T")
    actual = problem._get_process_data("heater_hot_side", x_canonical, "temperature")
    assert actual[0] is not None and actual[1] is not None
    np.testing.assert_array_equal(actual[0], expected[0])
    np.testing.assert_array_equal(actual[1], expected[1])


def test_source_quality_indicator_is_not_treated_as_a_physical_fraction():
    model = jp.Fluid("Water")
    inlet = model.get_state(jp.PT_INPUTS, 2e5, 500.0)
    result = expansion_process(model, float(inlet.enthalpy), 2e5, 1e5, efficiency=0.85)
    state = result["state_out"]
    exact = model.get_state(
        jp.HmassP_INPUTS,
        state.enthalpy,
        state.pressure,
        generalize_quality=True,
        supersaturation=True,
    )
    trimmed = model.get_state(jp.HmassP_INPUTS, state.enthalpy, state.pressure)
    assert not bool(state.is_two_phase)
    assert float(state.quality_mass) > 1.0
    assert float(trimmed.quality_mass) in (0.0, 1.0)
    np.testing.assert_allclose(
        state.quality_mass, exact.quality_mass, rtol=1e-12, atol=1e-12
    )
    assert state.supports_property("quality_volume") and not bool(
        state.property_valid("quality_volume")
    )
    with pytest.raises(jp.PropertyUnavailableError):
        _ = state.quality_volume
    assert np.isfinite(state.temperature) and int(state.status) == 0


def test_bicubic_diagnostic_constraint_fails_before_reference_preparation(
    tmp_path, monkeypatch
):
    config = deepcopy(load_case("sco2_recuperated")["problem_formulation"])
    config["fixed_parameters"]["working_fluid"].update(
        method="bicubic",
        table_options=dict(
            h_min=600e3,
            h_max=850e3,
            p_min=5e6,
            p_max=12e6,
            N_h=4,
            N_p=4,
            N_p_sat=4,
            coarse_step=1,
            table_dir=str(tmp_path / "tables"),
        ),
    )
    config["constraints"] = [
        dict(
            variable="$components.expander.state_out.superheating", type=">", value=0.0
        )
    ]

    def unexpected(*args):
        raise AssertionError(
            "unsupported constraint must fail before preparing source references"
        )

    monkeypatch.setattr("thermopt.optimize_cycle.prepare_cycle_references", unexpected)
    with pytest.raises(jp.UnsupportedPropertyError) as caught:
        th.ThermodynamicCycleProblem(config, out_dir=str(tmp_path))
    assert caught.value.provider == "bicubic"
    assert caught.value.request == "superheating"


def test_blended_heat_transfer_example_uses_public_state_stack(tmp_path, monkeypatch):
    script = (
        Path(__file__).resolve().parents[1]
        / "case_studies/HX_development/blended_calculation.py"
    )
    monkeypatch.chdir(tmp_path)
    result = runpy.run_path(str(script), run_name="__main__")
    assert isinstance(result["states"], jp.FluidState)
    assert result["states"].shape == (300,)
    assert np.all(np.isfinite(result["states"].temperature))
    assert np.all(
        (np.asarray(result["htc_blend_list"]) >= 500.0)
        & (np.asarray(result["htc_blend_list"]) <= 5000.0)
    )
    assert len(list(tmp_path.glob("*.png"))) == 2
