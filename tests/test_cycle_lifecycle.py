"""Configured models belong to one problem and are reused only while valid."""

from copy import deepcopy
import pickle

import jaxprop as jp
import numpy as np
import pytest
import yaml

import thermopt as th
from regression_support import (
    EXAMPLE_CASES,
    assert_snapshot,
    load_case,
    nearby_design,
    snapshot,
)


def configuration(case="sco2_recuperated"):
    return deepcopy(load_case(case)["problem_formulation"])


@pytest.fixture(scope="module")
def table_options(test_workspace):
    return dict(
        h_min=110e3,
        h_max=1.3e6,
        p_min=8e6,
        p_max=39e6,
        N_h=32,
        N_p=32,
        N_p_sat=4,
        coarse_step=1,
        table_dir=str(test_workspace / "lifecycle_tables"),
    )


@pytest.fixture(params=["coolprop", "perfect_gas", "bicubic"])
def model_configuration(request, table_options):
    config = configuration()
    settings = config["fixed_parameters"]["working_fluid"]
    settings["method"] = request.param
    if request.param == "bicubic":
        settings["table_options"] = deepcopy(table_options)
    elif request.param == "perfect_gas":
        settings["model_options"] = dict(T_ref=600.0, p_ref=1e5)
    return config


@pytest.mark.parametrize("case", sorted(EXAMPLE_CASES))
def test_every_topology_reuses_its_models_without_setup_during_fitness(
    case, tmp_path, monkeypatch
):
    problem = th.ThermodynamicCycleProblem(configuration(case), out_dir=str(tmp_path))
    expected = snapshot(problem, problem.x0)
    identities = {role: id(fluid) for role, fluid in problem.fluids.items()}

    def unexpected(*args, **kwargs):
        raise AssertionError(
            "Fluid construction/reference preparation belongs to setup"
        )

    monkeypatch.setattr(th.utilities, "fluid_from_configuration", unexpected)
    monkeypatch.setattr(problem, "_calculate_special_points", unexpected)
    monkeypatch.setattr("thermopt.cycles.setup.source_signature", unexpected)
    snapshot(problem, nearby_design(problem))
    assert_snapshot(snapshot(problem, problem.x0), expected)
    assert identities == {role: id(fluid) for role, fluid in problem.fluids.items()}
    assert all(fluid.identifier == role for role, fluid in problem.fluids.items())


def test_unchanged_reload_reuses_native_models_and_reference_points(
    model_configuration, tmp_path, monkeypatch
):
    problem = th.ThermodynamicCycleProblem(model_configuration, out_dir=str(tmp_path))
    models = dict(problem.fluids)
    reference = problem.reference_fluid
    points = dict(problem.params["working_fluid"])

    def unexpected(*args, **kwargs):
        raise AssertionError(
            "An unchanged reload must not reconstruct or reload models"
        )

    monkeypatch.setattr(th.utilities, "fluid_from_configuration", unexpected)
    problem.update_problem(deepcopy(model_configuration))
    assert all(problem.fluids[role] is fluid for role, fluid in models.items())
    assert problem.reference_fluid is reference
    assert all(
        problem.params["working_fluid"][key] is point for key, point in points.items()
    )


def test_nonfluid_changes_reuse_models_and_update_saved_configuration(tmp_path):
    config = configuration()
    problem = th.ThermodynamicCycleProblem(config, out_dir=str(tmp_path))
    previous = dict(problem.fluids)
    changed = deepcopy(config)
    changed["fixed_parameters"]["compressor"]["efficiency"] = 0.8
    changed["fixed_parameters"]["net_power"] = 2e6
    problem.update_problem(changed)
    problem.fitness(problem.x0)
    assert all(problem.fluids[key] is fluid for key, fluid in previous.items())
    assert problem.cycle_data["components"]["compressor"]["efficiency"] == 0.8
    assert problem.cycle_data["energy_analysis"]["net_cycle_power"] == pytest.approx(
        2e6
    )
    target = tmp_path / "reloaded.yaml"
    problem.save_current_configuration(target)
    saved = th.read_configuration_file(target)
    assert saved["fixed_parameters"]["compressor"]["efficiency"] == 0.8
    assert saved["fixed_parameters"]["net_power"] == 2e6
    assert config["fixed_parameters"]["compressor"]["efficiency"] != 0.8


@pytest.mark.parametrize(
    "change", ["fluid_name", "method", "table_domain", "calibration"]
)
def test_real_model_changes_rebuild_only_the_affected_role(
    change, table_options, tmp_path
):
    config = configuration()
    settings = config["fixed_parameters"]["working_fluid"]
    if change == "table_domain":
        settings.update(method="bicubic", table_options=deepcopy(table_options))
    if change == "calibration":
        settings.update(
            method="perfect_gas", model_options=dict(T_ref=600.0, p_ref=1e5)
        )
    problem = th.ThermodynamicCycleProblem(config, out_dir=str(tmp_path))
    previous = dict(problem.fluids)
    changed = deepcopy(config)
    settings = changed["fixed_parameters"]["working_fluid"]
    if change == "fluid_name":
        settings["name"] = "Nitrogen"
    elif change == "method":
        settings.update(
            method="perfect_gas", model_options=dict(T_ref=600.0, p_ref=1e5)
        )
    elif change == "table_domain":
        settings["table_options"]["h_max"] = 1.4e6
    else:
        settings["model_options"]["T_ref"] = 650.0
    problem.update_problem(changed)
    assert problem.fluid is not previous["working_fluid"]
    assert problem.fluid.configuration != previous["working_fluid"].configuration
    for role in ("heating_fluid", "cooling_fluid"):
        assert problem.fluids[role] is previous[role]


def test_changed_reference_temperatures_recompute_points_without_rebuilding_models(
    tmp_path,
):
    config = configuration()
    problem = th.ThermodynamicCycleProblem(config, out_dir=str(tmp_path))
    models = dict(problem.fluids)
    previous = problem.params["working_fluid"]
    changed = deepcopy(config)
    changed["fixed_parameters"]["special_points"]["maximum_temperature"] += 20.0
    problem.update_problem(changed)
    assert all(problem.fluids[role] is fluid for role, fluid in models.items())
    point = problem.params["working_fluid"]["gas_at_maximum_temperature"]
    assert point is not previous["gas_at_maximum_temperature"]
    assert point.temperature == pytest.approx(
        previous["gas_at_maximum_temperature"].temperature + 20.0
    )


@pytest.mark.parametrize(
    "failure", ["invalid_model", "unsupported_property", "invalid_expression"]
)
def test_rejected_reload_keeps_the_last_valid_problem_and_outputs(failure, tmp_path):
    config = configuration()
    problem = th.ThermodynamicCycleProblem(config, out_dir=str(tmp_path))
    original = snapshot(problem, problem.x0)
    models, initial = dict(problem.fluids), problem.x0.copy()
    changed = deepcopy(config)
    if failure == "invalid_model":
        changed["fixed_parameters"]["working_fluid"]["backend"] = "not-a-backend"
        error = jp.InvalidConfigurationError
    elif failure == "unsupported_property":
        changed["fixed_parameters"]["working_fluid"]["method"] = "perfect_gas"
        changed["constraints"] = list(changed["constraints"]) + [
            dict(
                variable="$components.expander.state_out.quality_mass",
                type=">",
                value=0.0,
            )
        ]
        error = jp.UnsupportedPropertyError
    else:
        changed["design_variables"]["compressor_inlet_enthalpy"][
            "value"
        ] = "$missing.reference"
        error = TypeError
    with pytest.raises(error):
        problem.update_problem(changed)
    np.testing.assert_array_equal(problem.x0, initial)
    assert all(problem.fluids[key] is fluid for key, fluid in models.items())
    assert_snapshot(snapshot(problem, initial), original)


def test_caller_mutation_does_not_change_a_configured_problem(tmp_path):
    config = configuration()
    problem = th.ThermodynamicCycleProblem(config, out_dir=str(tmp_path))
    expected = snapshot(problem, problem.x0)
    config["fixed_parameters"]["compressor"]["efficiency"] = 0.5
    assert_snapshot(snapshot(problem, problem.x0), expected)


def test_independent_problems_own_independent_reference_models(tmp_path):
    first = th.ThermodynamicCycleProblem(configuration(), out_dir=str(tmp_path / "one"))
    second = th.ThermodynamicCycleProblem(
        configuration(), out_dir=str(tmp_path / "two")
    )
    for role in first.fluids:
        assert first.fluids[role] is not second.fluids[role]
    assert first.reference_fluid is not second.reference_fluid
    expected = snapshot(first, first.x0)
    snapshot(second, nearby_design(second))
    assert_snapshot(snapshot(first, first.x0), expected)


@pytest.mark.parametrize(
    "operation",
    [deepcopy, lambda value: pickle.loads(pickle.dumps(value))],
    ids=["deepcopy", "pickle"],
)
def test_copied_problems_restore_independent_models_and_identical_results(
    operation, tmp_path
):
    original = th.ThermodynamicCycleProblem(configuration(), out_dir=str(tmp_path))
    expected = snapshot(original, original.x0)
    restored = operation(original)
    assert restored.reference_fluid is not original.reference_fluid
    assert all(
        restored.fluids[role] is not model for role, model in original.fluids.items()
    )
    snapshot(restored, nearby_design(restored))
    assert_snapshot(snapshot(restored, restored.x0), expected)
    assert_snapshot(snapshot(original, original.x0), expected)


def test_cycle_setup_does_not_use_the_shared_study_fluid_cache(tmp_path, monkeypatch):
    def unexpected(*args, **kwargs):
        raise AssertionError("Configured cycle references must be problem-owned")

    monkeypatch.setattr(th.utilities, "get_property_fluid", unexpected)
    problem = th.ThermodynamicCycleProblem(configuration(), out_dir=str(tmp_path))
    assert np.all(np.isfinite(problem.fitness(problem.x0)))


@pytest.mark.parametrize("method", ["plot_cycle", "save_data_to_excel"])
def test_successful_reload_cannot_report_previous_configuration_results(
    method, tmp_path
):
    config = configuration()
    problem = th.ThermodynamicCycleProblem(config, out_dir=str(tmp_path))
    problem.fitness(problem.x0)
    changed = deepcopy(config)
    changed["fixed_parameters"]["net_power"] = 2e6
    problem.update_problem(changed)
    with pytest.raises(RuntimeError, match="Evaluate fitness.*updating"):
        getattr(problem, method)()
    assert np.all(np.isfinite(problem.fitness(problem.x0)))
    assert problem.cycle_data["energy_analysis"]["net_cycle_power"] == pytest.approx(
        2e6
    )


def test_pseudo_pure_reference_interpretation_is_retained(tmp_path):
    config = configuration()
    config["fixed_parameters"]["working_fluid"]["name"] = "Air"
    problem = th.ThermodynamicCycleProblem(config, out_dir=str(tmp_path))
    points = problem.reference_fluid.references
    assert points.triple_liquid.classification == "pseudo_pure"
    assert points.triple_liquid.interpretation == "lower_saturation_endpoint"
    assert (
        problem.params["working_fluid"]["triple_point_liquid"].reference_id
        == points.reference_id
    )


def test_backend_change_rebuilds_only_affected_models_when_refprop_is_available(
    tmp_path,
):
    try:
        jp.Fluid("CO2", backend="REFPROP")
    except jp.InvalidConfigurationError:
        pytest.skip("REFPROP is not available in this installation")
    config = configuration()
    problem = th.ThermodynamicCycleProblem(config, out_dir=str(tmp_path))
    models = dict(problem.fluids)
    changed = deepcopy(config)
    changed["fixed_parameters"]["working_fluid"]["backend"] = "REFPROP"
    problem.update_problem(changed)
    assert problem.fluid.backend == problem.reference_fluid.backend == "REFPROP"
    assert problem.fluid is not models["working_fluid"]
    assert all(
        problem.fluids[role] is models[role]
        for role in ("heating_fluid", "cooling_fluid")
    )
    assert (
        problem.params["working_fluid"]["critical_point"].reference_id
        == problem.reference_fluid.reference_id
    )


def test_source_reference_change_invalidates_models_and_keeps_old_models_pinned(
    tmp_path,
):
    from source_reference_control import set_energy_reference

    set_energy_reference("CO2", "DEF")
    try:
        config = configuration()
        problem = th.ThermodynamicCycleProblem(config, out_dir=str(tmp_path))
        models = dict(problem.fluids)
        original_reference = problem.reference_fluid
        original = problem.fluid.get_state(jp.PT_INPUTS, 9e6, 320.0)
        set_energy_reference("CO2", "ASHRAE")
        problem.update_problem(config)
        current = problem.fluid.get_state(jp.PT_INPUTS, 9e6, 320.0)
        assert problem.fluid is not models["working_fluid"]
        assert problem.reference_fluid is not original_reference
        assert current.reference_id != original.reference_id
        assert abs(float(current.enthalpy - original.enthalpy)) > 1000.0
        assert problem.reference_fluid.reference_id == current.reference_id
        assert (
            problem.params["working_fluid"]["critical_point"].reference_id
            == current.reference_id
        )
        np.testing.assert_allclose(
            models["working_fluid"].get_state(jp.PT_INPUTS, 9e6, 320.0).enthalpy,
            original.enthalpy,
            atol=1e-8,
        )
        assert all(
            problem.fluids[role] is models[role]
            for role in ("heating_fluid", "cooling_fluid")
        )
    finally:
        set_energy_reference("CO2", "DEF")
