"""Manager configuration and user-facing export workflows, without a solve.

All cases are frozen test inputs; exports live in the test workspace. Numerical
readback checks protect the public YAML/Excel/report interface before backend
changes, rather than merely checking that a file was created.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import yaml

import thermopt as th

from regression_support import assert_snapshot, load_case, nearby_design, snapshot


@pytest.fixture
def cycle_manager(tmp_path):
    config = load_case("sco2_recuperated")
    config["solver_options"]["callbacks"] = {}
    config["solver_options"]["print_convergence"] = False
    config["solver_options"]["plot_convergence"] = False
    config["solver_options"]["max_iterations"] = int(config["solver_options"]["max_iterations"])
    config["problem_formulation"]["plot_settings"]["fluid"] = {
        "plot_saturation_line": False,
        "plot_critical_point": False,
    }
    path = tmp_path / "input.yaml"
    path.write_text(yaml.safe_dump(th.convert_numpy_to_python(config, precision=14)), encoding="utf-8")
    return th.ThermodynamicCycleOptimization(str(path), out_dir=str(tmp_path / "results"))


def test_single_configuration_update_rebuilds_problem_and_changes_compressor(cycle_manager):
    manager = cycle_manager
    previous_problem, previous_solver = manager.problem, manager.solver
    previous_outlet = float(previous_problem.cycle_data["components"]["compressor"]["state_out"].h)
    previous_efficiency = float(previous_problem.cycle_data["energy_analysis"]["cycle_efficiency"])

    manager.set_config_value("problem_formulation.fixed_parameters.compressor.efficiency", 0.80)

    assert manager.problem is not previous_problem
    assert manager.solver is not previous_solver
    assert manager.solver.problem is manager.problem
    compressor = manager.problem.cycle_data["components"]["compressor"]
    assert compressor["efficiency"] == pytest.approx(0.80)
    assert float(compressor["state_out"].h) > previous_outlet
    assert float(manager.problem.cycle_data["energy_analysis"]["cycle_efficiency"]) < previous_efficiency
    assert np.all(np.isfinite(manager.problem.fitness(manager.problem.x0)))


def test_deferred_configuration_update_requires_reload_to_change_cycle(cycle_manager):
    manager = cycle_manager
    previous_problem, previous_solver = manager.problem, manager.solver
    original_power = float(previous_problem.cycle_data["energy_analysis"]["net_cycle_power"])
    original_flow = float(previous_problem.cycle_data["energy_analysis"]["mass_flow_working_fluid"])

    manager.set_config_value("problem_formulation.fixed_parameters.net_power", 2 * original_power, reload=False)

    assert manager.problem is previous_problem
    assert manager.solver is previous_solver
    manager.problem.fitness(manager.problem.x0)
    assert manager.problem.cycle_data["energy_analysis"]["net_cycle_power"] == pytest.approx(original_power)
    manager.load_config(manager.config)
    assert manager.problem is not previous_problem
    assert manager.problem.cycle_data["energy_analysis"]["net_cycle_power"] == pytest.approx(2 * original_power)
    assert manager.problem.cycle_data["energy_analysis"]["mass_flow_working_fluid"] == pytest.approx(2 * original_flow)


def test_batch_configuration_updates_are_applied_with_one_reload(cycle_manager, monkeypatch):
    manager = cycle_manager
    reloads = []
    original_load = manager.load_config

    def observed_load(config):
        reloads.append(config)
        return original_load(config)

    monkeypatch.setattr(manager, "load_config", observed_load)
    manager.set_config_values({
        "problem_formulation.fixed_parameters.net_power": 2e6,
        "problem_formulation.fixed_parameters.compressor.efficiency": 0.82,
        "problem_formulation.fixed_parameters.expander.efficiency": 0.88,
    })

    assert len(reloads) == 1
    data = manager.problem.cycle_data
    assert data["energy_analysis"]["net_cycle_power"] == pytest.approx(2e6)
    assert data["components"]["compressor"]["efficiency"] == pytest.approx(0.82)
    assert data["components"]["expander"]["efficiency"] == pytest.approx(0.88)
    assert np.all(np.isfinite(manager.problem.fitness(manager.problem.x0)))


def test_constraint_update_changes_every_heater_node_without_adding_duplicates(cycle_manager):
    manager = cycle_manager
    name = "$components.heater.temperature_difference"
    count_before = len(manager.config["problem_formulation"]["constraints"])
    manager.set_constraint(name, type=">", value=25.0, normalize=False)

    constraints = manager.config["problem_formulation"]["constraints"]
    assert len(constraints) == count_before
    assert sum(item["variable"] == name for item in constraints) == 1
    assert manager.problem.get_nic() == 75
    gap = np.asarray(manager.problem.cycle_data["components"]["heater"]["temperature_difference"])
    np.testing.assert_allclose(np.asarray(manager.problem.c_ineq)[:25], 25.0 - gap)
    rows = [item for item in manager.problem.constraint_data if item["name"].startswith(name)]
    assert len(rows) == 25
    assert all(item["normalize"] is None and item["target"] == 25.0 for item in rows)


def test_adding_a_scalar_constraint_updates_solver_dimensions_and_residual(cycle_manager):
    manager = cycle_manager
    manager.set_constraint("$energy_analysis.cycle_efficiency", type=">", value=0.4, normalize=True)
    efficiency = float(manager.problem.cycle_data["energy_analysis"]["cycle_efficiency"])
    assert manager.problem.get_nic() == 76
    residual = manager.problem.fitness(manager.problem.x0)
    assert residual.shape == (77,)
    assert residual[-1] == pytest.approx((0.4 - efficiency) / 0.4)
    assert manager.solver.problem is manager.problem
    # The re-created solver must also know about the extra inequality.
    assert manager.solver.N_ineq == 76


def test_saved_configuration_reloads_current_design_with_same_cycle_outputs(cycle_manager, tmp_path):
    problem = cycle_manager.problem
    point = nearby_design(problem)
    before = snapshot(problem, point)
    path = tmp_path / "current_configuration.yaml"
    problem.save_current_configuration(str(path))
    saved = th.read_configuration_file(str(path))

    # Save the evaluated design, not the initial values in the source fixture.
    np.testing.assert_allclose(
        [saved["design_variables"][name]["value"] for name in problem.variable_names],
        point, rtol=0, atol=1e-9,
    )
    reloaded = th.ThermodynamicCycleProblem(saved, out_dir=str(tmp_path / "reloaded"))
    assert reloaded.variable_names == problem.variable_names
    assert_snapshot(snapshot(reloaded, reloaded.x0), before)
    np.testing.assert_allclose(reloaded.get_bounds(), problem.get_bounds(), rtol=0, atol=1e-9)


def test_excel_export_roundtrips_all_endpoint_properties_and_energy(cycle_manager, tmp_path):
    problem = cycle_manager.problem
    point = nearby_design(problem)
    problem.fitness(point)
    path = tmp_path / "cycle.xlsx"
    problem.save_data_to_excel(str(path))
    workbook = pd.read_excel(path, sheet_name=None, engine="openpyxl")
    assert set(workbook) == {"cycle_states", "energy_analysis"}
    states = workbook["cycle_states"].set_index("state")
    assert states.index.is_unique
    properties = {
        "temperature": ("T", "K"), "pressure": ("p", "Pa"),
        "density": ("rho", "kg/m3"), "quality": ("Q", "-"),
        "compressibility_factor": ("Z", "-"), "internal_energy": ("u", "J/kg"),
        "enthalpy": ("h", "J/kg"), "entropy": ("s", "J/kg/K"),
        "isobaric_heat_capacity": ("cp", "J/kg/K"),
        "isochoric_heat_capacity": ("cv", "J/kg/K"),
        "heat_capacity_ratio": ("gamma", "-"), "speed_of_sound": ("a", "m/s"),
        "dynamic_viscosity": ("mu", "Pa*s"), "thermal_conductivity": ("k", "W/m/K"),
        "superheating_degree": ("superheating", "K"),
        "subcooling_degree": ("subcooling", "K"),
    }
    assert set(states.columns) == {"fluid_name", *properties}
    assert states.loc["units", "fluid_name"] == "-"
    for column, (_, unit) in properties.items():
        assert states.loc["units", column] == unit

    expected_names = {"units"}
    for name, component in problem.cycle_data["components"].items():
        streams = ((f"{name}_{side}", component[side]) for side in ["hot_side", "cold_side"]) if component["type"] == "heat_exchanger" else [(name, component)]
        for stream_name, stream in streams:
            for suffix, endpoint in [("in", "state_in"), ("out", "state_out")]:
                row_name = f"{stream_name}_{suffix}"
                expected_names.add(row_name)
                expected = stream[endpoint]
                assert states.loc[row_name, "fluid_name"] == expected["fluid_name"]
                for column, (alias, _) in properties.items():
                    np.testing.assert_allclose(
                        float(states.loc[row_name, column]), float(expected[alias]),
                        rtol=1e-12, atol=1e-9, equal_nan=True,
                        err_msg=f"{row_name}/{column}",
                    )
    assert set(states.index) == expected_names
    energy = workbook["energy_analysis"].set_index("Parameter")["Value"]
    assert set(energy.index) == set(problem.cycle_data["energy_analysis"])
    for name, expected in problem.cycle_data["energy_analysis"].items():
        assert float(energy[name]) == pytest.approx(float(expected), rel=1e-12, abs=1e-8)


@pytest.mark.plotting
def test_iteration_callbacks_export_the_evaluated_design(cycle_manager):
    manager = cycle_manager
    point = nearby_design(manager.problem)
    normalized = manager.problem.scale_physical_to_normalized(point)
    # Evaluate through the solver adapter, as occurs before iteration observers.
    manager.solver.fitness(normalized)
    manager.save_config_callback(normalized, 4)
    saved = th.read_configuration_file(str(Path(manager.optimization_dir) / "iteration_004.yaml"))
    np.testing.assert_allclose(
        [saved["design_variables"][name]["value"] for name in manager.problem.variable_names],
        point, rtol=0, atol=1e-8,
    )

    manager.plot_cycle_callback(normalized, 4)
    assert manager.problem.figure._suptitle.get_text() == "Optimization iteration: 004"
    manager.problem.figure.set_size_inches(3, 2)
    manager.save_plot_callback(normalized, 4)
    image = Path(manager.optimization_dir) / "iteration_004.png"
    assert image.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert np.std(plt.imread(image)[..., :3]) > 0.01

    manager.save_report_callback(normalized, 4)
    text = (Path(manager.optimization_dir) / "iteration_004.txt").read_text(encoding="utf-8")
    assert "compressor_inlet_pressure" in text
    assert "heater.temperature_difference" in text
    assert "KKT" in text
