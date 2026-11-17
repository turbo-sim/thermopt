"""Check native-state outputs and rerunnable snapshots at the I/O boundary."""

from copy import deepcopy
import json
from pathlib import Path

import jax
import jax.numpy as jnp
import jaxprop as jp
import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd
import pysolver_view as psv
import pytest
import yaml

import thermopt as th
from regression_support import load_case, nearby_design, snapshot, assert_snapshot

jax.config.update("jax_enable_x64", True)


@pytest.fixture(scope="module")
def output_table_options(test_workspace):
    return dict(
        h_min=110e3,
        h_max=1.3e6,
        p_min=8e6,
        p_max=39e6,
        N_h=32,
        N_p=32,
        N_p_sat=4,
        coarse_step=1,
        table_dir=str(test_workspace / "output_tables"),
    )


@pytest.fixture(params=["coolprop", "bicubic"])
def output_manager(request, output_table_options, tmp_path):
    config = load_case("sco2_recuperated")
    config["solver_options"].update(
        library="pysolver",
        callbacks={},
        print_convergence=False,
        plot_convergence=False,
        update_on="iteration",
        max_iterations=2,
        extra_options=dict(first_step_fraction=0.01, optimality_tol=1e-4),
    )
    formulation = config["problem_formulation"]
    formulation["plot_settings"].update(
        fluid=dict(plot_saturation_line=False, plot_critical_point=False),
        diagrams=[dict(x_prop="s", y_prop="T"), dict(x_prop="h", y_prop="p")],
        pinch_point_diagram=True,
    )
    if request.param == "bicubic":
        formulation["fixed_parameters"]["working_fluid"].update(
            method="bicubic",
            table_options=output_table_options,
        )
    source = tmp_path / "input.yaml"
    source.write_text(
        yaml.safe_dump(
            th.convert_numpy_to_python(config, precision=14), sort_keys=False
        ),
        encoding="utf-8",
    )
    return th.ThermodynamicCycleOptimization(source, out_dir=str(tmp_path / "results"))


def test_manager_snapshot_reloads_complete_configuration_and_current_design(
    output_manager, tmp_path
):
    manager = output_manager
    point = nearby_design(manager.problem)
    before = snapshot(manager.problem, point)
    original = deepcopy(manager.config)
    # A reporting precision must never round a solver tolerance to zero.
    manager.config["solver_options"]["tolerance"] = np.float64(1e-14)
    filename = tmp_path / "rerunnable.yaml"
    manager.save_current_configuration(filename)
    saved = th.read_configuration_file(filename)
    assert set(saved) == {"solver_options", "problem_formulation"}
    assert float(saved["solver_options"]["tolerance"]) == 1e-14
    formulation = saved["problem_formulation"]
    assert list(formulation["design_variables"]) == manager.problem.variable_names
    np.testing.assert_allclose(
        [entry["value"] for entry in formulation["design_variables"].values()],
        point,
        rtol=0,
        atol=1e-12,
    )
    np.testing.assert_equal(
        formulation["fixed_parameters"]["working_fluid"],
        original["problem_formulation"]["fixed_parameters"]["working_fluid"],
    )
    assert "critical_point" not in formulation["fixed_parameters"]["working_fluid"]
    reloaded = th.ThermodynamicCycleOptimization(
        filename, out_dir=str(tmp_path / "reloaded")
    )
    assert_snapshot(snapshot(reloaded.problem, reloaded.problem.x0), before)
    np.testing.assert_allclose(
        reloaded.problem.get_bounds(), manager.problem.get_bounds()
    )
    assert manager.problem.fluid.reference_id == reloaded.problem.fluid.reference_id
    # Saving must not replace the manager's configured initial design.
    np.testing.assert_equal(
        manager.config["problem_formulation"], original["problem_formulation"]
    )


def test_excel_has_numeric_energy_cells_and_never_evaluates_properties(
    output_manager, tmp_path, monkeypatch
):
    problem = output_manager.problem
    problem.fitness(nearby_design(problem))
    expected = dict(problem.cycle_data["energy_analysis"])

    def unexpected(*args, **kwargs):
        raise AssertionError("output must consume returned states without new flashes")

    monkeypatch.setattr(jp.FluidCoolProp, "get_state", unexpected)
    monkeypatch.setattr(jp.FluidBicubic, "get_state", unexpected)
    filename = tmp_path / "native.xlsx"
    problem.save_data_to_excel(filename)
    book = openpyxl.load_workbook(filename, data_only=True)
    try:
        for name, value in book["energy_analysis"].iter_rows(min_row=2):
            assert value.data_type == "n", name.value
            assert value.value == pytest.approx(
                float(expected[name.value]), rel=1e-12, abs=1e-10
            )
        sheet = book["cycle_states"]
        headers = {cell.value: cell.column for cell in sheet[1]}
        compressor = next(
            row
            for row in sheet.iter_rows(min_row=3)
            if row[0].value == "compressor_out"
        )
        assert compressor[headers["temperature"] - 1].data_type == "n"
        if problem.fluid.provider == "bicubic":
            assert compressor[headers["superheating_degree"] - 1].value is None
    finally:
        book.close()


@pytest.mark.plotting
def test_native_cycle_plots_reuse_artists_and_source_backgrounds(
    output_manager, monkeypatch
):
    problem = output_manager.problem
    problem.plot_settings["fluid"].update(
        plot_saturation_line=True, plot_critical_point=True
    )

    def draw():
        with monkeypatch.context() as patch:
            if problem.fluid.provider == "bicubic":
                # Saturation backgrounds extend far below this table's 8 MPa floor.
                # JaxProp must use its source EOS rather than extrapolate the table.
                def unexpected(*args, **kwargs):
                    raise AssertionError(
                        "phase backgrounds must not flash through the table"
                    )

                patch.setattr(jp.FluidBicubic, "get_state", unexpected)
            problem.plot_cycle()
        for index, (x, y) in enumerate([("s", "T"), ("h", "p")]):
            states = problem.cycle_data["components"]["compressor"]["states"]
            line = problem.graphics["process_lines"][index]["compressor"]
            np.testing.assert_allclose(line.get_xdata(), states[x])
            np.testing.assert_allclose(line.get_ydata(), states[y])
        for name, artists in problem.graphics["pinch_point_lines"][2].items():
            hot, cold = artists["hot_line"], artists["cold_line"]
            np.testing.assert_allclose(hot.get_xdata(), cold.get_xdata(), atol=1e-6)
            np.testing.assert_allclose(
                np.asarray(hot.get_ydata()) - cold.get_ydata(),
                problem.cycle_data["components"][name]["temperature_difference"],
            )

    draw()
    figure = problem.figure
    artist = problem.graphics["process_lines"][0]["compressor"]
    counts = [len(axis.lines) for axis in figure.axes]
    problem.fitness(nearby_design(problem))
    draw()
    assert problem.figure is figure
    assert problem.graphics["process_lines"][0]["compressor"] is artist
    assert [len(axis.lines) for axis in figure.axes] == counts
    plt.close(figure)
    draw()
    assert problem.figure is not figure
    assert problem.graphics["process_lines"][0]["compressor"] is not artist


@pytest.mark.optimization
def test_final_outputs_restore_solver_design_after_reporting_probe(
    output_manager, tmp_path
):
    manager = output_manager
    manager.run_optimization()
    final = manager.solver.x_final.copy()
    expected = snapshot(manager.problem, final)
    manager.problem.fitness(nearby_design(manager.problem))
    manager.save_results()
    filename = Path(manager.out_dir) / "optimal_solution.yaml"
    restored = th.ThermodynamicCycleOptimization(
        filename, out_dir=str(tmp_path / "reloaded")
    )
    np.testing.assert_allclose(restored.problem.x0, final, rtol=0, atol=1e-12)
    assert_snapshot(snapshot(restored.problem, restored.problem.x0), expected)
    workbook = pd.read_excel(
        Path(manager.out_dir) / "optimal_solution.xlsx", sheet_name="energy_analysis"
    )
    energy = workbook.set_index("Parameter")["Value"]
    assert float(energy["cycle_efficiency"]) == pytest.approx(
        float(manager.problem.cycle_data["energy_analysis"]["cycle_efficiency"]),
        rel=1e-12,
    )


@pytest.mark.parametrize("storage", [np.asarray, jnp.asarray])
def test_numeric_output_conversion_preserves_precision_and_json_csv_types(
    storage, tmp_path
):
    payload = dict(
        temperature=storage(300.125),
        enabled=storage(True),
        steps=storage(3),
        curve=storage([1.0, 2.0]),
        tolerance=storage(1e-14),
        unavailable=None,
    )
    native = th.convert_numpy_to_python(payload, precision=None)
    assert native == dict(
        temperature=300.125,
        enabled=True,
        steps=3,
        curve=[1.0, 2.0],
        tolerance=1e-14,
        unavailable=None,
    )
    assert type(native["steps"]) is int and type(native["enabled"]) is bool
    assert json.loads(json.dumps(native, allow_nan=False)) == native
    table = pd.DataFrame(
        [
            dict(
                temperature=native["temperature"],
                tolerance=native["tolerance"],
                unavailable=native["unavailable"],
            )
        ]
    )
    filename = tmp_path / "numeric.csv"
    table.to_csv(filename, index=False)
    restored = pd.read_csv(filename)
    assert restored.loc[0, "tolerance"] == pytest.approx(1e-14, rel=1e-12, abs=0)
    assert restored.loc[0, "temperature"] == 300.125
    assert pd.isna(restored.loc[0, "unavailable"])


@pytest.mark.parametrize(
    "error_type",
    [
        jp.InvalidStateError,
        jp.OutOfDomainError,
        jp.CalculationError,
        jp.ConvergenceError,
        jp.UnsupportedPropertyError,
        jp.InvalidConfigurationError,
        ValueError,
        KeyError,
    ],
)
def test_solver_boundary_translates_only_recognized_domain_failures(
    output_manager, monkeypatch, error_type
):
    import thermopt.optimize_cycle as module

    problem = output_manager.problem
    if issubclass(error_type, jp.JaxPropError):
        failure = error_type(
            request="HmassP_INPUTS",
            reason="controlled test failure",
            provider=problem.fluid.provider,
            fluid_name=problem.fluid.name,
        )
    else:
        failure = error_type("application defect")

    def failed_evaluation(*args, **kwargs):
        raise failure

    monkeypatch.setitem(
        module.CYCLE_TOPOLOGIES, problem.cycle_topology, failed_evaluation
    )
    with pytest.raises(error_type) as direct:
        problem.fitness(problem.x0)
    assert direct.value is failure
    normalized = problem.scale_physical_to_normalized(problem.x0)
    recoverable = isinstance(failure, (jp.InvalidStateError, jp.CalculationError))
    with pytest.raises(psv.DomainError if recoverable else error_type) as solver_error:
        output_manager.solver.fitness(normalized)
    if recoverable:
        assert solver_error.value.__cause__ is failure
        assert solver_error.value.metadata["provider"] == problem.fluid.provider
        assert solver_error.value.metadata["request"] == "HmassP_INPUTS"
    else:
        assert solver_error.value is failure


def test_real_invalid_trial_retains_native_failure_and_marks_solver_domain(
    output_manager,
):
    manager = output_manager
    problem = manager.problem
    point = problem.x0.copy()
    point[problem.variable_names.index("compressor_inlet_pressure")] = -1.0
    with pytest.raises((jp.InvalidStateError, jp.CalculationError)) as direct:
        problem.fitness(point)
    with pytest.raises(psv.DomainError) as translated:
        manager.solver.fitness(problem.scale_physical_to_normalized(point))
    assert isinstance(translated.value.__cause__, type(direct.value))
    assert translated.value.metadata["provider"] == problem.fluid.provider


def test_unsolved_manager_rejects_final_export_without_writing_files(output_manager):
    manager = output_manager
    with pytest.raises(RuntimeError, match="Run optimization"):
        manager.save_results()
    assert list(Path(manager.out_dir).glob("optimal_solution.*")) == []


def test_iteration_outputs_restore_the_observed_design(output_manager):
    manager = output_manager
    point = nearby_design(manager.problem)
    normalized = manager.problem.scale_physical_to_normalized(point)
    # The callback receives solver coordinates; compare with their exact
    # physical mapping, including the final floating-point rounding step.
    point = manager.problem.scale_normalized_to_physical(normalized)
    # An optimizer can have visited a derivative probe before notifying an
    # accepted iteration. The callback argument is the design to export.
    manager.save_config_callback(normalized, 3)
    saved = th.read_configuration_file(
        Path(manager.optimization_dir) / "iteration_003.yaml"
    )
    np.testing.assert_allclose(
        [
            item["value"]
            for item in saved["problem_formulation"]["design_variables"].values()
        ],
        point,
        rtol=0,
        atol=1e-12,
    )
    manager.problem.fitness(manager.problem.x0)
    manager.plot_cycle_callback(normalized, 3)
    np.testing.assert_array_equal(list(manager.problem.x0_dict.values()), point)
    manager.save_report_callback(normalized, 3)
    np.testing.assert_array_equal(list(manager.problem.x0_dict.values()), point)


@pytest.mark.optimization
def test_save_plot_callback_creates_its_figure_when_plotting_flag_is_disabled(
    output_manager,
):
    manager = output_manager
    manager.config["solver_options"]["callbacks"] = {
        "save_plot": True,
        "plot_cycle": False,
    }
    manager.run_optimization()
    assert list(Path(manager.optimization_dir).glob("iteration_*.png"))
