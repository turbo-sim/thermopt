"""Exercise actual JaxProp models through unchanged eager cycle/component code."""

import copy
import warnings

import jaxprop as jp
import numpy as np
import openpyxl
import pytest

import thermopt as th
from thermopt import utilities
from thermopt.components import compression_process, expansion_process, heat_exchanger
from regression_support import load_case


@pytest.fixture(scope="module")
def table_options(test_workspace):
    return dict(
        h_min=np.float64(110e3),
        h_max=1.3e6,
        p_min=8e6,
        p_max=39e6,
        N_h=np.int64(32),
        N_p=32,
        N_p_sat=4,
        coarse_step=1,
        table_dir=str(test_workspace / "bicubic_tables"),
    )


def test_problem_uses_native_bicubic_and_reuses_models(
    table_options, tmp_path, monkeypatch
):
    config = copy.deepcopy(load_case("sco2_recuperated")["problem_formulation"])
    config["fixed_parameters"]["working_fluid"].update(
        method="bicubic", table_options=table_options
    )
    with warnings.catch_warnings(record=True) as messages:
        warnings.simplefilter("always")
        problem = th.ThermodynamicCycleProblem(config, out_dir=str(tmp_path))
    assert not any(
        "array is being set as static" in str(message.message) for message in messages
    )
    assert type(problem.fluid) is jp.FluidBicubic
    assert type(problem.reference_fluid) is jp.FluidCoolProp
    assert (
        problem.params["working_fluid"]["gas_at_maximum_temperature"].p
        < table_options["p_min"]
    )
    identities = {role: id(fluid) for role, fluid in problem.fluids.items()}

    def unexpected_construction(*args, **kwargs):
        raise AssertionError("configured cycle evaluation must reuse its fluid objects")

    monkeypatch.setattr(utilities, "fluid_from_configuration", unexpected_construction)
    first = np.asarray(problem.fitness(problem.x0.copy()))
    repeated = np.asarray(problem.fitness(problem.x0.copy()))
    np.testing.assert_array_equal(first, repeated)
    assert np.all(np.isfinite(first))
    assert identities == {role: id(fluid) for role, fluid in problem.fluids.items()}
    assert (
        problem.cycle_data["components"]["compressor"]["state_out"].provider
        == "bicubic"
    )
    for component in problem.cycle_data["components"].values():
        if component["type"] == "heat_exchanger":
            np.testing.assert_allclose(
                component["temperature_difference"],
                component["hot_side"]["states"].T - component["cold_side"]["states"].T,
            )


@pytest.mark.parametrize(
    "process,p_in,t_in,p_out",
    [
        (compression_process, 9e6, 400.0, 18e6),
        (expansion_process, 18e6, 600.0, 9e6),
    ],
)
def test_native_bicubic_components_match_source(
    process, p_in, t_in, p_out, table_options
):
    source = jp.Fluid("CO2")
    table = jp.Fluid("CO2", method="bicubic", table_options=table_options)
    inlet = source.get_state(jp.PT_INPUTS, p_in, t_in)
    exact = process(source, float(inlet.h), p_in, p_out, efficiency=0.85)
    actual = process(table, float(inlet.h), p_in, p_out, efficiency=0.85)
    assert actual["state_out"].provider == "bicubic"
    np.testing.assert_allclose(
        actual["state_out"].T, exact["state_out"].T, atol=0.3, rtol=0
    )
    np.testing.assert_allclose(
        actual["specific_work"], exact["specific_work"], rtol=0.005
    )
    with pytest.raises(jp.OutOfDomainError):
        process(table, 2e6, p_in, p_out, efficiency=0.85)
    with pytest.raises(jp.UnsupportedOperationError, match="saturation"):
        heat_exchanger(
            table,
            900e3,
            800e3,
            18e6,
            18e6,
            table,
            400e3,
            500e3,
            9e6,
            9e6,
            include_saturation_nodes=True,
        )


def test_case_model_selection_is_explicit():
    assert (
        type(
            utilities.fluid_from_configuration(
                dict(name="Air", method="coolprop", backend="HEOS")
            )
        )
        is jp.FluidCoolProp
    )


def test_bicubic_excel_contains_numbers_and_blank_unsupported_fields(
    table_options, tmp_path
):
    config = copy.deepcopy(load_case("sco2_recuperated")["problem_formulation"])
    config["fixed_parameters"]["working_fluid"].update(
        method="bicubic", table_options=table_options
    )
    problem = th.ThermodynamicCycleProblem(config, out_dir=str(tmp_path))
    problem.fitness(problem.x0)
    filename = tmp_path / "native_bicubic.xlsx"
    problem.save_data_to_excel(str(filename))
    workbook = openpyxl.load_workbook(filename, data_only=True)
    sheet = workbook["cycle_states"]
    columns = {cell.value: cell.column for cell in sheet[1]}
    row = next(
        row for row in sheet.iter_rows(min_row=3) if row[0].value == "compressor_out"
    )
    state = problem.cycle_data["components"]["compressor"]["state_out"]
    for column, name in (("temperature", "T"), ("pressure", "p"), ("enthalpy", "h")):
        cell = row[columns[column] - 1]
        assert cell.data_type == "n"
        np.testing.assert_allclose(cell.value, float(state[name]), rtol=1e-12)
    assert row[columns["superheating_degree"] - 1].value is None
    assert row[columns["subcooling_degree"] - 1].value is None
