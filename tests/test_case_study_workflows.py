"""Protect native study ownership, saved-input boundaries and cache provenance."""

from copy import deepcopy
import importlib.util
from pathlib import Path
import subprocess
import sys
import json
import runpy
from types import SimpleNamespace

import numpy as np
import pytest
import yaml

import thermopt as th
from case_studies import study_io
from regression_support import load_case

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("complete", [False, True])
def test_study_reader_selects_the_documented_formulation(tmp_path, complete):
    config = load_case("orc_subcritical")
    payload = config if complete else config["problem_formulation"]
    source = tmp_path / "saved.yaml"
    source.write_text(
        yaml.safe_dump(
            th.convert_numpy_to_python(payload, precision=None), sort_keys=False
        ),
        encoding="utf-8",
    )
    result = study_io.read_formulation(source)
    assert list(result["design_variables"]) == list(
        config["problem_formulation"]["design_variables"]
    )
    np.testing.assert_equal(result, config["problem_formulation"])


def test_study_cache_identity_tracks_code_dependencies_and_source_references(
    tmp_path, monkeypatch
):
    config = load_case("orc_subcritical")
    code = tmp_path / "study.py"
    code.write_text("# first source\n", encoding="utf-8")
    initial = study_io.calculation_identity(config, [code])
    assert study_io.calculation_identity(config, [code]) == initial
    code.write_text("# changed source\n", encoding="utf-8")
    assert study_io.calculation_identity(config, [code]) != initial
    code.write_text("# first source\n", encoding="utf-8")
    version = study_io.version
    monkeypatch.setattr(
        study_io,
        "version",
        lambda name: "future.test" if name == "jaxprop" else version(name),
    )
    assert study_io.calculation_identity(config, [code]) != initial
    monkeypatch.setattr(study_io, "version", version)
    signature = study_io.source_signature
    monkeypatch.setattr(
        study_io,
        "source_signature",
        lambda *args: dict(signature(*args), test_reference_revision=1),
    )
    assert study_io.calculation_identity(config, [code]) != initial


@pytest.mark.parametrize("family", ["transcritical", "wet2dry"])
def test_study_models_caches_and_stale_results_are_owned_and_setup_only(
    family, tmp_path
):
    code = r"""
from copy import deepcopy
from pathlib import Path
import json
import sys
from unittest.mock import patch
import numpy as np
import thermopt as th
sys.path.insert(0, sys.argv[1])
directory = Path(sys.argv[2])
if sys.argv[3] == 'transcritical':
    import run_sensitivity as run
    from fast_evaluator import FastProblem
    config = run.configuration('butane', 50, .8, .8)
    first = FastProblem(deepcopy(config['problem_formulation']), out_dir=str(directory/'first'))
    second = FastProblem(deepcopy(config['problem_formulation']), out_dir=str(directory/'second'))
    baseline = first.fitness(first.x0)
    assert first._water_pt_cache and first._water_temperature_cache
    assert first._water_pt_cache is not second._water_pt_cache
    water = first.fluids['heating_fluid']
    assert water is not second.fluids['heating_fluid']
    with patch('thermopt.cycles.setup.source_signature', side_effect=AssertionError('hot source probe')), patch('thermopt.utilities.fluid_from_configuration', side_effect=AssertionError('hot construction')):
        np.testing.assert_array_equal(first.fitness(first.x0), baseline)
    cache = first._water_pt_cache
    first.update_problem(deepcopy(config['problem_formulation']))
    assert first._water_pt_cache is cache
    changed = deepcopy(config['problem_formulation'])
    changed['fixed_parameters']['heating_fluid']['name'] = 'water' if water.name != 'water' else 'Water'
    first.update_problem(changed)
    assert first.fluids['heating_fluid'] is not water
    assert first._water_pt_cache == {} and first._water_pt_cache is not cache
    point = run.case_name('butane', 50, .8, .8)
    invoke = lambda: run.solve_point('butane', 50, .8, .8)
else:
    import study as run
    import yaml
    config = run.configuration('butane', 'partial', 100, .8, seed='independent')
    first = th.ThermodynamicCycleProblem(deepcopy(config['problem_formulation']), out_dir=str(directory/'first'))
    second = th.ThermodynamicCycleProblem(deepcopy(config['problem_formulation']), out_dir=str(directory/'second'))
    assert first.fluids['working_fluid'] is not second.fluids['working_fluid']
    baseline = first.fitness(first.x0)
    with patch('thermopt.cycles.setup.source_signature', side_effect=AssertionError('hot source probe')), patch('thermopt.utilities.fluid_from_configuration', side_effect=AssertionError('hot construction')):
        np.testing.assert_array_equal(first.fitness(first.x0), baseline)
    original_initializers = run.INITIALIZERS
    run.INITIALIZERS = directory/'initializers'
    initializer = run.INITIALIZERS/'sensitivity_butane'/run.point_name('partial',100,.8)/'optimal_solution.yaml'
    initializer.parent.mkdir(parents=True)
    initializer.write_text(yaml.safe_dump({'design_variables':{'heat_source_exit_temperature':{'value':380.}}}))
    first_seed = run.configuration('butane','partial',100,.8,seed='saved')
    initializer.write_text(yaml.safe_dump({'design_variables':{'heat_source_exit_temperature':{'value':381.}}}))
    second_seed = run.configuration('butane','partial',100,.8,seed='saved')
    assert run.signature(first_seed) != run.signature(second_seed)
    run.INITIALIZERS = original_initializers
    point = 'sensitivity_butane/' + run.point_name('partial', 100, .8)
    invoke = lambda: run.solve_point('butane', 'partial', 100, .8)
signature = run.signature(config)
initial_only = deepcopy(config)
next(iter(initial_only['problem_formulation']['design_variables'].values()))['value'] += 1
assert run.signature(initial_only) == signature
solver_change = deepcopy(config)
solver_change['solver_options']['tolerance'] = 1e-8
assert run.signature(solver_change) != signature
run.RESULTS = directory/'fresh_results'
target = run.RESULTS/point
target.mkdir(parents=True)
(target/'result.json').write_text(json.dumps({'signature':'historical-signature'}))
try:
    invoke()
except RuntimeError as error:
    assert 'directory' in str(error)
else:
    raise AssertionError('a stale result must not be silently reused')
from pysolver_view import DomainError
run.RESULTS = directory/'domain_failure'
run.MAX_ATTEMPTS = 1
factory = 'Cycle' if sys.argv[3] == 'transcritical' else 'create_native_cycle'
with patch.object(run, factory, side_effect=DomainError('invalid seed')):
    result = invoke()
assert result['row']['status'] == 'unconverged'
assert not result['row']['success'] and result['row']['attempts'] == 1
import jaxprop as props
for index, error_type in enumerate((props.InvalidConfigurationError, props.UnsupportedInputPairError)):
    run.RESULTS = directory/f'setup_failure_{index}'
    failure = error_type(request='study setup', reason='unsupported setup probe')
    with patch.object(run, factory, side_effect=failure):
        try:
            invoke()
        except error_type as error:
            assert error is failure
        else:
            raise AssertionError('unsupported/configuration errors must not become failed starts')
"""
    folder = ROOT / (
        "case_studies/rankine_transcritical_vs_trilateral"
        if family == "transcritical"
        else "case_studies/rankine_wet2dry/scripts"
    )
    result = subprocess.run(
        [sys.executable, "-c", code, str(folder), str(tmp_path), family],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=90,
    )
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize(
    "source,variant",
    [
        (
            "sCO2_recuperated/case_sCO2_recuperated.yaml",
            "sCO2_recuperated/case_sCO2_recuperated_bicubic.yaml",
        ),
        (
            "sCO2_comparison/case_sCO2_simple.yaml",
            "sCO2_comparison/case_sCO2_simple_bicubic.yaml",
        ),
        (
            "refrigeration_air/case_air_refrigeration_recuperated.yaml",
            "refrigeration_air/case_air_refrigeration_recuperated_bicubic.yaml",
        ),
        (
            "sCO2_recompression/case_sCO2_recompression_coolprop.yaml",
            "sCO2_recompression/case_sCO2_recompression.yaml",
        ),
    ],
)
def test_explicit_native_table_cases_preserve_source_physics(source, variant):
    directory = ROOT / "case_studies"
    exact = th.read_configuration_file(directory / source)
    table = th.read_configuration_file(directory / variant)
    settings = table["problem_formulation"]["fixed_parameters"]["working_fluid"]
    assert settings["method"] == "bicubic"
    assert settings["table_options"]["N_h"] == settings["table_options"]["N_p"] == 200
    assert "results" in settings["table_options"]["table_dir"]
    settings["method"] = "coolprop"
    settings.pop("table_options")
    np.testing.assert_equal(table, exact)


def test_ptes_saved_input_table_programs_use_current_native_exports(
    tmp_path, monkeypatch
):
    formulation = load_case("ptes")["problem_formulation"]
    problem = th.ThermodynamicCycleProblem(
        formulation, out_dir=str(tmp_path / "evaluation")
    )
    problem.fitness(problem.x0)
    inputs, outputs = tmp_path / "inputs", tmp_path / "results"
    inputs.mkdir()
    problem.save_data_to_excel(inputs / "optimal_solution.xlsx")
    th.save_to_pickle(
        SimpleNamespace(problem=problem), str(inputs / "solver.pkl"), timestamp=False
    )
    monkeypatch.setenv("THERMOPT_INPUT_DIR", str(inputs))
    monkeypatch.setenv("THERMOPT_RESULTS_DIR", str(outputs))
    directory = ROOT / "case_studies/PTES_CO2/case_baseline"
    for filename in (
        "create_energy_analysis_table.py",
        "create_thermodynamic_states_table.py",
    ):
        runpy.run_path(str(directory / filename), run_name="__main__")
    energy = (outputs / "table_energy_analysis.txt").read_text(encoding="utf-8")
    states = (outputs / "component_states_summary.txt").read_text(encoding="utf-8")
    expected = float(problem.cycle_data["energy_analysis"]["roundtrip_efficiency"])
    assert f"{expected:.3f}" in energy and "PTES roundtrip efficiency" in energy
    assert "Compressor in" in states and "Pressure (bar)" in states


def test_case_inventory_covers_every_owned_program_and_declares_scope():
    from test_property_boundary import project_sources

    recorded = json.loads((ROOT / "case_studies/EXECUTION_INVENTORY.json").read_text())
    expected = {
        path.relative_to(ROOT / "case_studies").as_posix()
        for path in project_sources(ROOT / "case_studies")
        if path.suffix == ".py"
    }
    assert {item["path"] for item in recorded["programs"]} == expected
    for item in recorded["programs"]:
        assert item["category"] and item["execution_budget"] and item["acceptance"]
        if item["category"] == "historical-figure-preparation":
            assert item["execution_budget"] == "No automatic run"
        if item["category"] == "saved-results-analysis":
            assert item["inputs"]
