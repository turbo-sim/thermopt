"""Evaluate migrated study helpers in isolated interpreter namespaces."""

from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize(
    "folder,comparison,seed",
    [
        ("examples/rankine_transcritical_vs_trilateral", True, "independent"),
        ("examples/rankine_wet2dry/scripts", False, "independent"),
        ("examples/rankine_wet2dry/scripts", False, "saved"),
    ],
)
def test_study_configuration_and_property_evaluation(
    folder, comparison, seed, tmp_path
):
    directory = ROOT / folder
    code = """
import sys
from copy import deepcopy
from pathlib import Path
import numpy as np
import thermopt as th
sys.path.insert(0, sys.argv[1])
if sys.argv[3] == 'True':
    import run_sensitivity as run
    from fast_evaluator import FastProblem
    config = run.configuration('butane', 50, .8, .8)
    native = run.ComparisonProblem(deepcopy(config['problem_formulation']), out_dir=sys.argv[2])
    fast = FastProblem(deepcopy(config['problem_formulation']), out_dir=sys.argv[2])
    np.testing.assert_allclose(fast.fitness(fast.x0), native.fitness(native.x0), rtol=2e-6, atol=2e-6)
    assert run.reference_heat() > 0
else:
    import study as run
    config = run.configuration('butane', 'partial', 100, .8, seed=sys.argv[4])
    for settings in config['problem_formulation']['fixed_parameters'].values():
        if isinstance(settings, dict) and 'name' in settings and 'backend' in settings:
            assert settings['method'] == 'coolprop'
            assert 'exceptions' not in settings
    problem = th.ThermodynamicCycleProblem(config['problem_formulation'], out_dir=sys.argv[2])
    assert np.all(np.isfinite(problem.fitness(problem.x0)))
    assert run.fixed_reference_heat_J_per_kg() > 0
"""
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            code,
            str(directory),
            str(tmp_path),
            str(comparison),
            seed,
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
