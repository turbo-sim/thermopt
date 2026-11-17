"""Exercise watcher control/file updates headlessly, with a real input thread."""

import threading
from pathlib import Path

import numpy as np
import pytest
import yaml

import thermopt as th
from regression_support import load_case


def live_problem(tmp_path):
    config = load_case("sco2_recuperated")
    config["problem_formulation"]["plot_settings"].update(
        diagrams=[dict(x_prop="s", y_prop="T")],
        pinch_point_diagram=False,
        fluid=dict(plot_saturation_line=False, plot_critical_point=False),
    )
    source = tmp_path / "live.yaml"
    source.write_text(
        yaml.safe_dump(
            th.convert_numpy_to_python(config, precision=None), sort_keys=False
        ),
        encoding="utf-8",
    )
    problem = th.ThermodynamicCycleProblem(
        config["problem_formulation"], out_dir=str(tmp_path)
    )
    return config, source, problem


@pytest.mark.plotting
def test_watcher_reloads_a_design_writes_requested_report_and_pickles_after_exit(
    tmp_path, monkeypatch
):
    config, source, problem = live_problem(tmp_path)
    stop = threading.Event()
    monkeypatch.setattr("builtins.input", lambda: stop.wait(timeout=5))
    original_plot = problem.plot_cycle
    seen = []

    def observe():
        original_plot()
        seen.append(
            float(
                problem.configuration["design_variables"][
                    "heat_source_exit_temperature"
                ]["value"]
            )
        )
        if len(seen) == 1:
            config["problem_formulation"]["design_variables"][
                "heat_source_exit_temperature"
            ]["value"] += 1
            source.write_text(
                yaml.safe_dump(
                    th.convert_numpy_to_python(config, precision=None), sort_keys=False
                ),
                encoding="utf-8",
            )
        else:
            stop.set()

    with monkeypatch.context() as patch:
        patch.setattr(problem, "plot_cycle", observe)
        problem.plot_cycle_realtime(
            str(source), update_interval=0.001, write_report=True
        )
    stop.set()
    problem.input_thread.join(timeout=1)
    assert not problem.input_thread.is_alive()
    assert seen[-1] == seen[0] + 1
    report = (tmp_path / "initial_guess_report.txt").read_text(encoding="utf-8")
    assert (
        "heat_source_exit_temperature" in report
        and "heater.temperature_difference" in report
    )
    filename = th.save_to_pickle(
        problem, tmp_path / "after_watcher.pkl", timestamp=False
    )
    restored = th.load_from_pickle(filename)
    assert not hasattr(restored, "input_thread")
    np.testing.assert_allclose(
        restored.fitness(restored.x0), problem.fitness(problem.x0)
    )


@pytest.mark.plotting
def test_watcher_eof_still_produces_an_exportable_initial_figure(tmp_path, monkeypatch):
    _, source, problem = live_problem(tmp_path)
    import matplotlib.pyplot as plt

    original_plot = problem.plot_cycle
    draws = []

    def bounded_draw():
        original_plot()
        draws.append(1)
        if len(draws) >= 3:
            plt.close(problem.figure)

    def eof():
        raise EOFError

    monkeypatch.setattr("builtins.input", eof)
    with monkeypatch.context() as patch:
        patch.setattr(problem, "plot_cycle", bounded_draw)
        problem.plot_cycle_realtime(str(source), update_interval=0.001)
    problem.input_thread.join(timeout=1)
    assert problem.figure is not None
    assert len(draws) == 1
    problem.figure.savefig(tmp_path / "initial.png")
    assert (tmp_path / "initial.png").read_bytes().startswith(b"\x89PNG")
