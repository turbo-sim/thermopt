"""Check checkpoint eligibility and declared ratios without stopwatch assertions."""

import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location(
    "thermopt_phase_benchmark", Path(__file__).resolve().parents[1] / "benchmark/run.py"
)
protocol = importlib.util.module_from_spec(spec)
spec.loader.exec_module(protocol)


def record(timings, signature="same", passed=True):
    return dict(
        comparison_signature=signature,
        numerical_checks_passed=passed,
        timings_seconds=timings,
    )


def test_known_ratios_flag_only_cost_above_declared_limit():
    ratios, flags = protocol.compare_checkpoint(
        record({"source": 1.05, "table": 1.2}), record({"source": 1.0, "table": 1.0})
    )
    assert ratios == {"source": 1.05, "table": 1.2}
    assert flags == ["table"]


def test_exact_budget_and_improvement_are_accepted():
    _, flags = protocol.compare_checkpoint(
        record({"a": 1.10, "b": 0.8}), record({"a": 1.0, "b": 1.0})
    )
    assert flags == []


@pytest.mark.parametrize(
    "baseline",
    [record({"a": 1.0}, signature="different"), record({"a": 1.0}, passed=False)],
)
def test_changed_protocol_or_failed_numerics_cannot_be_a_reference(baseline):
    with pytest.raises(ValueError, match="incompatible"):
        protocol.compare_checkpoint(record({"a": 1.0}), baseline)


def test_missing_workload_does_not_become_a_partial_pass():
    with pytest.raises(ValueError, match="workload"):
        protocol.compare_checkpoint(record({"a": 1.0, "b": 1.0}), record({"a": 1.0}))


@pytest.mark.parametrize("value", [0.0, -1.0, float("nan"), float("inf")])
def test_invalid_timings_are_rejected(value):
    with pytest.raises(ValueError, match="finite and positive"):
        protocol.compare_checkpoint(record({"a": value}), record({"a": 1.0}))
