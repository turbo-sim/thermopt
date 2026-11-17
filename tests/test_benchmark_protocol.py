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


@pytest.mark.parametrize("requested", [0, -1, 0b10000, 0b10001])
def test_invalid_affinity_cannot_expand_or_empty_the_process_mask(requested):
    with pytest.raises(ValueError, match="nonzero subset"):
        protocol.validate_cpu_mask(requested, 0b1111)


@pytest.mark.parametrize("requested", [1, 0b1010, 0b1111])
def test_affinity_accepts_nonempty_subsets(requested):
    protocol.validate_cpu_mask(requested, 0b1111)


def test_default_affinity_observation_does_not_change_scheduling():
    assert protocol.cpu_affinity() == protocol.cpu_affinity()


@pytest.mark.parametrize("method", ["coolprop", "bicubic"])
def test_benchmark_can_select_source_or_table_from_bicubic_case_yaml(method, tmp_path):
    from benchmark.native_backends import configuration

    config, source = configuration("sco2_recompression", method, 100, False, tmp_path)
    assert source.is_file() and source.parent.parent.name == "case_studies"
    fluid = config["problem_formulation"]["fixed_parameters"]["working_fluid"]
    assert fluid["method"] == method
    if method == "coolprop":
        assert "table_options" not in fluid
    else:
        assert fluid["table_options"]["N_h"] == fluid["table_options"]["N_p"] == 100


def test_linux_affinity_keeps_the_inherited_worker_mask(monkeypatch):
    # Exercise the Linux branch on any host, including explicit inherited masks.
    from types import SimpleNamespace
    mask = {2, 5}
    def change(pid, requested):
        mask.clear()
        mask.update(requested)
    monkeypatch.setattr(protocol, "os", SimpleNamespace(name="posix", sched_getaffinity=lambda pid: set(mask), sched_setaffinity=change))
    assert protocol.cpu_affinity()["actual_mask"] == hex((1 << 2) | (1 << 5))
    assert protocol.cpu_affinity(1 << 5)["actual_mask"] == hex(1 << 5)
    with pytest.raises(ValueError, match="nonzero subset"):
        protocol.cpu_affinity(1 << 1)


def test_unsupported_affinity_is_explicit_and_observation_is_optional(monkeypatch):
    from types import SimpleNamespace
    monkeypatch.setattr(protocol, "os", SimpleNamespace(name="posix"))
    assert protocol.cpu_affinity() == dict(requested_mask=None, actual_mask=None, system_mask=None)
    with pytest.raises(RuntimeError, match="unavailable"):
        protocol.cpu_affinity(1)


def test_benchmark_does_not_normalize_retired_constructor_options(tmp_path):
    from benchmark import native_backends
    assert native_backends.ROOT == Path(__file__).resolve().parents[1]
    assert protocol.WORKLOAD == native_backends.ROOT / "benchmark/native_backends.py"
    # The maintained stream recipes must already use the current JaxProp schema.
    config, _ = native_backends.configuration("sco2_recuperated", "coolprop", 100, False, tmp_path)
    assert all("exceptions" not in value for key, value in config["problem_formulation"]["fixed_parameters"].items()
               if key.endswith("fluid"))
