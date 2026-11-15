"""Independent published-curve comparisons and headless reproduction checks."""

import json
from pathlib import Path

import jaxprop as props
import matplotlib.pyplot as plt
import numpy as np
import pytest

from development.turbomachinery import paper_figures as figures
from development.turbomachinery.extract_paper_data import axis_coordinates, svg_paths
from thermopt.components.turbomachinery_nondimensional import (
    compressor_reference_coefficients, compressor_mach_correction,
    compressor_relative_velocity_ratio, polytropic_head,
)


@pytest.fixture(scope="module")
def paper_data():
    """Evaluate a bounded plot grid; source comparisons always use all PDF vertices."""
    return figures.generate_data(points=15)


def test_reference_data_have_provenance_and_seven_complete_figures(paper_data):
    """Keep publication graphics attributable and distinct from raw validation data."""
    reference = paper_data["reference"]
    assert reference["source"]["doi"] == "10.1016/j.enconman.2025.120085"
    assert reference["source"]["license"] == "CC BY 4.0"
    assert len(reference["source"]["sha256"]) == 64
    assert set(reference["figures"]) == set(map(str, range(1, 8)))
    counts = {"1": [25, 25, 25], "2": [25, 25, 25], "4": [24, 24, 24],
              "5": [25, 25, 25], "6": [25], "7": [25]}
    for number, specification in reference["figures"].items():
        for source in specification["series"]:
            assert len(source["x"]) == len(source["y"]) >= 2
            assert np.all(np.diff(source["x"]) > 0)
            assert np.all(np.isfinite(source["y"]))
        if number in counts:
            assert [len(s["x"]) for s in specification["series"]
                    if s["kind"] == "published_meanline"] == counts[number]


@pytest.mark.parametrize("number", [1, 2, 4, 5, 6, 7])
def test_production_correlations_reproduce_independent_published_vectors(number, paper_data):
    """Compare the implementation with saved PDF curves, never with regenerated targets."""
    comparisons = [row for row in paper_data["comparisons"] if row["figure"] == number]
    assert len(comparisons) == 1
    assert comparisons[0]["failed_model_points"] == 0
    assert comparisons[0]["max_abs_error"] < 2e-4
    assert comparisons[0]["status"] == "agrees_with_published_curve"


def test_pressure_ratio_reconstruction_is_close_under_documented_assumptions(paper_data):
    """Use the declared approximate Figure 3 criterion, retaining measured errors."""
    comparisons = [row for row in paper_data["comparisons"] if row["figure"] == 3]
    assert len(comparisons) == 4
    for row in comparisons:
        assert row["failed_model_points"] == 0
        assert row["clipped_points_excluded"] >= 1
        assert row["max_abs_error"] < 5e-4
        assert row["status"] == "close_under_stated_assumptions"


def test_corrected_velocity_ratio_matches_independent_figure4_vertices(paper_data):
    """Protect the correction against Figure 4, independently of the erroneous Eq. 23."""
    coefficients = paper_data["figure4_published_quadratic"]
    np.testing.assert_allclose(coefficients, [-7.8, 3., .38], rtol=2e-5, atol=1e-5)
    source = next(series for series in paper_data["reference"]["figures"]["4"]["series"]
                  if series.get("compared"))
    calculated = compressor_relative_velocity_ratio(np.asarray(source["x"]))
    np.testing.assert_allclose(calculated, source["y"], rtol=0, atol=1e-5)
    row = next(row for row in paper_data["comparisons"] if row["figure"] == 4)
    assert row["status"] == "agrees_with_published_curve"
    assert row["max_abs_error"] < 1e-5


@pytest.mark.parametrize("pressure_ratio", [1.1, 3., 5., 10.])
def test_mach_only_reconstruction_satisfies_jaxprop_thermodynamic_efficiency(pressure_ratio):
    """Check the Figure 3 scalar solve against FluidPerfectGas endpoint states, Eq. 9."""
    fluid = props.FluidPerfectGas("Air", T_ref=298.)
    inlet = fluid.get_state(props.PT_INPUTS, 1e5, 298.)
    gamma = float(inlet.cp / inlet.cv)
    flow = .07
    efficiency = figures.mach_only_efficiency(flow, pressure_ratio, gamma)
    temperature_ratio = pressure_ratio**((gamma - 1) / (gamma * efficiency))
    outlet = fluid.get_state(props.PT_INPUTS, inlet.p * pressure_ratio, inlet.T * temperature_ratio)
    work = float(outlet.h - inlet.h)
    assert polytropic_head(inlet, outlet) / work == pytest.approx(efficiency, abs=1e-11)
    coefficient, head = compressor_reference_coefficients(flow)
    mach = np.sqrt(work / coefficient) / float(inlet.a)
    predicted = head / coefficient + compressor_mach_correction(flow, mach)
    assert predicted == pytest.approx(efficiency, abs=1e-11)


def test_svg_extraction_composes_transforms_and_keeps_glyphs_out(tmp_path):
    """Verify the extraction geometry independently of the actual paper or Poppler."""
    source = tmp_path / "drawing.svg"
    source.write_text('<svg xmlns="http://www.w3.org/2000/svg">'
                      '<defs><path d="M 1 1 L 2 2"/></defs>'
                      '<g transform="matrix(2,0,0,3,10,20)">'
                      '<path transform="matrix(1,0,0,1,4,5)" style="fill:none" d="M 0 0 L 1 2"/>'
                      '</g></svg>', encoding="utf-8")
    paths = svg_paths(source)
    assert len(paths) == 1
    np.testing.assert_allclose(paths[0]["endpoints"], [[18, 35], [20, 41]])
    specification = dict(frame=[10, 20, 30, 50], limits=[0, 1, 0, 3])
    x, y = axis_coordinates(paths[0]["endpoints"], specification)
    np.testing.assert_allclose(x, [.4, .5])
    np.testing.assert_allclose(y, [1.5, .9])


@pytest.mark.plotting
def test_all_seven_figures_and_comparison_reports_export(paper_data, tmp_path):
    """Exercise the complete headless plotting path and validate its output contracts."""
    output = figures.save_outputs(paper_data, tmp_path)
    for number in range(1, 8):
        png, svg = [output / f"figure_{number:02d}.{suffix}" for suffix in ("png", "svg")]
        assert png.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
        assert "<svg" in svg.read_text(encoding="utf-8")
        assert png.stat().st_size > 10000
    assert (output / "overview.png").stat().st_size > 10000
    report = json.loads((output / "comparison.json").read_text(encoding="utf-8"))
    assert len(report["comparisons"]) == 10
    assert not any(row["status"] == "discrepancy" for row in report["comparisons"])
    discussion = (output / "comparison.md").read_text(encoding="utf-8")
    assert "Equation 23" in discussion and "error" in discussion.lower()
    assert len((output / "calculated_curves.csv").read_text(encoding="utf-8").splitlines()) == 151
    assert not plt.get_fignums()
