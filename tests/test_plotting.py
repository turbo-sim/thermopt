"""Headless rendering and numerical artist checks for thermopt figures.

These assert the plotted physics/data and update behavior, rather than pixels,
so changes in installed fonts or Matplotlib antialiasing do not break the suite.
All cycle properties are obtained through thermopt's jaxprop-backed evaluator.
"""

import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
import pytest

import thermopt as th
from thermopt.utilities import graphics

from regression_support import load_case

pytestmark = pytest.mark.plotting


@pytest.fixture(autouse=True)
def isolated_graphics():
    """Do not leave open windows or global plotting settings for other tests."""
    with matplotlib.rc_context():
        # Portable test font; application defaults may name a Windows font.
        matplotlib.rcParams["font.serif"] = ["DejaVu Serif"]
        yield
    plt.close("all")


@pytest.fixture
def cycle_problem(tmp_path):
    config = load_case("sco2_recuperated")["problem_formulation"]
    config["plot_settings"]["fluid"] = {
        "plot_saturation_line": False,
        "plot_critical_point": False,
    }
    config["plot_settings"]["diagrams"] = [
        {"x_prop": "s", "y_prop": "T", "N": 12},
        {"x_prop": "h", "y_prop": "p", "y_scale": "log", "N": 12},
    ]
    config["plot_settings"]["pinch_point_diagram"] = True
    problem = th.ThermodynamicCycleProblem(config, out_dir=str(tmp_path))
    problem.fitness(problem.x0.copy())
    return problem


def assert_cycle_artists_match_data(problem):
    """The curves must depict the current evaluator outputs, not stale data."""
    for index, (prop_x, prop_y) in enumerate((("s", "T"), ("h", "p"))):
        process = problem.graphics["process_lines"][index]["compressor"]
        states = problem.cycle_data["components"]["compressor"]["states"]
        np.testing.assert_allclose(process.get_xdata(), states[prop_x])
        np.testing.assert_allclose(process.get_ydata(), states[prop_y])
        endpoints = problem.graphics["state_points"][index]["compressor"]
        np.testing.assert_allclose(endpoints.get_xdata(), np.asarray(states[prop_x])[[0, -1]])
        np.testing.assert_allclose(endpoints.get_ydata(), np.asarray(states[prop_y])[[0, -1]])

    exchangers = problem.graphics["pinch_point_lines"][2]
    assert set(exchangers) == {"heater", "recuperator", "cooler"}
    for name, artists in exchangers.items():
        component = problem.cycle_data["components"][name]
        hot, cold = artists["hot_line"], artists["cold_line"]
        # Both streams must use the same heat-duty coordinate, including when
        # their mass flows differ, and retain the evaluator's temperature gap.
        np.testing.assert_allclose(hot.get_xdata(), cold.get_xdata(), atol=1e-6)
        np.testing.assert_allclose(
            np.asarray(hot.get_ydata()) - np.asarray(cold.get_ydata()),
            component["temperature_difference"],
        )
        for line in (hot, cold):
            assert len(line.get_xdata()) == len(component["temperature_difference"])
            assert np.all(np.isfinite(line.get_xydata()))


def test_cycle_diagrams_render_properties_and_update_existing_artists(cycle_problem, tmp_path):
    problem = cycle_problem
    # Include the real jaxprop saturation and critical-point background. The
    # small phase-line sampling keeps this integration check inexpensive.
    problem.plot_settings["fluid"].update(
        plot_saturation_line=True, plot_critical_point=True
    )
    problem.plot_cycle()
    figure = problem.figure
    assert len(figure.axes) == 3
    assert "Entropy" in figure.axes[0].get_xlabel()
    assert "Temperature" in figure.axes[0].get_ylabel()
    assert "Enthalpy" in figure.axes[1].get_xlabel()
    assert "Pressure" in figure.axes[1].get_ylabel()
    assert figure.axes[1].get_yscale() == "log"
    assert "Heat flow rate" in figure.axes[2].get_xlabel()
    assert problem.fluid.sat_liq is not None
    assert problem.fluid.sat_vap is not None
    assert_cycle_artists_match_data(problem)

    line_counts = [len(axis.lines) for axis in figure.axes]
    artists = {
        (index, name): line
        for index, lines in problem.graphics["process_lines"].items()
        for name, line in lines.items()
    }
    old_source_temperature = np.asarray(
        problem.graphics["process_lines"][0]["heater_hot_side"].get_ydata()
    ).copy()
    point = problem.x0.copy()
    point[problem.variable_names.index("heat_source_exit_temperature")] += 2.0
    problem.fitness(point)
    problem.plot_cycle()

    assert problem.figure is figure
    assert [len(axis.lines) for axis in figure.axes] == line_counts
    assert all(problem.graphics["process_lines"][index][name] is artist
               for (index, name), artist in artists.items())
    assert not np.array_equal(
        problem.graphics["process_lines"][0]["heater_hot_side"].get_ydata(),
        old_source_temperature,
    )
    assert_cycle_artists_match_data(problem)

    output = tmp_path / "cycle.png"
    figure.savefig(output, dpi=80)
    pixels = plt.imread(output)
    assert pixels.shape[0] > 100 and pixels.shape[1] > 100
    assert float(np.std(pixels[..., :3])) > 0.01
    assert output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")


def test_cycle_plot_recreates_a_closed_figure(cycle_problem):
    problem = cycle_problem
    problem.plot_cycle()
    previous = problem.figure
    previous_line = problem.graphics["process_lines"][0]["compressor"]
    plt.close(previous)
    problem.plot_cycle()
    assert problem.figure is not previous
    assert problem.graphics["process_lines"][0]["compressor"] is not previous_line
    assert len(problem.figure.axes) == 3
    assert_cycle_artists_match_data(problem)


@pytest.mark.parametrize("prop_x", ["h", "s"])
def test_source_and_sink_temperature_curves_share_working_fluid_abscissa(cycle_problem, prop_x):
    problem = cycle_problem
    for component, side, other in (
        ("heater", "hot_side", "cold_side"),
        ("cooler", "cold_side", "hot_side"),
    ):
        x, y, _ = problem._get_process_data(f"{component}_{side}", prop_x, "T")
        exchanger = problem.cycle_data["components"][component]
        np.testing.assert_allclose(x, exchanger[other]["states"][prop_x])
        np.testing.assert_allclose(y, exchanger[side]["states"].T)


def test_external_stream_artists_hide_when_diagram_excludes_them(cycle_problem):
    problem = cycle_problem
    figure, axis = plt.subplots()
    visible = {"x_prop": "s", "y_prop": "T"}
    excluded = {"x_prop": "h", "y_prop": "p"}
    problem._plot_cycle_process("heater_hot_side", visible, axis, ax_index=0)
    line = problem.graphics["process_lines"][0]["heater_hot_side"]
    endpoints = problem.graphics["state_points"][0]["heater_hot_side"]
    assert line.get_visible() and endpoints.get_visible()
    problem._plot_cycle_process("heater_hot_side", excluded, axis, ax_index=0)
    assert not line.get_visible() and not endpoints.get_visible()
    problem._plot_cycle_process("heater_hot_side", visible, axis, ax_index=0)
    assert line.get_visible() and endpoints.get_visible()
    assert len(axis.lines) == 2


@pytest.mark.parametrize("mode,scale", [("multiply", 2.5), ("add", -3.0)])
@pytest.mark.parametrize("coordinate", ["x", "y"])
def test_graphics_scaling_applies_to_lines_polygons_and_contours(mode, scale, coordinate):
    figure, axis = plt.subplots()
    line, = axis.plot(np.array([1.0, 2.0]), np.array([3.0, 4.0]))
    polygon = Polygon([[1.0, 2.0], [2.0, 2.0], [2.0, 3.0]])
    axis.add_patch(polygon)
    x, y = np.meshgrid(np.linspace(0, 2, 5), np.linspace(0, 2, 5))
    contours = axis.contour(x, y, x + y, levels=[0.5, 1.5])
    index = 0 if coordinate == "x" else 1
    original_line = line.get_xydata().copy()
    original_polygon = polygon.get_xy().copy()
    original_paths = [path.vertices.copy() for path in contours.get_paths()]

    getattr(graphics, f"scale_graphics_{coordinate}")(figure, scale, mode=mode)

    def expected(values):
        changed = values.copy()
        changed[:, index] = values[:, index] * scale if mode == "multiply" else values[:, index] + scale
        return changed

    np.testing.assert_allclose(line.get_xydata(), expected(original_line))
    np.testing.assert_allclose(polygon.get_xy(), expected(original_polygon))
    for path, previous in zip(contours.get_paths(), original_paths, strict=True):
        np.testing.assert_allclose(path.vertices, expected(previous))


def test_figure_export_produces_raster_and_vector_documents(tmp_path):
    figure, axis = plt.subplots()
    axis.plot([0, 1, 2], [0, 2, 1], label="process")
    axis.set(xlabel="Heat duty [W]", ylabel="Temperature [K]")
    axis.legend()
    basename = tmp_path / "process"
    graphics.savefig_in_formats(figure, str(basename), formats=[".png", ".svg", ".pdf"], dpi=60)
    assert (tmp_path / "process.png").read_bytes().startswith(b"\x89PNG")
    assert "<svg" in (tmp_path / "process.svg").read_text(encoding="utf-8")
    assert (tmp_path / "process.pdf").read_bytes().startswith(b"%PDF")
    assert all(path.stat().st_size > 500 for path in tmp_path.iterdir())


def test_figure_export_rejects_unknown_formats_before_writing(tmp_path):
    figure, _ = plt.subplots()
    with pytest.raises(ValueError, match="Unsupported file formats"):
        graphics.savefig_in_formats(figure, str(tmp_path / "process"), formats=[".png", ".xyz"])
    assert list(tmp_path.iterdir()) == []


def test_plot_options_accept_a_custom_palette_and_restore_global_state():
    original_size = matplotlib.rcParams["font.size"]
    with matplotlib.rc_context():
        graphics.set_plot_options(fontsize=17, grid=False, margin=0.1,
                                  color_order=["red", "blue"], linewidth=2)
        figure, axis = plt.subplots()
        line, = axis.plot([0, 1], [0, 1])
        assert matplotlib.rcParams["font.size"] == 17
        assert matplotlib.rcParams["axes.grid"] is False
        assert axis.margins() == (0.1, 0.1)
        assert line.get_color() == "red"
        assert line.get_linewidth() == 2
    assert matplotlib.rcParams["font.size"] == original_size
