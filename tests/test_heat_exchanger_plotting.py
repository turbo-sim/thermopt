"""Check that T-Q diagrams pair temperatures at the same heat duty."""

import unittest

from matplotlib.figure import Figure
import numpy as np

from thermopt.components.basic_components import compute_component_energy_flows
from thermopt.optimize_cycle import ThermodynamicCycleProblem


class TestHeatExchangerPlotting(unittest.TestCase):
    def test_nonuniform_nodes_align_when_creating_and_updating_diagram(self):
        problem = ThermodynamicCycleProblem.__new__(ThermodynamicCycleProblem)
        problem.graphics = {"pinch_point_lines": {}}
        figure = Figure()
        axes = figure.subplots()

        # Different asymmetric grids exercise both creation and artist updates.
        for coordinate in ([0.0, 0.15, 0.65, 1.0], [0.0, 0.35, 0.8, 1.0]):
            with self.subTest(coordinate=coordinate):
                x = np.asarray(coordinate)
                hot = {
                    "mass_flow": 4.0,
                    "state_in": {"h": 150.0},
                    "states": {
                        "h": 100.0 + 50.0 * x,
                        "T": 350.0 + 20.0 * x,
                        "identifier": "working_fluid",
                    },
                }
                cold = {
                    "mass_flow": 2.0,
                    "state_in": {"h": 200.0},
                    "states": {
                        "h": 200.0 + 100.0 * x,
                        "T": 300.0 + 10.0 * x,
                        "identifier": "cooling_fluid",
                    },
                }
                components = {
                    "cooler": {
                        "type": "heat_exchanger",
                        "hot_side": hot,
                        "cold_side": cold,
                    }
                }
                compute_component_energy_flows(components)
                original_hot_duty = hot["heat_flow"].copy()
                problem.cycle_data = {"components": components}

                problem._plot_pinch_point_diagram(axes, 0)

                artists = problem.graphics["pinch_point_lines"][0]["cooler"]
                np.testing.assert_allclose(artists["hot_line"].get_xdata(), 200.0 * x)
                np.testing.assert_allclose(artists["cold_line"].get_xdata(), 200.0 * x)
                np.testing.assert_array_equal(
                    artists["hot_line"].get_ydata(), hot["states"]["T"]
                )
                np.testing.assert_array_equal(
                    artists["cold_line"].get_ydata(), cold["states"]["T"]
                )
                np.testing.assert_array_equal(hot["heat_flow"], original_hot_duty)
                np.testing.assert_allclose(artists["hot_start"].get_xdata(), [0.0])
                np.testing.assert_allclose(artists["hot_end"].get_xdata(), [200.0])


if __name__ == "__main__":
    unittest.main()
