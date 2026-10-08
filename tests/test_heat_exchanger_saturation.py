"""Fixed-size HX snapping, pressure mapping, and baseline isolation checks."""

import unittest
from unittest.mock import patch

import numpy as np
import equinox as eqx
import jaxprop as props

from thermopt.components.basic_components import heat_exchanger
from thermopt.utilities.optimization_utils import evaluate_constraints


class TestHeatExchangerSaturation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.water = props.Fluid("Water", "HEOS")
        cls.organic = props.Fluid("Cyclopentane", "HEOS")

    def args(self, drop=0, condensing=False):
        p = 1e5
        hl = float(self.organic.get_state(props.PQ_INPUTS, p, 0).h)
        hv = float(self.organic.get_state(props.PQ_INPUTS, p, 1).h)
        low, high = hl - .35*(hv-hl), hv + .35*(hv-hl)
        if condensing:
            return (self.organic, high, low, p, p*(1-drop),
                    self.water, 1e5, 2e5, 1e6, 1e6*(1-drop))
        return (self.water, 6e5, 2e5, 1e6, 1e6*(1-drop),
                self.organic, low, high, p, p*(1-drop))

    @staticmethod
    def duty(side):
        return ((np.asarray(side["states"].h) - float(side["state_in"].h))
                / float(side["state_out"].h - side["state_in"].h))

    def assert_endpoints_and_balance(self, base, out, counter):
        xh, xc = self.duty(out["hot_side"]), self.duty(out["cold_side"])
        np.testing.assert_allclose(1-xh if counter else xh, xc, rtol=0, atol=1e-9)
        self.assertTrue(np.all(np.diff(xc) > 0))
        self.assertEqual(set(out), set(base))
        self.assertEqual(len(out["temperature_difference"]), len(base["temperature_difference"]))
        for side in ("hot_side", "cold_side"):
            for prop in ("h", "p", "T"):
                for endpoint in ("state_in", "state_out"):
                    np.testing.assert_array_equal(out[side][endpoint][prop], base[side][endpoint][prop])
                np.testing.assert_array_equal(np.asarray(out[side]["states"][prop])[[0, -1]],
                                              np.asarray(base[side]["states"][prop])[[0, -1]])

    def test_snap_inside_nodes_and_map_opposite_pressure(self):
        for counter in (False, True):
            for condensing in (False, True):
                for drop in (0, .01):
                    with self.subTest(counter=counter, condensing=condensing, drop=drop):
                        args = self.args(drop, condensing)
                        base = heat_exchanger(*args, num_steps=12, counter_current=counter)
                        out = heat_exchanger(*args, num_steps=12, counter_current=counter,
                                             include_saturation_nodes=True)
                        phase = "hot_side" if condensing else "cold_side"
                        other = "cold_side" if condensing else "hot_side"
                        before, after = base[phase]["states"], out[phase]["states"]
                        wet = ((np.asarray(before.is_two_phase) == 1)
                               & (np.asarray(before.Q) > 0) & (np.asarray(before.Q) < 1))
                        expected = {}
                        for left, right in zip(range(11), range(1, 12)):
                            if wet[left] == wet[right]:
                                continue
                            inside, outside = (left, right) if wet[left] else (right, left)
                            if inside in (0, 11):
                                continue
                            bubble = self.organic.get_state(props.PQ_INPUTS, float(before.p[outside]), 0)
                            expected[inside] = 0 if before.h[outside] < bubble.h else 1
                        self.assertEqual(len(expected), 2)
                        changed = np.flatnonzero(np.abs(np.asarray(after.h) - np.asarray(before.h)) > 1e-5)
                        self.assertEqual(set(changed), set(expected))
                        for index, quality in expected.items():
                            self.assertAlmostEqual(float(after.Q[index]), quality, places=10)
                            self.assertEqual(float(after.is_two_phase[index]), 1)
                            saturated = self.organic.get_state(props.PQ_INPUTS, float(after.p[index]), quality)
                            self.assertAlmostEqual(float(after.h[index]), float(saturated.h), delta=1e-5)
                            self.assertAlmostEqual(float(after.T[index]), float(saturated.T), delta=1e-8)
                        untouched = [i for i in range(12) if i not in changed]
                        for side in ("hot_side", "cold_side"):
                            for prop in ("h", "p", "T"):
                                np.testing.assert_array_equal(np.asarray(out[side]["states"][prop])[untouched],
                                                              np.asarray(base[side]["states"][prop])[untouched])
                            fraction = self.duty(out[side])
                            path_p = (float(out[side]["state_in"].p) + fraction
                                      * float(out[side]["state_out"].p - out[side]["state_in"].p))
                            np.testing.assert_allclose(out[side]["states"].p, path_p,
                                                       rtol=1e-6 if side == phase else 1e-9, atol=1e-3)
                        if drop:
                            self.assertTrue(np.all(np.abs(np.asarray(out[other]["states"].p)[changed]
                                                          - np.asarray(base[other]["states"].p)[changed]) > 1e-3))
                        self.assert_endpoints_and_balance(base, out, counter)
                        # Sample the original continuous h,p path directly at each
                        # corrected duty, avoiding dense-grid interpolation errors.
                        side = out[phase]
                        fraction = self.duty(side)[changed]
                        path_p = (float(side["state_in"].p) + fraction
                                  * float(side["state_out"].p - side["state_in"].p))
                        reference = self.organic.get_state(props.HmassP_INPUTS,
                                                          np.asarray(after.h)[changed], path_p)
                        np.testing.assert_allclose(np.asarray(after.T)[changed], reference.T,
                                                   rtol=0, atol=2e-4)

    def test_eos_call_budget(self):
        args = self.args(.01)
        with patch.object(self.water, "get_state", wraps=self.water.get_state) as hot, \
             patch.object(self.organic, "get_state", wraps=self.organic.get_state) as cold:
            for count in (12, 24):
                hot.reset_mock()
                cold.reset_mock()
                heat_exchanger(*args, num_steps=count, include_saturation_nodes=True)
                # Two baseline calls plus at most four PQ and one matching HP
                # call per crossing; no iteration or all-node saturation scan.
                self.assertLessEqual(hot.call_count + cold.call_count, 12)

    def test_disabled_skips_patch_and_is_bitwise_uniform(self):
        args = self.args(.01)
        with patch("thermopt.components.basic_components._snap_saturation_nodes") as snap, \
             patch.object(self.water, "get_state", wraps=self.water.get_state) as hot, \
             patch.object(self.organic, "get_state", wraps=self.organic.get_state) as cold:
            default = heat_exchanger(*args, num_steps=12)
            disabled = heat_exchanger(*args, num_steps=12, include_saturation_nodes=False)
            snap.assert_not_called()
            self.assertEqual(hot.call_count, 2)
            self.assertEqual(cold.call_count, 2)
        for side, offset in (("hot_side", 0), ("cold_side", 5)):
            fluid, hi, ho, pi, po = args[offset:offset+5]
            expected = fluid.get_state(props.HmassP_INPUTS, np.linspace(hi, ho, 12), np.linspace(pi, po, 12))
            if side == "hot_side":
                expected = expected.flipped()
            for prop in ("h", "p", "T", "Q", "is_two_phase"):
                np.testing.assert_array_equal(default[side]["states"][prop], expected[prop])
                np.testing.assert_array_equal(disabled[side]["states"][prop], expected[prop])

    def test_endpoint_dew_crossing_ignores_false_liquid_quality(self):
        p = 5e4
        liquid = self.organic.get_state(props.PQ_INPUTS, p, 0)
        vapor = self.organic.get_state(props.PQ_INPUTS, p, 1)
        latent = float(vapor.h - liquid.h)
        args = (self.organic, float(vapor.h) + .025*latent, float(liquid.h) - .025*latent,
                p, p, self.water, 8e4, 1e5, 1e6, 1e6)
        get_state = self.organic.get_state

        def false_quality(*args, **kwargs):
            state = get_state(*args, **kwargs)
            quality = np.where(np.asarray(state.is_two_phase) == 1, np.asarray(state.Q), 0.)
            return eqx.tree_at(lambda value: value.Q, state, quality)

        for counter in (False, True):
            with self.subTest(counter=counter), patch.object(self.organic, "get_state", side_effect=false_quality):
                base = heat_exchanger(*args, num_steps=10, counter_current=counter)
                out = heat_exchanger(*args, num_steps=10, counter_current=counter, include_saturation_nodes=True)
                index = 8 if counter else 1
                self.assertEqual(float(base["hot_side"]["state_in"].Q), 0)
                self.assertEqual(float(base["hot_side"]["state_in"].is_two_phase), 0)
                self.assertEqual(float(out["hot_side"]["states"].Q[index]), 1)
                self.assert_endpoints_and_balance(base, out, counter)

    def test_saturated_endpoints_are_not_duplicated(self):
        p = 1e5
        hl = float(self.organic.get_state(props.PQ_INPUTS, p, 0).h)
        hv = float(self.organic.get_state(props.PQ_INPUTS, p, 1).h)
        args = (self.water, 6e5, 2e5, 1e6, 1e6, self.organic, hl, hv, p, p)
        base = heat_exchanger(*args, num_steps=8)
        out = heat_exchanger(*args, num_steps=8, include_saturation_nodes=True)
        for side in ("hot_side", "cold_side"):
            for prop in ("h", "p", "T"):
                np.testing.assert_array_equal(out[side]["states"][prop], base[side]["states"][prop])

    def test_condenser_pinch_matches_dense_reference(self):
        p = 5e4
        hl = float(self.organic.get_state(props.PQ_INPUTS, p, 0).h)
        hv = float(self.organic.get_state(props.PQ_INPUTS, p, 1).h)
        latent = hv - hl
        args = (self.organic, hv + .025*latent, hl - .005*latent, p, p,
                self.water, 8e4, 1.05e5, 1e6, 1e6)
        coarse = heat_exchanger(*args, num_steps=10)
        snapped = heat_exchanger(*args, num_steps=10, include_saturation_nodes=True)
        dense = heat_exchanger(*args, num_steps=1001)
        coarse_min = float(min(coarse["temperature_difference"]))
        snapped_min = float(min(snapped["temperature_difference"]))
        dense_min = float(min(dense["temperature_difference"]))
        self.assertGreater(coarse_min - dense_min, .4)
        self.assertAlmostEqual(snapped_min, dense_min, delta=.01)

    def test_no_crossing_and_supercritical(self):
        for p in (1e6, 1.1*self.water.abstract_state.p_critical()):
            args = (self.water, 3e5, 3e5, p, p, self.water, 1e5, 1e5, p, p)
            base = heat_exchanger(*args, num_steps=5)
            out = heat_exchanger(*args, num_steps=5, include_saturation_nodes=True)
            np.testing.assert_array_equal(out["temperature_difference"], base["temperature_difference"])

    def test_coarse_competing_crossings_preserve_order(self):
        for n in (2, 3, 4):
            for both in (False, True):
                with self.subTest(n=n, both=both):
                    args = self.args()
                    if both:
                        _, low, high, p, _ = args[5:]
                        args = (self.organic, high, low, p, p, self.organic, low, high, p, p)
                    base = heat_exchanger(*args, num_steps=n)
                    out = heat_exchanger(*args, num_steps=n, include_saturation_nodes=True)
                    self.assert_endpoints_and_balance(base, out, True)
                    self.assertNotIn("minimum_temperature_difference", out)

    def test_constraint_count(self):
        for enabled in (False, True):
            out = heat_exchanger(*self.args(.01), num_steps=12, include_saturation_nodes=enabled)
            _, residuals, _ = evaluate_constraints({"components": {"heater": out}}, [
                {"variable": "$components.heater.temperature_difference", "type": ">", "value": 5}])
            self.assertEqual(len(residuals), 12)


if __name__ == "__main__":
    unittest.main()
