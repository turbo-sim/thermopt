"""Numerical contracts for the smooth transitions used by cycle models.

The p-norm checks use positive inputs, its mathematically valid domain for
noninteger powers and negative-power minimum approximations. Known defects
are strict expected failures so a future fix prompts removal of the marker.
"""

import numpy as np
import pytest

from thermopt import math as tm


@pytest.mark.parametrize("method", ["boltzmann", "logsumexp", "p-norm"])
def test_smooth_extrema_are_symmetric_broadcast_elementwise_and_converge(method):
    first = np.array([[1.0], [3.0]])
    second = np.array([[2.0, 4.0, 5.0]])
    maximum = np.maximum(first, second)
    minimum = np.minimum(first, second)
    for function, exact in [(tm.smooth_maximum, maximum), (tm.smooth_minimum, minimum)]:
        diffuse = function(first, second, method=method, alpha=2)
        sharp = function(first, second, method=method, alpha=40)
        assert sharp.shape == (2, 3)
        assert np.all(np.isfinite(sharp))
        np.testing.assert_allclose(sharp, function(second, first, method=method, alpha=40))
        assert np.all(np.abs(sharp - exact) <= np.abs(diffuse - exact) + 1e-12)
        np.testing.assert_allclose(sharp, exact, rtol=0.02, atol=0.02)


@pytest.mark.parametrize("method", ["boltzmann", "logsumexp"])
def test_exponential_smooth_extrema_handle_large_signed_offsets(method):
    first = np.array([-1e9, 0.0, 1e9])
    second = first + 2
    # Direct unshifted exponentials would overflow for these states.
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        for function in [tm.smooth_maximum, tm.smooth_minimum]:
            result = function(first, second, method=method, alpha=100)
            assert np.all(np.isfinite(result))
            np.testing.assert_allclose(
                result - first,
                np.full(first.shape, float(function(0.0, 2.0, method=method, alpha=100))),
                atol=1e-6,
            )


@pytest.mark.parametrize("method", ["boltzmann", "logsumexp", "p-norm"])
def test_smooth_extrema_equal_inputs_have_known_bias(method):
    value, alpha = 4.0, 8.0
    if method == "boltzmann":
        expected_max = expected_min = value
    elif method == "logsumexp":
        expected_max = value + np.log(2) / alpha
        expected_min = value - np.log(2) / alpha
    else:
        expected_max = value * 2 ** (1 / alpha)
        expected_min = value * 2 ** (-1 / alpha)
    assert float(tm.smooth_maximum(value, value, method=method, alpha=alpha)) == pytest.approx(expected_max)
    assert float(tm.smooth_minimum(value, value, method=method, alpha=alpha)) == pytest.approx(expected_min)


@pytest.mark.parametrize("function", [tm.smooth_maximum, tm.smooth_minimum])
def test_logsumexp_matches_independent_stable_numpy_reference(function):
    first = np.array([[-10.0], [2.0]])
    second = np.array([[-9.0, 0.0, 5.0]])
    alpha = 3 if function is tm.smooth_maximum else -3
    expected = np.logaddexp(alpha * first, alpha * second) / alpha
    np.testing.assert_allclose(
        function(first, second, method="logsumexp", alpha=3), expected, atol=1e-12
    )


def test_exponential_extrema_have_expected_upper_lower_bounds():
    first = np.array([-4.0, 1.0, 5.0])
    second = np.array([-2.0, 1.0, 3.0])
    lower, upper = np.minimum(first, second), np.maximum(first, second)
    alpha = 5
    boltzmann_min = tm.smooth_minimum(first, second, alpha=alpha)
    boltzmann_max = tm.smooth_maximum(first, second, alpha=alpha)
    assert np.all(lower - 1e-12 <= boltzmann_min)
    assert np.all(boltzmann_min <= boltzmann_max + 1e-12)
    assert np.all(boltzmann_max <= upper + 1e-12)
    log_min = tm.smooth_minimum(first, second, method="logsumexp", alpha=alpha)
    log_max = tm.smooth_maximum(first, second, method="logsumexp", alpha=alpha)
    assert np.all(log_min <= lower + 1e-12)
    assert np.all(log_min >= lower - np.log(2) / alpha - 1e-12)
    assert np.all(log_max >= upper - 1e-12)
    assert np.all(log_max <= upper + np.log(2) / alpha + 1e-12)


def test_pnorm_has_known_two_norm_values_and_minimum_bound():
    assert float(tm.smooth_maximum(3.0, 4.0, method="p-norm", alpha=2)) == pytest.approx(5)
    assert float(tm.smooth_minimum(3.0, 4.0, method="p-norm", alpha=2)) == pytest.approx(2.4)
    first = np.array([1.0, 2.0, 3.0])
    second = np.array([4.0, 2.0, 1.0])
    result = tm.smooth_minimum(first, second, method="p-norm", alpha=4)
    exact = np.minimum(first, second)
    assert np.all(result <= exact)
    assert np.all(result >= 2 ** (-1 / 4) * exact)


@pytest.mark.parametrize("function", [tm.smooth_maximum, tm.smooth_minimum])
def test_unsupported_extremum_method_has_useful_error(function):
    with pytest.raises(ValueError, match="Unsupported method.*unknown"):
        function(1.0, 2.0, method="unknown")


@pytest.mark.parametrize("method", ["quadratic", "hyperbolic", "logarithmic"])
def test_smooth_absolute_is_even_nonnegative_and_has_zero_slope_at_origin(method):
    # 'logarithmic' is the currently accepted log-cosh selector.
    values = np.array([-2.0, -0.25, 0.0, 0.25, 2.0])
    result = tm.smooth_abs(values, method=method, epsilon=0.1)
    assert result.shape == values.shape
    assert np.all(result >= 0)
    np.testing.assert_allclose(result, result[::-1])
    step = 1e-5
    slope = (
        tm.smooth_abs(step, method=method, epsilon=0.1)
        - tm.smooth_abs(-step, method=method, epsilon=0.1)
    ) / (2 * step)
    assert float(slope) == pytest.approx(0, abs=1e-12)


def test_quadratic_absolute_has_known_values_and_uniform_error_bound():
    epsilon = 0.04
    values = np.array([-0.3, 0.0, 0.3])
    result = tm.smooth_abs(values, epsilon=epsilon)
    np.testing.assert_allclose(result, [np.sqrt(0.13), 0.2, np.sqrt(0.13)])
    assert np.all(result >= np.abs(values))
    assert np.all(result - np.abs(values) <= np.sqrt(epsilon) + 1e-12)


def test_hyperbolic_absolute_approaches_magnitude_from_below():
    values = np.array([-1.0, -0.1, 0.0, 0.1, 1.0])
    broad = tm.smooth_abs(values, method="hyperbolic", epsilon=0.1)
    narrow = tm.smooth_abs(values, method="hyperbolic", epsilon=0.001)
    assert np.all(broad <= np.abs(values))
    assert np.all(np.abs(narrow - np.abs(values)) <= np.abs(broad - np.abs(values)))
    np.testing.assert_allclose(narrow, np.abs(values), atol=1e-12)


def test_logarithmic_absolute_has_known_value_at_one_smoothing_width():
    epsilon = 0.2
    expected = epsilon * np.log((np.e + 1 / np.e) / 2)
    assert float(tm.smooth_abs(epsilon, method="logarithmic", epsilon=epsilon)) == pytest.approx(expected)


def test_unsupported_absolute_method_is_rejected():
    with pytest.raises(ValueError, match="Unsupported method.*unknown"):
        tm.smooth_abs(1.0, method="unknown")


@pytest.mark.xfail(strict=True, raises=ValueError, reason="Existing API defect: documented 'logcosh' selector is rejected; code accepts 'logarithmic'.")
def test_documented_logcosh_absolute_selector_is_accepted():
    assert float(tm.smooth_abs(0.2, method="logcosh", epsilon=0.2)) == pytest.approx(0.2 * np.log(np.cosh(1)))


@pytest.mark.xfail(strict=True, raises=FloatingPointError, reason="Existing stability defect: direct cosh(x/epsilon) overflows at ordinary inputs and default epsilon.")
def test_logarithmic_absolute_remains_finite_at_default_smoothing_scale():
    with np.errstate(over="raise", invalid="raise"):
        result = tm.smooth_abs(np.array([-1.0, 0.0, 1.0]), method="logarithmic")
    assert np.all(np.isfinite(result))
    np.testing.assert_allclose(result, [1.0, 0.0, 1.0], atol=1e-5)


@pytest.mark.parametrize("function", [tm.sigmoid_smoothstep, tm.sigmoid_smootherstep])
def test_polynomial_sigmoids_clip_outside_interval_and_have_flat_endpoints(function):
    values = np.linspace(-1, 2, 101)
    result = function(values)
    assert np.all((result >= 0) & (result <= 1))
    assert np.all(np.diff(result) >= -1e-12)
    np.testing.assert_allclose(function([-1, 0, 0.5, 1, 2]), [0, 0, 0.5, 1, 1])
    step = 1e-6
    assert float((function(step) - function(0)) / step) == pytest.approx(0, abs=4e-6)
    assert float((function(1) - function(1 - step)) / step) == pytest.approx(0, abs=4e-6)


def test_hyperbolic_sigmoid_is_centered_symmetric_and_stable_far_from_center():
    center, width = 2.0, 0.25
    offsets = np.array([-100.0, -0.5, 0.0, 0.5, 100.0])
    result = tm.sigmoid_hyperbolic(center + offsets, x0=center, alpha=width)
    assert np.all(np.isfinite(result))
    assert np.all(np.diff(result) > 0)
    assert result[2] == pytest.approx(0.5)
    np.testing.assert_allclose(result + result[::-1], 1)
    np.testing.assert_allclose(result[[0, -1]], [0, 1])


def test_rational_sigmoid_has_clipped_endpoints_and_symmetric_transition():
    values = np.array([-1.0, 0.0, 0.25, 0.5, 0.75, 1.0, 2.0])
    result = tm.sigmoid_rational(values, n=2, m=2)
    np.testing.assert_allclose(result, [0, 0, 0.1, 0.5, 0.9, 1, 1])
    np.testing.assert_allclose(result + result[::-1], 1)
