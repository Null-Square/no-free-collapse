import itertools
import math

import numpy as np

from no_free_collapse.gram import feature_vector, monomial_masks, optimal_absolute_linear_readout_for_fixed_Q
from no_free_collapse.interactions import monomial
from no_free_collapse.symmetric import mean_field_asymptotic_rate, mean_field_full_parity_capacity, mean_field_full_parity_log_capacity, mean_field_gram


def _direct_capacity(n: int, kappa: float) -> tuple[float, float, int]:
    Q = mean_field_gram(n, kappa)
    masks = monomial_masks(n, 1)
    xs = list(itertools.product((-1, 1), repeat=n))
    features = np.asarray([feature_vector(x, masks) for x in xs])
    weights = np.asarray([monomial(x, (1 << n) - 1) / (1 << n) for x in xs])
    value, _, _ = optimal_absolute_linear_readout_for_fixed_Q(Q, features, weights)
    denominators = np.einsum("bi,ij,bj->b", features, Q, features).real
    condition = denominators.max() / denominators.min()
    return value, condition, np.linalg.matrix_rank(Q, tol=1e-10)


def test_mean_field_closed_form_matches_exact_gram_optimizer():
    for n in (4, 6, 8):
        for kappa in (2.0, 5.0, 10.0):
            direct, condition, rank = _direct_capacity(n, kappa)
            closed = mean_field_full_parity_capacity(n, kappa)
            assert abs(direct - closed) < 1e-10
            assert abs(condition - kappa) < 1e-10
            assert rank == n


def test_known_exact_values():
    assert math.isclose(mean_field_full_parity_capacity(4, 2.0), 0.075)
    assert math.isclose(mean_field_full_parity_capacity(6, 2.0), 1.0 / 104.0)
    assert math.isclose(mean_field_full_parity_capacity(6, 10.0), 0.125)


def test_log_capacity_avoids_underflow_and_matches_direct_value():
    for n in (10, 20, 40):
        value = mean_field_full_parity_capacity(n, 5.0)
        assert math.isclose(math.log(value), mean_field_full_parity_log_capacity(n, 5.0))


def test_finite_size_exponent_converges_toward_rate():
    kappa = 5.0
    rate = mean_field_asymptotic_rate(kappa)
    e40 = -mean_field_full_parity_log_capacity(40, kappa) / 40
    e160 = -mean_field_full_parity_log_capacity(160, kappa) / 160
    assert abs(e160 - rate) < abs(e40 - rate)
