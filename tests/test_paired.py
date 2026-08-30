import itertools
import math

import numpy as np

from no_free_collapse.gram import feature_vector, monomial_masks, optimal_absolute_linear_readout_for_fixed_Q
from no_free_collapse.interactions import monomial
from no_free_collapse.paired import matched_pair_asymptotic_rate, matched_pair_capacity, matched_pair_gram
from no_free_collapse.symmetric import mean_field_asymptotic_rate, mean_field_full_parity_capacity


def _direct(Q: np.ndarray, n: int) -> tuple[float, float, int]:
    masks = monomial_masks(n, 1)
    xs = list(itertools.product((-1, 1), repeat=n))
    features = np.asarray([feature_vector(x, masks) for x in xs])
    weights = np.asarray([monomial(x, (1 << n) - 1) / (1 << n) for x in xs])
    value, _, _ = optimal_absolute_linear_readout_for_fixed_Q(Q, features, weights)
    q = np.einsum("bi,ij,bj->b", features, Q, features).real
    return value, float(q.max() / q.min()), int(np.linalg.matrix_rank(Q, tol=1e-10))


def test_closed_form_matches_exact_gram_optimizer():
    for n in (4, 6, 8):
        for kappa in (1.5, 2.0, 5.0, 10.0):
            value, condition, rank = _direct(matched_pair_gram(n, kappa), n)
            assert abs(value - matched_pair_capacity(n, kappa)) < 1e-10
            assert abs(condition - kappa) < 1e-10
            assert rank == n // 2 + 1


def test_exact_n4_counterexample():
    assert math.isclose(matched_pair_capacity(4, 2.0), 1.0 / 8.0)
    assert matched_pair_capacity(4, 2.0) > mean_field_full_parity_capacity(4, 2.0)


def test_pair_family_beats_mean_field_at_low_conditioning():
    for n in (4, 6, 8, 10):
        assert matched_pair_capacity(n, 2.0) > mean_field_full_parity_capacity(n, 2.0)


def test_asymptotic_rate_crossover_exists():
    assert matched_pair_asymptotic_rate(2.0) < mean_field_asymptotic_rate(2.0)
    assert mean_field_asymptotic_rate(10.0) < matched_pair_asymptotic_rate(10.0)
