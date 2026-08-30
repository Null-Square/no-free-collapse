import itertools
import math

import numpy as np

from no_free_collapse.gram import feature_vector, monomial_masks, optimal_absolute_linear_readout_for_fixed_Q
from no_free_collapse.interactions import monomial
from no_free_collapse.symmetric import (
    decreasing_branch_capacity,
    decreasing_branch_full_parity_capacity,
    increasing_branch_capacity,
    mean_field_full_parity_capacity,
    mean_field_gram,
    symmetric_branch_gram,
)


def _direct_capacity(Q: np.ndarray, n: int) -> tuple[float, float]:
    masks = monomial_masks(n, 1)
    xs = list(itertools.product((-1, 1), repeat=n))
    features = np.asarray([feature_vector(x, masks) for x in xs])
    weights = np.asarray([monomial(x, (1 << n) - 1) / (1 << n) for x in xs])
    value, _, _ = optimal_absolute_linear_readout_for_fixed_Q(Q, features, weights)
    denominators = np.einsum("bi,ij,bj->b", features, Q, features).real
    return value, float(denominators.max() / denominators.min())


def test_increasing_branch_formula_matches_exact_optimizer():
    for n in (4, 6, 8):
        for kappa in (2.0, 5.0, 10.0):
            for lam in (0.0, 0.25, 0.5, 0.75, 1.0):
                Q = symmetric_branch_gram(n, kappa, lam, increasing=True)
                direct, condition = _direct_capacity(Q, n)
                closed = increasing_branch_capacity(n, kappa, lam)
                assert abs(direct - closed) < 1e-10
                assert abs(condition - kappa) < 1e-10


def test_decreasing_branch_formula_matches_exact_optimizer():
    for n in (4, 6, 8):
        for kappa in (2.0, 5.0, 10.0):
            for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
                lam = frac / kappa
                Q = symmetric_branch_gram(n, kappa, lam, increasing=False)
                direct, condition = _direct_capacity(Q, n)
                closed = decreasing_branch_capacity(n, kappa, lam)
                assert abs(direct - closed) < 1e-10
                assert abs(condition - kappa) < 1e-10


def test_mean_field_is_best_allocation_on_each_branch():
    for n in (4, 6, 8):
        for kappa in (1.2, 2.0, 5.0, 10.0):
            positive = mean_field_full_parity_capacity(n, kappa)
            assert increasing_branch_capacity(n, kappa, 0.5) < positive
            negative = decreasing_branch_full_parity_capacity(n, kappa)
            assert decreasing_branch_capacity(n, kappa, 0.5 / kappa) < negative


def test_increasing_branch_strictly_beats_decreasing_for_n_at_least_four():
    for n in (4, 6, 8, 10):
        for kappa in (1.1, 2.0, 5.0, 10.0, 100.0):
            assert mean_field_full_parity_capacity(n, kappa) > decreasing_branch_full_parity_capacity(n, kappa)


def test_mean_field_capacity_increases_with_condition_budget():
    for n in (4, 6, 8, 10):
        vals = [mean_field_full_parity_capacity(n, k) for k in (1.1, 1.5, 2.0, 5.0, 10.0)]
        assert all(a < b for a, b in zip(vals, vals[1:]))


def test_mean_field_constructor_is_positive_branch_zero_constant_endpoint():
    for n in (4, 6):
        for kappa in (2.0, 5.0):
            assert np.allclose(mean_field_gram(n, kappa), symmetric_branch_gram(n, kappa, 0.0, increasing=True))
