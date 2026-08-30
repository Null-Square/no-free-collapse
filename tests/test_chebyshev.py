import math

from no_free_collapse.chebyshev import (
    asymptotic_chebyshev_rate,
    born_chebyshev_interaction_bound,
    chebyshev_residual_bound,
)
from no_free_collapse.conditioning import delta_from_condition_number


def test_chebyshev_residual_matches_closed_form():
    for kappa in (1.1, 2.0, 5.0, 10.0):
        c = (kappa + 1.0) / (kappa - 1.0)
        for m in range(1, 6):
            expected = 1.0 / math.cosh(m * math.acosh(c))
            assert abs(chebyshev_residual_bound(m, kappa) - expected) < 1e-12


def test_chebyshev_improves_geometric_bound():
    for kappa in (1.1, 2.0, 5.0, 10.0):
        delta = delta_from_condition_number(kappa)
        for m in range(1, 8):
            assert chebyshev_residual_bound(m, kappa) <= delta**m + 1e-12


def test_born_bound_has_helstrom_half_factor():
    assert math.isclose(born_chebyshev_interaction_bound(4, 1, 2.0), 1.0 / 6.0)
    assert math.isclose(born_chebyshev_interaction_bound(6, 1, 2.0), 1.0 / 34.0)


def test_fixed_norm_eliminates_high_order_interaction():
    assert born_chebyshev_interaction_bound(5, 1, 1.0) == 0.0


def test_asymptotic_rate_is_positive_and_decreases_with_conditioning():
    rates = [asymptotic_chebyshev_rate(kappa) for kappa in (1.1, 2.0, 5.0, 10.0)]
    assert all(rate > 0 for rate in rates)
    assert all(a > b for a, b in zip(rates, rates[1:]))
