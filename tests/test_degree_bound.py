import numpy as np
import pytest

from no_free_collapse.born import born_score, polynomial_state, preparation_order, squared_norm
from no_free_collapse.constructions import parity_state_coefficients, projector_zero
from no_free_collapse.interactions import evaluate_on_hypercube, interaction_degree, max_abs_by_order, walsh_coefficients


@pytest.mark.parametrize("r,n", [(1, 6), (2, 7), (3, 8)])
def test_raw_born_score_has_degree_at_most_twice_preparation_order(r, n):
    rng = np.random.default_rng(20260830 + r)
    dim = 5
    masks = [m for m in range(1 << n) if m.bit_count() <= r]
    coefficients = {
        m: rng.normal(size=dim) + 1j * rng.normal(size=dim)
        for m in masks
        if rng.random() < 0.45
    }
    if not coefficients:
        coefficients[0] = np.ones(dim, dtype=np.complex128)

    a = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
    effect = a.conj().T @ a

    values = evaluate_on_hypercube(
        n,
        lambda x: born_score(polynomial_state(coefficients, x), effect),
    )
    coeffs = walsh_coefficients(values, n)

    assert preparation_order(coefficients) <= r
    assert interaction_degree(coeffs, atol=1e-8) <= 2 * r
    for order, magnitude in max_abs_by_order(coeffs).items():
        if order > 2 * r:
            assert magnitude < 1e-8


@pytest.mark.parametrize("k", range(1, 11))
def test_parity_construction_is_fixed_norm_and_reaches_tight_order(k):
    coefficients = parity_state_coefficients(k)
    effect = projector_zero()

    values = {}
    for x in __import__("itertools").product((-1, 1), repeat=k):
        state = polynomial_state(coefficients, x)
        assert squared_norm(state) == pytest.approx(1.0, abs=1e-12)
        values[x] = born_score(state, effect)
        assert values[x] == pytest.approx((1 + np.prod(x)) / 2, abs=1e-12)

    coeffs = walsh_coefficients(values, k)
    assert interaction_degree(coeffs, atol=1e-10) == k
    assert preparation_order(coefficients) == (k + 1) // 2
