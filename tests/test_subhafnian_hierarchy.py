"""Tests for the all-orders averaged subhafnian hierarchy and Roos hybrid."""
from math import comb

import numpy as np
import pytest

from no_free_collapse.subhafnian_hierarchy import (
    all_principal_hafnians,
    averaged_step_bound,
    closed_form_hybrid_bound,
    edge_density_upper_bound,
    hybrid_crossover_order,
    iterated_hybrid_bound,
    normalized_step_coefficients,
    odd_double_factorial,
    roos_G,
    spectral_diameter,
    subhafnian_energy,
)


def matching_projection(n: int) -> np.ndarray:
    a = np.eye(n) / 2
    for i in range(0, n, 2):
        a[i, i + 1] = a[i + 1, i] = 0.5 * (-1 if i % 4 else 1)
    return a


def balanced_rank_one(n: int) -> np.ndarray:
    s = np.where(np.arange(n) % 3 == 0, -1.0, 1.0)
    return np.outer(s, s) / n


@pytest.mark.parametrize("r", [2, 3, 4])
@pytest.mark.parametrize("family", ["rank_one", "rank_n_minus_one", "matching"])
def test_critical_families_sharpen_every_order(r, family):
    n = 4 * r - 2
    if family == "rank_one":
        a = balanced_rank_one(n)
    elif family == "rank_n_minus_one":
        a = np.eye(n) - balanced_rank_one(n)
    else:
        a = matching_projection(n)
    h = all_principal_hafnians(a)
    left, right = averaged_step_bound(a, r, hafnians=h)
    assert np.isclose(left, right, rtol=2e-12, atol=2e-13)
    assert comb(n, 2 * r) == comb(n, 2 * r - 2)
    assert np.isclose(subhafnian_energy(a, r, hafnians=h),
                      subhafnian_energy(a, r - 1, hafnians=h) / 4,
                      rtol=2e-12, atol=2e-13)


@pytest.mark.parametrize("n,r", [(6, 2), (8, 2), (10, 2), (10, 3), (12, 3), (14, 3), (14, 4)])
def test_random_real_symmetric_averaged_recurrence(n, r):
    rng = np.random.default_rng(1000 + 17 * n + r)
    for _ in range(5):
        a = rng.normal(size=(n, n)); a = (a + a.T) / 2
        h = all_principal_hafnians(a)
        left, right = averaged_step_bound(a, r, hafnians=h)
        assert left <= right * (1 + 3e-12) + 1e-12


@pytest.mark.parametrize("n,r", [(10, 3), (12, 3), (14, 4)])
def test_random_psd_hybrid_dominates_both_pure_chains(n, r):
    rng = np.random.default_rng(222 + n + r)
    for _ in range(5):
        q = np.linalg.qr(rng.normal(size=(n, n)))[0]
        vals = rng.uniform(0, 1, size=n)
        a = (q * vals) @ q.T
        h = all_principal_hafnians(a)
        actual, hybrid = iterated_hybrid_bound(a, r, hafnians=h)
        g1 = roos_G(a, 1, hafnians=h)
        d = spectral_diameter(a)
        pure_roos = g1**r
        pure_spectral = g1 * np.prod([d*d/(4*(2*s-1)**2) for s in range(2, r+1)])
        assert actual <= hybrid * (1 + 4e-12) + 1e-14
        assert hybrid <= pure_roos + 1e-16
        assert hybrid <= pure_spectral + 1e-16


@pytest.mark.parametrize("r", [2, 3, 4, 5])
def test_matching_projection_is_strictly_in_spectral_regime(r):
    g1 = 1 / (4 * (4 * r - 3))
    spectral = 1 / (4 * (2 * r - 1) ** 2)
    assert spectral < g1
    assert np.isclose(g1 / spectral, (2*r-1)**2/(4*r-3))


@pytest.mark.parametrize("r", [2, 3, 4, 5])
def test_balanced_rank_one_sits_at_top_order_crossover(r):
    n = 4 * r - 2
    g1 = 1 / n**2
    spectral = 1 / (4 * (2 * r - 1) ** 2)
    assert g1 == spectral


@pytest.mark.parametrize("n", range(4, 15))
def test_edge_density_upper_bound_for_random_symmetric_matrices(n):
    rng = np.random.default_rng(800 + n)
    for _ in range(20):
        a = rng.normal(size=(n, n)); a = (a + a.T) / 2
        h = all_principal_hafnians(a)
        d = spectral_diameter(a)
        if d == 0:
            continue
        theta = 4 * roos_G(a, 1, hafnians=h) / d**2
        assert theta <= edge_density_upper_bound(n) + 2e-12


def test_normalizations_and_invalid_inputs():
    assert odd_double_factorial(0) == 1
    assert odd_double_factorial(1) == 1
    assert odd_double_factorial(2) == 3
    assert odd_double_factorial(4) == 105
    with pytest.raises(ValueError): odd_double_factorial(-1)
    with pytest.raises(ValueError): averaged_step_bound(np.eye(8), 3)
    with pytest.raises(ValueError): normalized_step_coefficients(np.eye(6), 1)
    with pytest.raises(ValueError): spectral_diameter(np.eye(6, dtype=complex))
    with pytest.raises(ValueError): all_principal_hafnians(np.zeros((3, 4)))


@pytest.mark.parametrize("r", [2, 3])
def test_middle_operator_hafnian_recursion_identity(r):
    from itertools import combinations
    from no_free_collapse.matching_dilation import matching_operator
    n = 4 * r - 2
    rng = np.random.default_rng(9000 + r)
    a = rng.integers(-3, 4, size=(n, n)).astype(float)
    a = np.triu(a, 1); a = a + a.T
    h = all_principal_hafnians(a)
    features = list(combinations(range(n), 2 * r - 2))
    v = np.array([h[sum(1 << i for i in e)] for e in features], dtype=float)
    complement_h = np.array([
        h[((1 << n) - 1) ^ sum(1 << i for i in e)] for e in features
    ], dtype=float)
    assert np.array_equal(matching_operator(a) @ v, r * complement_h)


def test_global_double_counting_identity_before_interlacing():
    from itertools import combinations
    n, r = 11, 3
    local_n = 4 * r - 2
    rng = np.random.default_rng(4545)
    a = rng.integers(-2, 3, size=(n, n)).astype(float)
    a = np.triu(a, 1); a = a + a.T
    h = all_principal_hafnians(a)
    global_hi = subhafnian_energy(a, r, hafnians=h)
    global_lo = subhafnian_energy(a, r - 1, hafnians=h)
    summed_hi = 0.0
    summed_lo = 0.0
    for u in combinations(range(n), local_n):
        local = a[np.ix_(u, u)]
        lh = all_principal_hafnians(local)
        summed_hi += subhafnian_energy(local, r, hafnians=lh)
        summed_lo += subhafnian_energy(local, r - 1, hafnians=lh)
    assert summed_hi == comb(n - 2*r, 2*r - 2) * global_hi
    assert summed_lo == comb(n - 2*r + 2, 2*r) * global_lo


def test_closed_form_hybrid_has_one_switch_and_matches_product():
    rng = np.random.default_rng(7711)
    for n, r in ((10, 3), (14, 4)):
        for _ in range(10):
            q = np.linalg.qr(rng.normal(size=(n, n)))[0]
            vals = rng.uniform(0, 1, size=n)
            a = (q * vals) @ q.T
            h = all_principal_hafnians(a)
            actual_product, product_bound = iterated_hybrid_bound(a, r, hafnians=h)
            actual_closed, closed_bound, k = closed_form_hybrid_bound(a, r, hafnians=h)
            assert np.isclose(actual_product, actual_closed, rtol=1e-13, atol=1e-16)
            assert np.isclose(product_bound, closed_bound, rtol=2e-13, atol=1e-16)
            assert k == hybrid_crossover_order(a, r, hafnians=h)
            coeffs = [normalized_step_coefficients(a, s, hafnians=h) for s in range(2, r + 1)]
            winners = [c["spectral_wins"] for c in coeffs]
            assert winners == sorted(winners)


def test_edge_density_equality_for_half_rank_constant_diagonal_projection():
    for n in (6, 8, 10, 12, 14):
        a = matching_projection(n)
        h = all_principal_hafnians(a)
        theta = 4 * roos_G(a, 1, hafnians=h) / spectral_diameter(a)**2
        assert np.isclose(theta, edge_density_upper_bound(n), rtol=1e-14, atol=1e-15)


def test_hybrid_upper_sequence_is_log_subadditive():
    n, r = 14, 4
    a = matching_projection(n)
    h = all_principal_hafnians(a)
    upper = {1: roos_G(a, 1, hafnians=h)}
    for s in range(2, r + 1):
        _, upper[s], _ = closed_form_hybrid_bound(a, s, hafnians=h)
    for total in range(2, r + 1):
        for left in range(1, total):
            assert upper[total] <= upper[left] * upper[total-left] + 1e-18
