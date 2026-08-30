import math

import numpy as np

from no_free_collapse.projection_gradient import (
    SIX_VARIABLE_EDGES,
    complementary_four_hafnians,
    rank_three_projection_involution,
    six_variable_edge_vector,
    six_variable_gradient_energies,
    six_variable_perfect_matching_operator,
)


def _random_projection(rng: np.random.Generator, rank: int) -> np.ndarray:
    Q, _ = np.linalg.qr(rng.normal(size=(6, rank)))
    return Q @ Q.T


def _elementary(values: np.ndarray, degree: int) -> float:
    import itertools

    total = 0.0
    for subset in itertools.combinations(range(len(values)), degree):
        product = 1.0
        for i in subset:
            product *= float(values[i])
        total += product
    return total


def test_perfect_matching_operator_exact_gradient_identity():
    rng = np.random.default_rng(7)
    A = rng.normal(size=(6, 6))
    A = 0.5 * (A + A.T)
    a = six_variable_edge_vector(A)
    h = complementary_four_hafnians(A)
    T = six_variable_perfect_matching_operator(A)
    assert T.shape == (15, 15)
    assert np.allclose(T, T.T, atol=1e-12, rtol=0.0)
    assert np.allclose(T @ a, 2.0 * h, atol=1e-12, rtol=1e-12)


def test_gradient_energies_match_direct_norms():
    rng = np.random.default_rng(11)
    A = rng.normal(size=(6, 6))
    A = 0.5 * (A + A.T)
    q1, q2 = six_variable_gradient_energies(A)
    a = six_variable_edge_vector(A)
    h = complementary_four_hafnians(A)
    assert math.isclose(q1, float(a @ a), rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(q2, float(h @ h), rel_tol=0.0, abs_tol=1e-12)


def test_rank_one_projection_gradient_formula_and_bound():
    rng = np.random.default_rng(19)
    for _ in range(100):
        u = rng.normal(size=6)
        u /= np.linalg.norm(u)
        P = np.outer(u, u)
        d = u * u
        q1, q2 = six_variable_gradient_energies(P)
        e2 = _elementary(d, 2)
        e4 = _elementary(d, 4)
        assert math.isclose(q1, e2, rel_tol=0.0, abs_tol=2e-12)
        assert math.isclose(q2, 9.0 * e4, rel_tol=0.0, abs_tol=2e-12)
        assert q2 <= 0.25 * q1 + 2e-12


def test_balanced_rank_one_projection_attains_one_quarter_constant():
    u = np.full(6, 1.0 / math.sqrt(6.0))
    P = np.outer(u, u)
    q1, q2 = six_variable_gradient_energies(P)
    assert math.isclose(q2, 0.25 * q1, rel_tol=0.0, abs_tol=1e-12)


def test_projection_complement_preserves_gradient_energies():
    rng = np.random.default_rng(23)
    for rank in (1, 2, 3, 4, 5):
        P = _random_projection(rng, rank)
        q1, q2 = six_variable_gradient_energies(P)
        q1c, q2c = six_variable_gradient_energies(np.eye(6) - P)
        assert math.isclose(q1, q1c, rel_tol=0.0, abs_tol=2e-12)
        assert math.isclose(q2, q2c, rel_tol=0.0, abs_tol=2e-12)


def test_rank_three_involution_reduction_is_exact():
    rng = np.random.default_rng(29)
    for _ in range(25):
        P = _random_projection(rng, 3)
        K = rank_three_projection_involution(P)
        assert np.allclose(K, K.T, atol=1e-12, rtol=0.0)
        assert np.allclose(K @ K, np.eye(6), atol=2e-12, rtol=2e-12)
        assert math.isclose(float(np.trace(K)), 0.0, rel_tol=0.0, abs_tol=2e-12)

        q1p, q2p = six_variable_gradient_energies(P)
        q1k, q2k = six_variable_gradient_energies(K)
        # Off-diagonal entries scale by 1/2, while complementary 4x4
        # hafnians are quadratic in those entries.
        assert math.isclose(q1p, q1k / 4.0, rel_tol=0.0, abs_tol=2e-12)
        assert math.isclose(q2p, q2k / 16.0, rel_tol=0.0, abs_tol=2e-12)


def test_edge_order_has_all_fifteen_edges_once():
    assert len(SIX_VARIABLE_EDGES) == 15
    assert len(set(SIX_VARIABLE_EDGES)) == 15
