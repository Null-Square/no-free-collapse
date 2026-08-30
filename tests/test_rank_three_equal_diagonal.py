import itertools
import math

import numpy as np

from no_free_collapse.projection_gradient import six_variable_gradient_energies


def _hafnian4(matrix: np.ndarray) -> float:
    return float(
        matrix[0, 1] * matrix[2, 3]
        + matrix[0, 2] * matrix[1, 3]
        + matrix[0, 3] * matrix[1, 2]
    )


def _complementary_hafnian_energy(matrix: np.ndarray) -> float:
    total = 0.0
    for i, j in itertools.combinations(range(6), 2):
        keep = [k for k in range(6) if k not in (i, j)]
        value = _hafnian4(matrix[np.ix_(keep, keep)])
        total += value * value
    return total


def _random_orthogonal_3(rng: np.random.Generator) -> np.ndarray:
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    return Q


def test_zero_diagonal_involution_exact_fourth_power_identity():
    rng = np.random.default_rng(101)
    for _ in range(100):
        R = _random_orthogonal_3(rng)
        zero = np.zeros((3, 3), dtype=np.float64)
        K = np.block([[zero, R.T], [R, zero]])
        assert np.allclose(np.diag(K), 0.0, atol=1e-12, rtol=0.0)
        assert np.allclose(K @ K, np.eye(6), atol=2e-12, rtol=2e-12)

        q2 = _complementary_hafnian_energy(K)
        fourth = sum(K[i, j] ** 4 for i, j in itertools.combinations(range(6), 2))
        assert math.isclose(q2, fourth, rel_tol=0.0, abs_tol=5e-12)
        assert q2 <= 3.0 + 5e-12


def test_conference_involution_also_satisfies_identity():
    signs = [
        1, 1, -1, -1, 1,
        1, 1, -1, -1,
        -1, 1, -1,
        -1, -1,
        -1,
    ]
    S = np.zeros((6, 6), dtype=np.float64)
    cursor = 0
    for i in range(6):
        for j in range(i + 1, 6):
            S[i, j] = S[j, i] = signs[cursor]
            cursor += 1
    K = S / math.sqrt(5.0)
    assert np.allclose(K @ K, np.eye(6), atol=2e-12, rtol=2e-12)
    assert np.allclose(np.diag(K), 0.0, atol=1e-12, rtol=0.0)

    q2 = _complementary_hafnian_energy(K)
    fourth = sum(K[i, j] ** 4 for i, j in itertools.combinations(range(6), 2))
    assert math.isclose(q2, fourth, rel_tol=0.0, abs_tol=2e-12)
    assert math.isclose(q2, 3.0 / 5.0, rel_tol=0.0, abs_tol=2e-12)


def test_equal_diagonal_rank_three_projection_gradient_bound():
    rng = np.random.default_rng(131)
    for _ in range(100):
        R = _random_orthogonal_3(rng)
        zero = np.zeros((3, 3), dtype=np.float64)
        K = np.block([[zero, R.T], [R, zero]])
        P = 0.5 * (np.eye(6) + K)
        assert np.allclose(P @ P, P, atol=2e-12, rtol=2e-12)
        assert np.allclose(np.diag(P), 0.5, atol=1e-12, rtol=0.0)

        q1, q2 = six_variable_gradient_energies(P)
        assert math.isclose(q1, 0.75, rel_tol=0.0, abs_tol=3e-12)
        assert q2 <= 0.25 * q1 + 5e-12


def test_signed_perfect_matching_is_sharp_equality_case():
    K = np.zeros((6, 6), dtype=np.float64)
    signs = (1.0, -1.0, 1.0)
    for pair, sign in enumerate(signs):
        i = 2 * pair
        K[i, i + 1] = K[i + 1, i] = sign
    P = 0.5 * (np.eye(6) + K)

    q1, q2 = six_variable_gradient_energies(P)
    assert math.isclose(q1, 3.0 / 4.0, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(q2, 3.0 / 16.0, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(q2, 0.25 * q1, rel_tol=0.0, abs_tol=1e-12)
