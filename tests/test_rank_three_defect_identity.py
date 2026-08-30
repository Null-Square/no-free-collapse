import math

import numpy as np

from no_free_collapse.projection_gradient import (
    rank_three_projection_gradient_defect,
    rank_three_projection_stability_terms,
    six_variable_gradient_energies,
    six_variable_perfect_matching_operator,
)


def _random_rank_three_projection(seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.normal(size=(6, 3)))
    return Q @ Q.T


def _pair_projection() -> np.ndarray:
    P = np.zeros((6, 6), dtype=np.float64)
    for i in range(0, 6, 2):
        P[i : i + 2, i : i + 2] = 0.5
    return P


def _conference_involution() -> np.ndarray:
    signs = [
        1, 1, -1, -1, 1,
        1, 1, -1, -1,
        -1, 1, -1,
        -1, -1,
        -1,
    ]
    S = np.zeros((6, 6), dtype=np.float64)
    k = 0
    for i in range(6):
        for j in range(i + 1, 6):
            S[i, j] = S[j, i] = signs[k]
            k += 1
    assert np.allclose(S @ S, 5.0 * np.eye(6), atol=1e-12, rtol=0.0)
    return S / math.sqrt(5.0)


def test_rank_three_defect_identity_on_random_projections():
    for seed in range(50):
        P = _random_rank_three_projection(seed)
        q1, q2 = six_variable_gradient_energies(P)
        lhs = 0.25 * q1 - q2
        rhs = rank_three_projection_gradient_defect(P)
        assert np.isclose(lhs, rhs, atol=2e-12, rtol=2e-12)


def test_stability_terms_reconstruct_defect_exactly():
    for seed in range(20, 40):
        P = _random_rank_three_projection(seed)
        W, L, S2, S4 = rank_three_projection_stability_terms(P)
        rhs = W + 0.5 * L - 0.125 * (S2 + S2 * S2 - 10.0 * S4)
        assert np.isclose(
            rank_three_projection_gradient_defect(P), rhs, atol=2e-12, rtol=2e-12
        )


def test_equal_diagonal_pair_projection_is_zero_defect_equality_case():
    P = _pair_projection()
    W, L, S2, S4 = rank_three_projection_stability_terms(P)
    assert np.isclose(W, 0.0, atol=1e-12)
    assert np.isclose(L, 0.0, atol=1e-12)
    assert np.isclose(S2, 0.0, atol=1e-12)
    assert np.isclose(S4, 0.0, atol=1e-12)
    assert np.isclose(rank_three_projection_gradient_defect(P), 0.0, atol=1e-12)


def test_defect_target_is_equivalent_to_single_stability_inequality():
    for seed in range(60, 80):
        P = _random_rank_three_projection(seed)
        q1, q2 = six_variable_gradient_energies(P)
        W, L, S2, S4 = rank_three_projection_stability_terms(P)
        stability_margin = 8.0 * W + 4.0 * L + 10.0 * S4 - S2 - S2 * S2
        assert np.isclose(8.0 * (0.25 * q1 - q2), stability_margin, atol=2e-12, rtol=2e-12)


def test_full_matching_operator_norm_conjecture_is_false_exactly():
    K = _conference_involution()
    T = six_variable_perfect_matching_operator(K)
    op_norm = float(np.max(np.abs(np.linalg.eigvalsh(T))))
    assert math.isclose(op_norm, 3.0 / math.sqrt(5.0), rel_tol=0.0, abs_tol=2e-12)
    assert op_norm > 1.0

    # The false operator-norm strengthening does not contradict the actual
    # target, which only tests T on K's own edge vector.
    q1, q2 = six_variable_gradient_energies(K)
    assert q2 <= q1 + 1e-12
