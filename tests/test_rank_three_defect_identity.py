import numpy as np

from no_free_collapse.projection_gradient import (
    rank_three_projection_gradient_defect,
    rank_three_projection_stability_terms,
    six_variable_gradient_energies,
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
