import numpy as np

from no_free_collapse.rank_two_projection import (
    rank_two_harmonic_gradient_defect,
    rank_two_harmonic_invariants,
    rank_two_heavy_coordinate_decomposition,
    rank_two_heavy_coordinate_projection,
    rank_two_projection_direct_defect,
)


def _random_rank_two_projection(seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.normal(size=(6, 2)))
    return Q @ Q.T


def _random_orthonormal_pair(seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.normal(size=(5, 2)))
    return Q[:, 0], Q[:, 1]


def test_harmonic_defect_matches_direct_hafnian_evaluation():
    for seed in range(50):
        P = _random_rank_two_projection(seed)
        direct = rank_two_projection_direct_defect(P)
        harmonic = rank_two_harmonic_gradient_defect(P)
        assert np.isclose(harmonic, direct, atol=3e-12, rtol=3e-12)


def test_harmonic_invariants_are_latent_basis_independent_in_magnitude():
    P = _random_rank_two_projection(123)
    D2, D3, D4, A, B, S = rank_two_harmonic_invariants(P)
    d = np.diag(P)
    assert np.isclose(D2, np.sum(d**2), atol=1e-12)
    assert np.isclose(D3, np.sum(d**3), atol=1e-12)
    assert np.isclose(D4, np.sum(d**4), atol=1e-12)
    assert abs(A) <= D2 + 1e-12
    assert abs(B) <= D2 + 1e-12
    assert abs(S) <= 1.0 + 1e-12


def test_heavy_coordinate_decomposition_reconstructs_projection():
    for seed in range(20):
        p, q = _random_orthonormal_pair(seed)
        eps = 0.05 + 0.4 * (seed / 19.0)
        P = rank_two_heavy_coordinate_projection(eps, p, q)
        recovered_eps, recovered_p, recovered_q, index = rank_two_heavy_coordinate_decomposition(
            P, index=0
        )
        reconstructed = rank_two_heavy_coordinate_projection(
            recovered_eps, recovered_p, recovered_q
        )
        assert index == 0
        assert np.isclose(recovered_eps, eps, atol=2e-12, rtol=2e-12)
        assert np.allclose(reconstructed, P, atol=3e-11, rtol=3e-11)


def test_two_direction_harmonic_endpoint_has_unit_second_harmonic():
    # Coordinate plus a rank-one vector on the remaining five coordinates.
    q = np.ones(5, dtype=np.float64) / np.sqrt(5.0)
    p = np.zeros(5, dtype=np.float64)
    p[0] = 1.0
    # Orthogonalize q to p and renormalize, leaving a two-direction frame.
    q[0] = 0.0
    q /= np.linalg.norm(q)
    P = rank_two_heavy_coordinate_projection(0.3, p, q)
    *_, S = rank_two_harmonic_invariants(P)
    assert np.isclose(abs(S), 1.0, atol=2e-12, rtol=2e-12)
