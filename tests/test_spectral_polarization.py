import itertools

import numpy as np

from no_free_collapse.projection_gradient import six_variable_gradient_energies
from no_free_collapse.spectral_polarization import (
    nested_projection_polarized_defect,
    normalized_spectral_chain,
    spectral_chain_polarized_average,
    zeon_degree_four_product,
)


def _random_orthogonal(seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.normal(size=(6, 6)))
    return Q


def _projection(Q: np.ndarray, rank: int) -> np.ndarray:
    if rank == 0:
        return np.zeros((6, 6))
    if rank == 6:
        return np.eye(6)
    return Q[:, :rank] @ Q[:, :rank].T


def _random_contraction(seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.normal(size=(6, 6)))
    eigenvalues = rng.uniform(0.0, 1.0, size=6)
    return Q @ np.diag(eigenvalues) @ Q.T


def test_zeon_square_norm_is_four_q2():
    for seed in range(20):
        A = _random_contraction(seed)
        _, q2 = six_variable_gradient_energies(A)
        product = zeon_degree_four_product(A, A)
        assert np.isclose(product @ product, 4.0 * q2, atol=3e-11, rtol=3e-11)


def test_all_equal_polarization_recovers_projection_defect():
    for rank in range(7):
        for seed in range(5):
            Q = _random_orthogonal(100 * rank + seed)
            P = _projection(Q, rank)
            q1, q2 = six_variable_gradient_energies(P)
            mixed = nested_projection_polarized_defect([P, P, P, P])
            assert np.isclose(mixed, 0.25 * q1 - q2, atol=3e-11, rtol=3e-11)


def test_identity_argument_makes_mixed_q2_part_vanish():
    # Omega_I=0, so any polarized term containing I has no quartic q2 part.
    Q = _random_orthogonal(2026)
    Ps = [_projection(Q, rank) for rank in (1, 2, 4)]
    I = np.eye(6)
    value = nested_projection_polarized_defect([*Ps, I])
    q1_part = 0.0
    for i, j in itertools.combinations(range(3), 2):
        q1_part += sum(Ps[i][a, b] * Ps[j][a, b] for a, b in itertools.combinations(range(6), 2))
    q1_part /= 24.0
    assert np.isclose(value, q1_part, atol=2e-12, rtol=2e-12)
    assert value >= -2e-12


def test_normalized_spectral_chain_reconstructs_contraction():
    for seed in range(20):
        A = _random_contraction(seed)
        scale, delta, projections = normalized_spectral_chain(A)
        B = A / scale
        reconstructed = sum(delta[k] * projections[k] for k in range(6))
        assert np.isclose(np.sum(delta), 1.0, atol=2e-12)
        assert np.min(delta) >= -2e-12
        assert np.allclose(reconstructed, B, atol=3e-11, rtol=3e-11)


def test_spectral_chain_polarization_identity():
    for seed in range(12):
        direct, polarized = spectral_chain_polarized_average(_random_contraction(seed))
        assert np.isclose(direct, polarized, atol=2e-9, rtol=2e-9)


def test_nested_mixed_defect_random_diagnostic_is_nonnegative():
    # This deterministic random scan is evidence for the still-conjectural
    # nested mixed inequality, not a substitute for its proof.
    rng = np.random.default_rng(20260831)
    for seed in range(40):
        Q = _random_orthogonal(seed)
        ranks = sorted(rng.integers(1, 7, size=4).tolist())
        Ps = [_projection(Q, rank) for rank in ranks]
        value = nested_projection_polarized_defect(Ps)
        assert value >= -2e-11
