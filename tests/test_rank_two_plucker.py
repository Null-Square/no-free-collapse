import numpy as np

from no_free_collapse.rank_two_projection import (
    rank_two_plucker_weights,
    rank_two_projection_direct_defect,
    rank_two_projection_gradient_defect,
    rank_two_tensor_q2_formula,
)
from no_free_collapse.projection_gradient import six_variable_gradient_energies


def _random_rank_two_projection(seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    q, _ = np.linalg.qr(rng.normal(size=(6, 2)))
    return q @ q.T


def test_plucker_weights_have_exact_vertex_marginals_and_total_mass():
    for seed in range(50):
        p = _random_rank_two_projection(seed)
        m = rank_two_plucker_weights(p)
        assert np.all(m >= -1e-12)
        assert np.allclose(np.sum(m, axis=1), np.diag(p), atol=2e-12, rtol=2e-12)
        assert np.isclose(np.sum(np.triu(m, 1)), 1.0, atol=2e-12, rtol=0.0)


def test_tensor_parseval_q2_formula_matches_direct_hafnians():
    for seed in range(50, 100):
        p = _random_rank_two_projection(seed)
        _, q2 = six_variable_gradient_energies(p)
        assert np.isclose(rank_two_tensor_q2_formula(p), q2, atol=3e-12, rtol=3e-12)


def test_plucker_defect_formula_matches_direct_gradient_defect():
    for seed in range(100, 150):
        p = _random_rank_two_projection(seed)
        direct = rank_two_projection_direct_defect(p)
        reduced = rank_two_projection_gradient_defect(p)
        assert np.isclose(reduced, direct, atol=3e-12, rtol=3e-12)


def test_pairwise_plucker_interval_constraints():
    for seed in range(150, 175):
        p = _random_rank_two_projection(seed)
        d = np.diag(p)
        m = rank_two_plucker_weights(p)
        for i in range(6):
            for j in range(i + 1, 6):
                lower = max(0.0, d[i] + d[j] - 1.0)
                upper = d[i] * d[j]
                assert m[i, j] >= lower - 3e-12
                assert m[i, j] <= upper + 3e-12
