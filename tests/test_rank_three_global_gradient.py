import numpy as np

from no_free_collapse.projection_gradient import six_variable_gradient_energies


def _random_rank_three_projection(seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    q, _ = np.linalg.qr(rng.normal(size=(6, 3)))
    return q @ q.T


def _pair_projection() -> np.ndarray:
    p = np.zeros((6, 6), dtype=np.float64)
    for i in range(0, 6, 2):
        p[i : i + 2, i : i + 2] = 0.5
    return p


def test_global_rank_three_gradient_contraction_random_projections():
    for seed in range(200):
        p = _random_rank_three_projection(seed)
        q1, q2 = six_variable_gradient_energies(p)
        assert q2 <= 0.25 * q1 + 2e-12


def test_pair_projection_attains_sharp_one_quarter_constant():
    p = _pair_projection()
    q1, q2 = six_variable_gradient_energies(p)
    assert np.isclose(q1, 0.75, atol=1e-12, rtol=0.0)
    assert np.isclose(q2, 3.0 / 16.0, atol=1e-12, rtol=0.0)
    assert np.isclose(q2, 0.25 * q1, atol=1e-12, rtol=0.0)


def test_coordinate_projection_is_degenerate_equality():
    p = np.diag([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    q1, q2 = six_variable_gradient_energies(p)
    assert np.isclose(q1, 0.0, atol=1e-12)
    assert np.isclose(q2, 0.0, atol=1e-12)
