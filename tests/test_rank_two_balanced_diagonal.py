import numpy as np

from no_free_collapse.projection_gradient import six_variable_gradient_energies
from no_free_collapse.rank_two_projection import rank_two_projection_gradient_defect


def _random_rank_two_projection(seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    q, _ = np.linalg.qr(rng.normal(size=(6, 2)))
    return q @ q.T


def _regular_six_direction_projection() -> np.ndarray:
    angles = np.arange(6, dtype=np.float64) * np.pi / 6.0
    u = np.column_stack((np.cos(angles), np.sin(angles))) / np.sqrt(3.0)
    assert np.allclose(u.T @ u, np.eye(2), atol=1e-12, rtol=0.0)
    return u @ u.T


def test_regular_multidirectional_frame_is_balanced_and_contracts():
    p = _regular_six_direction_projection()
    assert np.allclose(np.diag(p), np.full(6, 1.0 / 3.0), atol=1e-12, rtol=0.0)
    q1, q2 = six_variable_gradient_energies(p)
    assert q2 <= 0.25 * q1 + 1e-12
    assert rank_two_projection_gradient_defect(p) > 0.0


def test_random_balanced_rank_two_projections_respect_one_quarter_bound():
    tested = 0
    seed = 0
    while tested < 250 and seed < 20000:
        p = _random_rank_two_projection(seed)
        seed += 1
        if float(np.max(np.diag(p))) > 0.5 + 1e-12:
            continue
        q1, q2 = six_variable_gradient_energies(p)
        assert q2 <= 0.25 * q1 + 3e-12
        assert rank_two_projection_gradient_defect(p) >= -3e-12
        tested += 1
    assert tested == 250


def test_balanced_boundary_example_with_three_half_diagonals():
    # A Parseval frame with diagonal pattern (1/2,1/2,1/2,1/6,1/6,1/6).
    # The first three doubled angles form an equilateral triple; the last
    # three repeat the same directions with smaller weights.
    angles = np.asarray([0.0, np.pi / 3.0, 2.0 * np.pi / 3.0] * 2)
    d = np.asarray([0.5, 0.5, 0.5, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0])
    u = np.column_stack((np.sqrt(d) * np.cos(angles), np.sqrt(d) * np.sin(angles)))
    # Rotate/scaling pattern is chosen so the doubled-angle weighted sum is zero.
    assert np.allclose(u.T @ u, np.eye(2), atol=1e-12, rtol=0.0)
    p = u @ u.T
    q1, q2 = six_variable_gradient_energies(p)
    assert q2 <= 0.25 * q1 + 2e-12
