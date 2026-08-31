import math

import numpy as np

from no_free_collapse.projection_gradient import six_variable_gradient_energies
from no_free_collapse.rank_two_dual import rank_two_quadratic_dual_certificate
from no_free_collapse.rank_two_projection import rank_two_projection_direct_defect


def _random_projection(rank: int, seed: int) -> np.ndarray:
    if rank == 0:
        return np.zeros((6, 6))
    if rank == 6:
        return np.eye(6)
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.normal(size=(6, rank)))
    return Q @ Q.T


def _low_edge_correction(x: float, y: float) -> float:
    s = x + y
    if not s < 0.5:
        return 0.0
    return 4.0 * x * y * (1.0 - s) * (1.0 - 2.0 * s)


def _seven_edge_symmetric_correction(S: float, Z: float) -> float:
    X = S - Z
    top = X / 3.0
    bottom = Z / 2.0
    return 6.0 * _low_edge_correction(top, bottom) + _low_edge_correction(bottom, bottom)


def test_low_edge_universal_rational_bound():
    s_star = (9.0 - math.sqrt(17.0)) / 16.0
    exact_max = s_star**2 * (1.0 - s_star) * (1.0 - 2.0 * s_star)
    assert exact_max < 13.0 / 512.0

    grid = np.linspace(0.0, 0.5, 401)
    sampled = 0.0
    for x in grid:
        for y in grid:
            if x + y < 0.5:
                sampled = max(sampled, _low_edge_correction(float(x), float(y)))
    assert sampled <= 13.0 / 512.0 + 1e-12


def test_middle_strip_constant_margins_are_positive():
    assert np.isclose(13.0 / 75.0 - 6.0 * 13.0 / 512.0, 403.0 / 19200.0)
    assert np.isclose(9.0 / 32.0 - 10.0 * 13.0 / 512.0, 7.0 / 256.0)
    assert np.isclose(13.0 / 75.0 - 9.0 / 64.0, 157.0 / 4800.0)


def test_seven_edge_jensen_reduction_has_exact_maximum():
    assert np.isclose(_seven_edge_symmetric_correction(1.0, 0.25), 9.0 / 64.0)

    # Enlarged domain used in the analytic proof.  The minimum feasible S is
    # max(1, Z+3/4), and the correction peaks at Z=1/4.
    zs = np.linspace(0.0, 0.5, 1001)
    values = []
    for Z in zs:
        S = max(1.0, float(Z) + 0.75)
        values.append(_seven_edge_symmetric_correction(S, float(Z)))
    assert max(values) <= 9.0 / 64.0 + 2e-7


def test_quadratic_dual_certificate_is_nonnegative_on_random_middle_diagonals():
    rng = np.random.default_rng(20260831)
    checked = 0
    for _ in range(20000):
        d = 2.0 * rng.dirichlet(np.ones(6))
        if np.max(d) > 1.0:
            continue
        d = np.sort(d)[::-1]
        if d[0] >= 0.5 and d[0] + d[1] <= 1.5:
            assert rank_two_quadratic_dual_certificate(d) >= -2e-12
            checked += 1
            if checked >= 500:
                break
    assert checked >= 250


def test_global_rank_two_gradient_contraction_randomly():
    for seed in range(250):
        P = _random_projection(2, seed)
        assert rank_two_projection_direct_defect(P) >= -8e-11


def test_rank_four_follows_by_projection_complementation():
    for seed in range(100):
        P = _random_projection(2, seed)
        Q = np.eye(6) - P
        q1p, q2p = six_variable_gradient_energies(P)
        q1q, q2q = six_variable_gradient_energies(Q)
        assert np.isclose(q1p, q1q, atol=2e-12, rtol=2e-12)
        assert np.isclose(q2p, q2q, atol=2e-12, rtol=2e-12)
        assert q2q <= 0.25 * q1q + 8e-11


def test_all_projection_ranks_obey_gradient_contraction_randomly():
    for rank in range(7):
        for seed in range(20):
            P = _random_projection(rank, 1000 * rank + seed)
            q1, q2 = six_variable_gradient_energies(P)
            assert q2 <= 0.25 * q1 + 1e-10
