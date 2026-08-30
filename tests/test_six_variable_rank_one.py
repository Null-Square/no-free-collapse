import numpy as np

from no_free_collapse.hafnian_bounds import hafnian, normalized_cube_max, twice_offdiag


def _candidate_sides(u: np.ndarray) -> tuple[float, float, float, float, float]:
    B = np.outer(u, u)
    C = twice_offdiag(B)
    M = normalized_cube_max(B)
    t = float(np.trace(B))
    s = M - t
    h = abs(float(hafnian(C).real))
    return 54.0 * h, t * s * M, h, t, s


def test_rank_one_hafnian_closed_form():
    rng = np.random.default_rng(20260830)
    for _ in range(20):
        u = rng.normal(size=6)
        B = np.outer(u, u)
        actual = abs(float(hafnian(twice_offdiag(B)).real))
        expected = 120.0 * float(np.prod(np.abs(u)))
        assert abs(actual - expected) <= 1e-10 * max(1.0, expected)


def test_candidate_inequality_for_random_rank_one_psd_matrices():
    rng = np.random.default_rng(57721)
    for _ in range(100):
        u = rng.normal(size=6)
        lhs, rhs, _, _, _ = _candidate_sides(u)
        assert lhs <= rhs + 1e-10 * max(1.0, rhs)


def test_equal_absolute_coordinates_saturate_candidate_inequality():
    # Sum |u_i|=1, so max_x x^T B x = 1.
    u = np.asarray([1.0, -1.0, 1.0, 1.0, -1.0, -1.0]) / 6.0
    lhs, rhs, h, t, s = _candidate_sides(u)
    assert abs(normalized_cube_max(np.outer(u, u)) - 1.0) < 1e-12
    assert abs(t - 1.0 / 6.0) < 1e-12
    assert abs(s - 5.0 / 6.0) < 1e-12
    assert abs(h - 5.0 / 1944.0) < 1e-12
    assert abs(lhs - rhs) < 1e-12


def test_degenerate_rank_one_case_is_trivial_equality():
    u = np.asarray([2.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    lhs, rhs, h, _, s = _candidate_sides(u)
    assert h == 0.0
    assert s == 0.0
    assert lhs == rhs == 0.0
