import numpy as np

from no_free_collapse.hafnian_bounds import (
    disjoint_pair_hafnian_directional_derivative,
    disjoint_pair_quadratic,
    disjoint_pair_tangent_invariants,
    hafnian,
    normalized_cube_max,
    twice_offdiag,
)


def _pair_basis_tangent(size: int, active_values: np.ndarray, null_values: np.ndarray) -> np.ndarray:
    """Build H with prescribed f_j^T H f_j and g_j^T H g_j values."""
    H = np.zeros((size, size), dtype=np.float64)
    for p, (a, b) in enumerate(zip(active_values, null_values, strict=True)):
        i = 2 * p
        f = np.zeros(size)
        g = np.zeros(size)
        f[i] = f[i + 1] = 1.0
        g[i] = 1.0
        g[i + 1] = -1.0
        H += 0.25 * a * np.outer(f, f)
        H += 0.25 * b * np.outer(g, g)
    return H


def test_pair_tangent_identity_for_arbitrary_symmetric_directions():
    rng = np.random.default_rng(20260830)
    for size in (4, 6, 8):
        for _ in range(10):
            raw = rng.normal(size=(size, size))
            H = 0.5 * (raw + raw.T)
            active, nullspace, matched = disjoint_pair_tangent_invariants(H)
            assert abs(4.0 * matched - (active - nullspace)) < 1e-11


def test_pair_hafnian_derivative_matches_finite_difference():
    rng = np.random.default_rng(314159)
    epsilon = 1e-7
    for size in (4, 6, 8):
        raw = rng.normal(size=(size, size))
        H = 0.5 * (raw + raw.T)
        B = disjoint_pair_quadratic(size)
        base = hafnian(twice_offdiag(B)).real
        numerical = (hafnian(twice_offdiag(B + epsilon * H)).real - base) / epsilon
        analytic = disjoint_pair_hafnian_directional_derivative(H)
        assert abs(numerical - analytic) < 2e-5


def test_constructed_feasible_tangents_have_nonpositive_derivative():
    epsilon = 1e-4
    for size in (4, 6, 8):
        L = size // 2
        active_values = np.linspace(-0.3, 0.1, L)
        null_values = np.linspace(0.02, 0.08, L)
        H = _pair_basis_tangent(size, active_values, null_values)

        active, nullspace, matched = disjoint_pair_tangent_invariants(H)
        derivative = disjoint_pair_hafnian_directional_derivative(H)
        expected = 0.5 * (size ** (1 - L)) * (active - nullspace)

        assert active <= 0.0
        assert nullspace >= 0.0
        assert matched <= 0.0
        assert abs(derivative - expected) < 1e-14
        assert derivative <= 0.0

        B_epsilon = disjoint_pair_quadratic(size) + epsilon * H
        assert np.linalg.eigvalsh(B_epsilon)[0] >= -1e-12
        assert normalized_cube_max(B_epsilon) <= 1.0 + 1e-12
