import numpy as np

from no_free_collapse.hafnian_bounds import (
    disjoint_pair_basis_blocks,
    disjoint_pair_hafnian_directional_derivative,
    disjoint_pair_quadratic,
    hafnian,
    normalized_cube_max,
    six_variable_critical_hafnian_quadratic_coefficient,
    six_variable_critical_path_second_order_bound,
    twice_offdiag,
)


def _pair_basis_matrix() -> np.ndarray:
    P = np.zeros((6, 6), dtype=np.float64)
    scale = 1.0 / np.sqrt(2.0)
    for p in range(3):
        i = 2 * p
        P[i, p] = scale
        P[i + 1, p] = scale
        P[i, 3 + p] = scale
        P[i + 1, 3 + p] = -scale
    return P


def _critical_direction(d: np.ndarray, X: np.ndarray) -> np.ndarray:
    P = _pair_basis_matrix()
    pair_coordinates = np.block(
        [[np.diag(d), X], [X.T, np.zeros((3, 3), dtype=np.float64)]]
    )
    return P @ pair_coordinates @ P.T


def _exact_schur_path(d: np.ndarray, X: np.ndarray, epsilon: float) -> np.ndarray:
    """PSD path realizing a critical tangent with the minimal nullspace curvature."""
    P = _pair_basis_matrix()
    A = np.eye(3) / 6.0 + epsilon * np.diag(d)
    cross = epsilon * X
    N = cross.T @ np.linalg.inv(A) @ cross
    pair_coordinates = np.block([[A, cross], [cross.T, N]])
    return P @ pair_coordinates @ P.T


def test_pair_basis_blocks_round_trip_critical_structure():
    d = np.asarray([0.04, -0.01, -0.03])
    X = np.arange(9, dtype=np.float64).reshape(3, 3) / 100.0
    H = _critical_direction(d, X)
    FF, FG, GG = disjoint_pair_basis_blocks(H)
    assert np.max(np.abs(FF - np.diag(d))) < 1e-12
    assert np.max(np.abs(FG - X)) < 1e-12
    assert np.max(np.abs(GG)) < 1e-12
    assert abs(disjoint_pair_hafnian_directional_derivative(H)) < 1e-14


def test_exact_critical_quadratic_hafnian_coefficient():
    rng = np.random.default_rng(271828)
    B_star = disjoint_pair_quadratic(6)
    phi_star = hafnian(twice_offdiag(B_star)).real

    for _ in range(10):
        d = rng.normal(size=3)
        d -= d.mean()
        X = rng.normal(size=(3, 3))
        H = _critical_direction(d, X)

        # Phi(B_star+epsilon H) is cubic in epsilon.  Evaluating at +/-1
        # isolates the epsilon^2 coefficient exactly.
        phi_plus = hafnian(twice_offdiag(B_star + H)).real
        phi_minus = hafnian(twice_offdiag(B_star - H)).real
        exact_coefficient = 0.5 * (phi_plus + phi_minus - 2.0 * phi_star)
        formula = six_variable_critical_hafnian_quadratic_coefficient(H)
        assert abs(exact_coefficient - formula) < 1e-10


def test_every_nonzero_critical_direction_has_strict_negative_second_order_bound():
    rng = np.random.default_rng(161803)
    for _ in range(20):
        d = rng.normal(size=3)
        d -= d.mean()
        X = rng.normal(size=(3, 3))
        H = _critical_direction(d, X)
        assert six_variable_critical_path_second_order_bound(H) < -1e-12


def test_minimal_schur_curvature_path_attains_second_order_bound_and_is_feasible():
    d = np.asarray([0.03, -0.01, -0.02])
    X = np.asarray(
        [
            [0.015, -0.010, 0.005],
            [0.007, -0.012, 0.009],
            [-0.004, 0.006, 0.011],
        ]
    )
    H = _critical_direction(d, X)
    expected = six_variable_critical_path_second_order_bound(H)
    assert expected < 0.0

    epsilon = 1e-4
    B_plus = _exact_schur_path(d, X, epsilon)
    B_minus = _exact_schur_path(d, X, -epsilon)
    B_star = disjoint_pair_quadratic(6)

    for B in (B_plus, B_minus):
        assert np.linalg.eigvalsh(B)[0] >= -1e-12
        assert normalized_cube_max(B) <= 1.0 + 1e-12

    phi_star = hafnian(twice_offdiag(B_star)).real
    phi_plus = hafnian(twice_offdiag(B_plus)).real
    phi_minus = hafnian(twice_offdiag(B_minus)).real
    numerical = (phi_plus + phi_minus - 2.0 * phi_star) / (2.0 * epsilon**2)
    assert abs(numerical - expected) < 1e-7
