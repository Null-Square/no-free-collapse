"""PSD diagonal-completion certificates for fixed off-diagonal interactions."""

from __future__ import annotations

import itertools

import numpy as np


def _validate_offdiag_coefficients(C: np.ndarray) -> np.ndarray:
    C = np.asarray(C, dtype=np.float64)
    if C.ndim != 2 or C.shape[0] != C.shape[1]:
        raise ValueError("C must be square")
    if not np.allclose(C, C.T, atol=1e-12, rtol=1e-12):
        raise ValueError("C must be symmetric")
    if not np.allclose(np.diag(C), 0.0, atol=1e-12, rtol=0.0):
        raise ValueError("C must have zero diagonal")
    return C


def offdiag_cube_extrema(C: np.ndarray) -> tuple[float, float]:
    """Return the exact min and max of the off-diagonal quadratic on a cube."""
    C = _validate_offdiag_coefficients(C)
    n = C.shape[0]
    worst = np.inf
    best = -np.inf
    for x_tuple in itertools.product((-1.0, 1.0), repeat=n):
        x = np.asarray(x_tuple, dtype=np.float64)
        value = 0.5 * float(x @ C @ x)
        worst = min(worst, value)
        best = max(best, value)
    return float(worst), float(best)


def offdiag_cube_max(C: np.ndarray) -> float:
    """Return ``max_x sum_{i<j} C_ij x_i x_j`` exactly on a small cube."""
    return offdiag_cube_extrema(C)[1]


def six_variable_range_hafnian_bound(C: np.ndarray) -> float:
    """Sharp moment-method range bound for a six-variable edge quadratic.

    Let

        r_C(x) = sum_{i<j} C_ij x_i x_j,
        a = -min_x r_C(x),
        s =  max_x r_C(x).

    Then, without any PSD assumption,

        |haf(C)| <= a*s*(a+s)/48.

    The proof quotients the global-sign symmetry to five Boolean variables.
    The two parity classes have equal first and second moments, and the sharp
    third-moment interval for a zero-mean variable in [-a,s] gives the factor
    1/48.
    """
    C = _validate_offdiag_coefficients(C)
    if C.shape != (6, 6):
        raise ValueError("the range hafnian bound is specific to size 6")
    minimum, maximum = offdiag_cube_extrema(C)
    a = max(0.0, -minimum)
    s = max(0.0, maximum)
    return a * s * (a + s) / 48.0


def completion_matrix(C: np.ndarray, diagonal: np.ndarray) -> np.ndarray:
    """Return the symmetric matrix with off-diagonal coefficients ``C/2``.

    If ``B=completion_matrix(C, d)``, then on the Boolean cube

        x^T B x = sum(d) + sum_{i<j} C_ij x_i x_j.
    """
    C = _validate_offdiag_coefficients(C)
    diagonal = np.asarray(diagonal, dtype=np.float64)
    n = C.shape[0]
    if diagonal.shape != (n,):
        raise ValueError("diagonal must have shape (n,)")
    return np.diag(diagonal) + 0.5 * C


def completion_primal_trace(
    C: np.ndarray,
    diagonal: np.ndarray,
    *,
    atol: float = 1e-10,
) -> float:
    """Verify a PSD diagonal completion and return its trace.

    This is a certificate checker, not an SDP solver.  A feasible diagonal gives
    an upper bound on the minimum completion trace ``tau(C)``.
    """
    B = completion_matrix(C, diagonal)
    if float(np.linalg.eigvalsh(B)[0]) < -atol:
        raise ValueError("the supplied diagonal does not give a PSD completion")
    return float(np.trace(B))


def completion_dual_value(
    C: np.ndarray,
    Y: np.ndarray,
    *,
    atol: float = 1e-10,
) -> float:
    """Verify an elliptope dual certificate and return its objective value.

    The minimum-trace PSD completion

        tau(C) = min sum_i d_i  subject to diag(d) + C/2 >= 0

    has the dual

        tau(C) = max -sum_{i<j} C_ij Y_ij
                 subject to Y >= 0, diag(Y)=1.

    A feasible ``Y`` therefore gives a rigorous lower bound on ``tau(C)``.
    """
    C = _validate_offdiag_coefficients(C)
    Y = np.asarray(Y, dtype=np.float64)
    n = C.shape[0]
    if Y.shape != (n, n):
        raise ValueError("Y must have the same shape as C")
    if not np.allclose(Y, Y.T, atol=atol, rtol=atol):
        raise ValueError("Y must be symmetric")
    if not np.allclose(np.diag(Y), 1.0, atol=atol, rtol=0.0):
        raise ValueError("Y must have unit diagonal")
    if float(np.linalg.eigvalsh(Y)[0]) < -atol:
        raise ValueError("Y must be positive semidefinite")
    return -float(np.sum(np.triu(C * Y, k=1)))
