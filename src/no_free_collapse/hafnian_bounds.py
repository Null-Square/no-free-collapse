"""Hafnian bounds induced by bounded quadratic forms on the Boolean cube."""

from __future__ import annotations

import math

import numpy as np


def hafnian(matrix: np.ndarray) -> complex:
    """Exact recursive hafnian for small even matrices."""
    matrix = np.asarray(matrix)
    n = matrix.shape[0]
    if matrix.shape != (n, n):
        raise ValueError("matrix must be square")
    if n % 2:
        return 0.0
    if n == 0:
        return 1.0
    total = 0.0 + 0.0j
    for j in range(1, n):
        keep = [k for k in range(n) if k not in (0, j)]
        total += matrix[0, j] * hafnian(matrix[np.ix_(keep, keep)])
    return total


def twice_offdiag(matrix: np.ndarray) -> np.ndarray:
    """Return 2 times the off-diagonal part of a square matrix."""
    matrix = np.asarray(matrix)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("matrix must be square")
    return 2.0 * (matrix - np.diag(np.diag(matrix)))


def monic_chebyshev_error(power: int) -> float:
    """Best uniform error approximating t**power by degree power-1 on [0,1].

    Equivalently this is the minimum sup norm of a monic degree-``power``
    polynomial on [0,1], attained by a shifted/scaled Chebyshev polynomial.
    """
    if power < 1:
        raise ValueError("power must be positive")
    return 2.0 ** (1 - 2 * power)


def bounded_quadratic_hafnian_bound(size: int) -> float:
    """Universal hafnian bound for a bounded PSD quadratic form.

    Let ``size=2L`` and B be real PSD with

        max_{x in {+/-1}^size} x^T B x <= 1.

    Then for C=2*offdiag(B),

        |haf(C)| <= 2^(1-size) / L!.
    """
    if size < 2 or size % 2:
        raise ValueError("size must be a positive even integer")
    L = size // 2
    return monic_chebyshev_error(L) / math.factorial(L)


def disjoint_pair_hafnian_value(size: int) -> float:
    """Value size^(-size/2) attained by equal disjoint rank-one pairs."""
    if size < 2 or size % 2:
        raise ValueError("size must be a positive even integer")
    return size ** (-size / 2)


def disjoint_pair_quadratic(size: int) -> np.ndarray:
    """PSD extremizer for size=4 and a lower-bound witness for larger sizes.

    The returned B satisfies max_x x^T B x = 1 and has one rank-one 2x2 block
    for each matched pair.  Its scaled off-diagonal hafnian equals
    ``size**(-size/2)``.
    """
    if size < 2 or size % 2:
        raise ValueError("size must be a positive even integer")
    L = size // 2
    block = 1.0 / (4.0 * L)
    B = np.zeros((size, size), dtype=np.float64)
    for p in range(L):
        i = 2 * p
        B[i, i] = block
        B[i + 1, i + 1] = block
        B[i, i + 1] = block
        B[i + 1, i] = block
    return B


def disjoint_pair_basis_blocks(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Represent a symmetric matrix in orthonormal pair-sum/difference coordinates.

    The returned blocks are ``(FF, FG, GG)``, where the columns of F are
    ``(e_{2j}+e_{2j+1})/sqrt(2)`` and the columns of G are
    ``(e_{2j}-e_{2j+1})/sqrt(2)``.
    """
    matrix = np.asarray(matrix, dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("matrix must be square")
    n = matrix.shape[0]
    if n < 2 or n % 2:
        raise ValueError("matrix must have positive even dimension")
    if not np.allclose(matrix, matrix.T, atol=1e-12, rtol=1e-12):
        raise ValueError("matrix must be symmetric")

    L = n // 2
    F = np.zeros((n, L), dtype=np.float64)
    G = np.zeros((n, L), dtype=np.float64)
    scale = 1.0 / math.sqrt(2.0)
    for p in range(L):
        i = 2 * p
        F[i, p] = scale
        F[i + 1, p] = scale
        G[i, p] = scale
        G[i + 1, p] = -scale
    return F.T @ matrix @ F, F.T @ matrix @ G, G.T @ matrix @ G


def disjoint_pair_tangent_invariants(H: np.ndarray) -> tuple[float, float, float]:
    """Return the three first-order quantities in the pair local-optimum proof.

    For pair vectors ``f_j=e_{2j}+e_{2j+1}`` and
    ``g_j=e_{2j}-e_{2j+1}``, this returns

    ``(sum_j f_j^T H f_j, sum_j g_j^T H g_j,
       sum_j H[2j,2j+1])``.

    They satisfy the exact identity

        4 * matched_edge_sum = active_cube_average - nullspace_trace.

    Along every differentiable feasible path through the disjoint-pair point,
    the first quantity is non-positive and the second is non-negative.
    """
    H = np.asarray(H, dtype=np.float64)
    if H.ndim != 2 or H.shape[0] != H.shape[1]:
        raise ValueError("H must be square")
    n = H.shape[0]
    if n < 2 or n % 2:
        raise ValueError("H must have positive even dimension")
    if not np.allclose(H, H.T, atol=1e-12, rtol=1e-12):
        raise ValueError("H must be symmetric")

    active_cube_average = 0.0
    nullspace_trace = 0.0
    matched_edge_sum = 0.0
    for p in range(n // 2):
        i = 2 * p
        active_cube_average += H[i, i] + H[i + 1, i + 1] + 2.0 * H[i, i + 1]
        nullspace_trace += H[i, i] + H[i + 1, i + 1] - 2.0 * H[i, i + 1]
        matched_edge_sum += H[i, i + 1]
    return active_cube_average, nullspace_trace, matched_edge_sum


def disjoint_pair_hafnian_directional_derivative(H: np.ndarray) -> float:
    """Directional derivative of the scaled-offdiagonal hafnian at B_star.

    For ``n=2L`` and

        B_star = (1/(2n)) * diag(J_2, ..., J_2),
        Phi(B) = haf(2*offdiag(B)),

    the derivative in a symmetric direction H is

        D Phi(B_star)[H] = 2 n^(1-L) sum_j H[2j,2j+1].

    Hence any feasible tangent direction has non-positive derivative by the
    active-cube and PSD-nullspace inequalities encoded in
    :func:`disjoint_pair_tangent_invariants`.
    """
    H = np.asarray(H, dtype=np.float64)
    n = H.shape[0] if H.ndim == 2 else 0
    active, nullspace, matched = disjoint_pair_tangent_invariants(H)
    del active, nullspace
    L = n // 2
    return 2.0 * (n ** (1 - L)) * matched


def six_variable_critical_hafnian_quadratic_coefficient(H: np.ndarray) -> float:
    """Exact epsilon^2 coefficient for a critical n=6 pair perturbation.

    A first-order critical tangent has, in pair-sum/difference coordinates,

        H = [[D, X], [X^T, 0]],

    where D is diagonal and trace(D)=0.  For

        Phi(B)=haf(2*offdiag(B)),
        B(epsilon)=B_star+epsilon*H,

    this returns the exact coefficient of ``epsilon**2`` in Phi:

        -||diag(D)||_2^2/12 - ||offdiag(X)||_F^2/3.
    """
    H = np.asarray(H, dtype=np.float64)
    if H.shape != (6, 6):
        raise ValueError("H must have shape (6, 6)")
    FF, FG, GG = disjoint_pair_basis_blocks(H)
    offdiag_FF = FF - np.diag(np.diag(FF))
    if np.max(np.abs(GG)) > 1e-9:
        raise ValueError("critical H must have zero GG block")
    if np.max(np.abs(offdiag_FF)) > 1e-9 or abs(float(np.trace(FF))) > 1e-9:
        raise ValueError("critical H must have diagonal trace-zero FF block")

    d = np.diag(FF)
    offdiag_X = FG - np.diag(np.diag(FG))
    return -float(d @ d) / 12.0 - float(np.sum(offdiag_X**2)) / 3.0


def six_variable_critical_path_second_order_bound(H: np.ndarray) -> float:
    """Upper bound on the full epsilon^2 coefficient of any feasible C^2 path.

    For a twice-differentiable feasible path

        B(epsilon)=B_star+epsilon*H+(epsilon^2/2)K+o(epsilon^2)

    with a first-order critical nonzero H, PSD curvature forces
    ``K_GG >= 12 X^T X`` while active cube constraints give
    ``trace(K_FF)<=0``.  Consequently the K contribution is at most
    ``-||X||_F^2/6``.  Adding the exact pure-H coefficient gives this bound,
    which is strictly negative for every nonzero critical H.
    """
    pure = six_variable_critical_hafnian_quadratic_coefficient(H)
    _, X, _ = disjoint_pair_basis_blocks(H)
    return pure - float(np.sum(X**2)) / 6.0


def full_parity_power_coefficient(B: np.ndarray) -> float:
    """Exhaustively compute the full-parity coefficient of q(x)^L.

    Here q(x)=x^T B x and B has even dimension 2L.  This routine is only for
    small CPU verification; it enumerates the Boolean cube.
    """
    B = np.asarray(B, dtype=np.float64)
    n = B.shape[0]
    if B.shape != (n, n) or n % 2:
        raise ValueError("B must be square with even dimension")
    L = n // 2
    total = 0.0
    count = 1 << n
    for mask in range(count):
        x = np.asarray([1.0 if mask & (1 << i) else -1.0 for i in range(n)])
        parity = float(np.prod(x))
        q = float(x @ B @ x)
        total += parity * (q**L)
    return total / count


def normalized_cube_max(B: np.ndarray) -> float:
    """Exact max of x^T B x on a small Boolean cube."""
    B = np.asarray(B, dtype=np.float64)
    n = B.shape[0]
    if B.shape != (n, n):
        raise ValueError("B must be square")
    best = -math.inf
    for mask in range(1 << n):
        x = np.asarray([1.0 if mask & (1 << i) else -1.0 for i in range(n)])
        best = max(best, float(x @ B @ x))
    return best
