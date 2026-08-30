"""Six-variable hafnian-gradient and perfect-matching operator utilities.

The functions here expose exact finite-dimensional identities used in the
six-variable PSD extremal problem.  They are deliberately small and
CPU-verifiable; no optimization package is required.
"""

from __future__ import annotations

import itertools

import numpy as np

from .hafnian_bounds import hafnian


SIX_VARIABLE_EDGES: tuple[tuple[int, int], ...] = tuple(itertools.combinations(range(6), 2))


def _validate_six_symmetric(matrix: np.ndarray) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=np.float64)
    if matrix.shape != (6, 6):
        raise ValueError("matrix must have shape (6, 6)")
    if not np.allclose(matrix, matrix.T, atol=1e-12, rtol=1e-12):
        raise ValueError("matrix must be symmetric")
    return matrix


def six_variable_edge_vector(matrix: np.ndarray) -> np.ndarray:
    """Return the 15 off-diagonal entries in lexicographic edge order."""
    matrix = _validate_six_symmetric(matrix)
    return np.asarray([matrix[i, j] for i, j in SIX_VARIABLE_EDGES], dtype=np.float64)


def complementary_four_hafnians(matrix: np.ndarray) -> np.ndarray:
    """Return the 15 complementary 4x4 hafnians.

    The entry indexed by edge ``(i,j)`` is

        haf(matrix with rows/columns i,j deleted).

    Diagonal entries of ``matrix`` do not affect these hafnians.
    """
    matrix = _validate_six_symmetric(matrix)
    values = []
    for i, j in SIX_VARIABLE_EDGES:
        keep = [k for k in range(6) if k not in (i, j)]
        values.append(float(np.real(hafnian(matrix[np.ix_(keep, keep)]))))
    return np.asarray(values, dtype=np.float64)


def six_variable_perfect_matching_operator(matrix: np.ndarray) -> np.ndarray:
    """Return the 15x15 weighted perfect-matching operator T(matrix).

    Rows and columns are indexed by the 15 edges of K_6.  If two edges e,f
    intersect, T[e,f]=0.  If they are disjoint, there is a unique third edge g
    completing them to a perfect matching of the six vertices, and

        T[e,f] = matrix[g].

    For ``a=six_variable_edge_vector(matrix)`` and
    ``h=complementary_four_hafnians(matrix)``, the exact identity is

        T(matrix) @ a = 2 h.
    """
    matrix = _validate_six_symmetric(matrix)
    T = np.zeros((15, 15), dtype=np.float64)
    universe = set(range(6))
    for row, edge_e in enumerate(SIX_VARIABLE_EDGES):
        set_e = set(edge_e)
        for col, edge_f in enumerate(SIX_VARIABLE_EDGES):
            if not set_e.isdisjoint(edge_f):
                continue
            remaining = tuple(sorted(universe - set_e - set(edge_f)))
            if len(remaining) == 2:
                T[row, col] = matrix[remaining]
    return T


def six_variable_gradient_energies(matrix: np.ndarray) -> tuple[float, float]:
    """Return ``(q1,q2)`` for the six-variable hafnian gradient problem.

    Here

        q1 = sum_{i<j} matrix[i,j]^2,
        q2 = sum_{i<j} haf(matrix_{without i,j})^2.

    The conjectural global contraction for PSD matrices is

        q2 <= lambda_max(matrix)^2 * q1 / 4.

    This routine only evaluates the two sides; it does not assert the
    conjecture.
    """
    a = six_variable_edge_vector(matrix)
    h = complementary_four_hafnians(matrix)
    return float(a @ a), float(h @ h)


def rank_three_projection_involution(projection: np.ndarray, *, atol: float = 1e-10) -> np.ndarray:
    """Convert a rank-three orthogonal projection P to K=2P-I.

    For a valid rank-three projection, K is a symmetric orthogonal involution
    with trace zero.  The helper verifies the projection conditions before
    returning K.
    """
    P = _validate_six_symmetric(projection)
    if not np.allclose(P @ P, P, atol=atol, rtol=atol):
        raise ValueError("projection must satisfy P^2=P")
    if not np.isclose(np.trace(P), 3.0, atol=atol, rtol=0.0):
        raise ValueError("projection must have rank/trace three")
    K = 2.0 * P - np.eye(6)
    if not np.allclose(K @ K, np.eye(6), atol=10 * atol, rtol=10 * atol):
        raise ValueError("internal involution check failed")
    return K


def rank_three_projection_stability_terms(
    projection: np.ndarray,
    *,
    atol: float = 1e-10,
) -> tuple[float, float, float, float]:
    """Return the exact stability variables ``(W,L,S2,S4)`` for rank three.

    Let ``P`` be a rank-three orthogonal projection and write

        y_i = P_ii - 1/2,
        x_ij = P_ij^2  (i != j).

    The returned quantities are

        W  = sum_i sum_{j<k, j,k != i} x_ij x_ik,
        L  = sum_{i<j} x_ij (y_i-y_j)^2,
        S2 = sum_i y_i^2,
        S4 = sum_i y_i^4.

    Here W is the nonnegative wedge-spreading energy and

        L = (1/2) ||[diag(y),P]||_F^2

    is the nonnegative diagonal-imbalance commutator energy.
    """
    P = _validate_six_symmetric(projection)
    if not np.allclose(P @ P, P, atol=atol, rtol=atol):
        raise ValueError("projection must satisfy P^2=P")
    if not np.isclose(np.trace(P), 3.0, atol=atol, rtol=0.0):
        raise ValueError("projection must have rank/trace three")

    y = np.diag(P) - 0.5
    x = P * P
    np.fill_diagonal(x, 0.0)

    W = 0.0
    for i in range(6):
        row = np.asarray([x[i, j] for j in range(6) if j != i])
        W += 0.5 * (float(np.sum(row)) ** 2 - float(row @ row))

    L = 0.0
    for i, j in SIX_VARIABLE_EDGES:
        L += float(x[i, j] * (y[i] - y[j]) ** 2)

    S2 = float(y @ y)
    S4 = float(np.sum(y**4))
    return float(W), float(L), S2, S4


def rank_three_projection_gradient_defect(
    projection: np.ndarray,
    *,
    atol: float = 1e-10,
) -> float:
    """Return the exact rank-three gradient defect ``q1/4-q2``.

    For every rank-three projection P the following identity holds:

        q1(P)/4 - q2(P)
          = W + L/2 - (S2 + S2^2 - 10 S4)/8,

    where ``(W,L,S2,S4)`` are returned by
    :func:`rank_three_projection_stability_terms`.

    Consequently the remaining rank-three gradient contraction

        q2(P) <= q1(P)/4

    is equivalent to the single stability inequality

        8 W + 4 L + 10 S4 >= S2 + S2^2.

    This function evaluates the proved identity; it does not assert the final
    stability inequality globally.
    """
    W, L, S2, S4 = rank_three_projection_stability_terms(projection, atol=atol)
    return W + 0.5 * L - 0.125 * (S2 + S2 * S2 - 10.0 * S4)
