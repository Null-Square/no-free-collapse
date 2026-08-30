"""Exact rank-two projection identities for the six-variable gradient problem."""

from __future__ import annotations

import numpy as np

from .projection_gradient import SIX_VARIABLE_EDGES, six_variable_gradient_energies


def _validate_rank_two_projection(projection: np.ndarray, *, atol: float = 1e-10) -> np.ndarray:
    P = np.asarray(projection, dtype=np.float64)
    if P.shape != (6, 6):
        raise ValueError("projection must have shape (6, 6)")
    if not np.allclose(P, P.T, atol=atol, rtol=atol):
        raise ValueError("projection must be symmetric")
    if not np.allclose(P @ P, P, atol=atol, rtol=atol):
        raise ValueError("projection must satisfy P^2=P")
    if not np.isclose(np.trace(P), 2.0, atol=atol, rtol=0.0):
        raise ValueError("projection must have rank/trace two")
    return P


def rank_two_plucker_weights(
    projection: np.ndarray,
    *,
    atol: float = 1e-10,
) -> np.ndarray:
    """Return squared Pluecker coordinates of a rank-two projection.

    If ``P=UU^T`` with ``U`` a ``6 x 2`` Parseval frame and row vectors
    ``u_i``, then

        m_ij = det(u_i,u_j)^2
             = P_ii P_jj - P_ij^2.

    The returned symmetric matrix has zero diagonal and satisfies the exact
    marginal identities

        sum_{j != i} m_ij = P_ii,
        sum_{i<j} m_ij = 1.

    Thus the squared Pluecker coordinates form a probability distribution on
    the 15 edges of K_6 whose vertex marginals are the diagonal of P.
    """
    P = _validate_rank_two_projection(projection, atol=atol)
    d = np.diag(P)
    m = np.outer(d, d) - P * P
    np.fill_diagonal(m, 0.0)
    if float(np.min(m)) < -20.0 * atol:
        raise ValueError("rank-two principal minors must be nonnegative")
    m = np.maximum(m, 0.0)
    return m


def rank_two_tensor_q2_formula(
    projection: np.ndarray,
    *,
    atol: float = 1e-10,
) -> float:
    """Evaluate the exact tensor-Parseval formula for ``q2`` at rank two.

    Write ``d_i=P_ii``, ``D_k=sum_i d_i^k`` and

        S = sum_{i<j} (d_i d_j + 2 P_ij^2)^2.

    For every rank-two orthogonal projection on R^6,

        q2 = 1 - (5/2) D2 + 3 D3 - (9/8) D4 + S/4.

    The identity follows by applying Parseval to the fourth Gaussian moment
    tensor of the underlying two-dimensional Parseval frame and partitioning
    ordered four-tuples by collision pattern.
    """
    P = _validate_rank_two_projection(projection, atol=atol)
    d = np.diag(P)
    D2 = float(np.sum(d**2))
    D3 = float(np.sum(d**3))
    D4 = float(np.sum(d**4))
    S = 0.0
    for i, j in SIX_VARIABLE_EDGES:
        S += float((d[i] * d[j] + 2.0 * P[i, j] ** 2) ** 2)
    return 1.0 - 2.5 * D2 + 3.0 * D3 - 1.125 * D4 + 0.25 * S


def rank_two_projection_gradient_defect(
    projection: np.ndarray,
    *,
    atol: float = 1e-10,
) -> float:
    """Evaluate the exact rank-two gradient defect ``q1/4-q2``.

    Let ``m_ij`` be the squared Pluecker coordinates and set

        Dk = sum_i d_i^k,
        C  = sum_{i<j} d_i d_j m_ij,
        M2 = sum_{i<j} m_ij^2.

    Then

        8 (q1/4-q2)
          = -6 + 19 D2 - 24 D3 + 18 D4 - 9 D2^2
            + 24 C - 8 M2.

    This is an identity, not a claim that the defect is globally nonnegative.
    """
    P = _validate_rank_two_projection(projection, atol=atol)
    d = np.diag(P)
    m = rank_two_plucker_weights(P, atol=atol)
    D2 = float(np.sum(d**2))
    D3 = float(np.sum(d**3))
    D4 = float(np.sum(d**4))
    C = 0.0
    M2 = 0.0
    for i, j in SIX_VARIABLE_EDGES:
        C += float(d[i] * d[j] * m[i, j])
        M2 += float(m[i, j] ** 2)
    return (
        -6.0
        + 19.0 * D2
        - 24.0 * D3
        + 18.0 * D4
        - 9.0 * D2 * D2
        + 24.0 * C
        - 8.0 * M2
    ) / 8.0


def rank_two_projection_direct_defect(projection: np.ndarray) -> float:
    """Return ``q1/4-q2`` by direct hafnian evaluation for cross-checking."""
    q1, q2 = six_variable_gradient_energies(projection)
    return 0.25 * q1 - q2
