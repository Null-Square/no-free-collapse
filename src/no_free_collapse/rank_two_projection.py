"""Exact rank-two projection identities for the six-variable gradient problem."""

from __future__ import annotations

import math

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


def _rank_two_parseval_frame(projection: np.ndarray) -> np.ndarray:
    """Return one orthonormal-column frame U with P=UU^T."""
    eigenvalues, eigenvectors = np.linalg.eigh(projection)
    order = np.argsort(eigenvalues)[-2:]
    return eigenvectors[:, order]


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


def rank_two_harmonic_invariants(
    projection: np.ndarray,
    *,
    atol: float = 1e-10,
) -> tuple[float, float, float, complex, complex, complex]:
    """Return the exact planar Fourier invariants ``(D2,D3,D4,A,B,S)``.

    Choose any Parseval frame ``P=UU^T`` and write its nonzero row vectors as

        u_i = sqrt(d_i) (cos(theta_i), sin(theta_i)),
        z_i = exp(2 i theta_i).

    Parseval isotropy is exactly ``sum_i d_i z_i = 0``.  Define

        A = sum_i d_i^2 z_i,
        B = sum_i d_i^2 z_i^2,
        S = (1/2) sum_i d_i z_i^2.

    A latent-plane rotation changes only the phases of ``A,B,S``; their
    magnitudes and the harmonic defect formula are basis independent.
    """
    P = _validate_rank_two_projection(projection, atol=atol)
    U = _rank_two_parseval_frame(P)
    d = np.diag(P)
    w = U[:, 0].astype(np.complex128) + 1j * U[:, 1]
    z = np.ones(6, dtype=np.complex128)
    active = d > atol
    z[active] = (w[active] ** 2) / d[active]

    isotropy = np.sum(d * z)
    if abs(isotropy) > 100.0 * atol:
        raise ValueError("internal Parseval harmonic check failed")

    D2 = float(np.sum(d**2))
    D3 = float(np.sum(d**3))
    D4 = float(np.sum(d**4))
    A = complex(np.sum((d**2) * z))
    B = complex(np.sum((d**2) * (z**2)))
    S = complex(0.5 * np.sum(d * (z**2)))
    return D2, D3, D4, A, B, S


def rank_two_harmonic_gradient_defect(
    projection: np.ndarray,
    *,
    atol: float = 1e-10,
) -> float:
    """Evaluate the exact harmonic formula for ``q1/4-q2`` at rank two.

    With ``(D2,D3,D4,A,B,S)`` from :func:`rank_two_harmonic_invariants`,

        8 (q1/4-q2)
          = -6 + 19 D2 - 24 D3 + 18 D4
            - (9/2) D2^2 - 4 |A|^2 - (1/2) |B|^2.

    The auxiliary second harmonic ``S`` does not enter this identity directly;
    it is useful because ``|S|=1`` characterizes the exact two-direction
    endpoint of the planar moment problem.
    """
    D2, D3, D4, A, B, _ = rank_two_harmonic_invariants(projection, atol=atol)
    return (
        -6.0
        + 19.0 * D2
        - 24.0 * D3
        + 18.0 * D4
        - 4.5 * D2 * D2
        - 4.0 * abs(A) ** 2
        - 0.5 * abs(B) ** 2
    ) / 8.0


def rank_two_heavy_coordinate_decomposition(
    projection: np.ndarray,
    *,
    index: int | None = None,
    atol: float = 1e-10,
) -> tuple[float, np.ndarray, np.ndarray, int]:
    """Return the canonical heavy-coordinate interpolation ``(eps,p,q,index)``.

    If ``t=P[index,index]`` lies strictly between zero and one, set
    ``eps=1-t``.  After choosing the latent basis so the selected row is on the
    first latent axis, every rank-two projection has the exact block form

        P = [[1-eps, sqrt(eps(1-eps)) p^T],
             [sqrt(eps(1-eps)) p, eps p p^T + q q^T]],

    where ``p,q`` are orthonormal vectors in ``R^5``.  If ``index`` is omitted,
    the largest diagonal entry is used.  The intended high-leverage regime is
    ``P[index,index] > 1/2``; the identity itself holds for any nondegenerate
    selected coordinate.
    """
    P = _validate_rank_two_projection(projection, atol=atol)
    d = np.diag(P)
    if index is None:
        index = int(np.argmax(d))
    if not 0 <= index < 6:
        raise ValueError("index must be between 0 and 5")

    t = float(d[index])
    eps = 1.0 - t
    if t <= atol or eps <= atol:
        raise ValueError("selected coordinate must have diagonal strictly between 0 and 1")

    others = [j for j in range(6) if j != index]
    scale = math.sqrt(t * eps)
    p = np.asarray(P[index, others] / scale, dtype=np.float64)
    block = P[np.ix_(others, others)]
    qq = block - eps * np.outer(p, p)
    eigenvalues, eigenvectors = np.linalg.eigh(qq)
    q = np.asarray(eigenvectors[:, int(np.argmax(eigenvalues))], dtype=np.float64)
    q /= np.linalg.norm(q)

    if not np.isclose(p @ p, 1.0, atol=100.0 * atol, rtol=100.0 * atol):
        raise ValueError("internal heavy-coordinate p normalization failed")
    if abs(float(p @ q)) > 100.0 * atol:
        raise ValueError("internal heavy-coordinate orthogonality failed")
    if not np.allclose(qq, np.outer(q, q), atol=200.0 * atol, rtol=200.0 * atol):
        raise ValueError("internal heavy-coordinate rank-one remainder failed")
    return eps, p, q, index


def rank_two_heavy_coordinate_projection(
    eps: float,
    p: np.ndarray,
    q: np.ndarray,
    *,
    atol: float = 1e-10,
) -> np.ndarray:
    """Construct the canonical rank-two projection from ``eps,p,q``.

    ``p`` and ``q`` must be orthonormal vectors in ``R^5`` and ``0<=eps<=1``.
    The first coordinate has leverage ``1-eps``.
    """
    eps = float(eps)
    p = np.asarray(p, dtype=np.float64)
    q = np.asarray(q, dtype=np.float64)
    if p.shape != (5,) or q.shape != (5,):
        raise ValueError("p and q must have shape (5,)")
    if eps < -atol or eps > 1.0 + atol:
        raise ValueError("eps must lie in [0,1]")
    eps = min(1.0, max(0.0, eps))
    if not np.isclose(p @ p, 1.0, atol=atol, rtol=atol):
        raise ValueError("p must be a unit vector")
    if not np.isclose(q @ q, 1.0, atol=atol, rtol=atol):
        raise ValueError("q must be a unit vector")
    if abs(float(p @ q)) > atol:
        raise ValueError("p and q must be orthogonal")

    P = np.zeros((6, 6), dtype=np.float64)
    P[0, 0] = 1.0 - eps
    star = math.sqrt(max(0.0, eps * (1.0 - eps))) * p
    P[0, 1:] = star
    P[1:, 0] = star
    P[1:, 1:] = eps * np.outer(p, p) + np.outer(q, q)
    return P


def rank_two_projection_direct_defect(projection: np.ndarray) -> float:
    """Return ``q1/4-q2`` by direct hafnian evaluation for cross-checking."""
    q1, q2 = six_variable_gradient_energies(projection)
    return 0.25 * q1 - q2
