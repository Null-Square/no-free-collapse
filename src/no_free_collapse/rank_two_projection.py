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
    """Return squared Pluecker coordinates of a rank-two projection."""
    P = _validate_rank_two_projection(projection, atol=atol)
    d = np.diag(P)
    m = np.outer(d, d) - P * P
    np.fill_diagonal(m, 0.0)
    if float(np.min(m)) < -20.0 * atol:
        raise ValueError("rank-two principal minors must be nonnegative")
    return np.maximum(m, 0.0)


def rank_two_tensor_q2_formula(
    projection: np.ndarray,
    *,
    atol: float = 1e-10,
) -> float:
    """Evaluate the exact tensor-Parseval formula for ``q2`` at rank two."""
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
    """Evaluate the exact rank-two gradient defect ``q1/4-q2``."""
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


def rank_two_stability_decomposition(
    projection: np.ndarray,
    *,
    atol: float = 1e-10,
) -> tuple[float, float, float]:
    """Return ``(F, commutator_energy, angle_mixing)`` in the exact defect split.

    For every rank-two projection, with ``D=diag(diag(P))`` and squared
    Pluecker weights ``m_ij``, the exact identity is

        8 (q1/4-q2) = F + 8 E_comm + 8 J,

    where

        F = -6 + 19 D2 - 32 D3 + 18 D4 - D2^2,
        E_comm = ||(I-P) D P||_F^2,
        J = sum_{i<j} m_ij P_ij^2.

    Both energies are nonnegative.  The balanced-diagonal theorem proves
    ``F>=0`` when every diagonal entry is at most one half; outside that
    region the two energies quantify exactly how geometry repairs a possibly
    negative scalar term.
    """
    P = _validate_rank_two_projection(projection, atol=atol)
    d = np.diag(P)
    D = np.diag(d)
    D2 = float(np.sum(d**2))
    D3 = float(np.sum(d**3))
    D4 = float(np.sum(d**4))
    F = -6.0 + 19.0 * D2 - 32.0 * D3 + 18.0 * D4 - D2 * D2
    E_comm = float(np.linalg.norm((np.eye(6) - P) @ D @ P, ord="fro") ** 2)
    m = rank_two_plucker_weights(P, atol=atol)
    J = 0.0
    for i, j in SIX_VARIABLE_EDGES:
        J += float(m[i, j] * P[i, j] ** 2)
    return F, E_comm, J


def rank_two_harmonic_invariants(
    projection: np.ndarray,
    *,
    atol: float = 1e-10,
) -> tuple[float, float, float, complex, complex, complex]:
    """Return the exact planar Fourier invariants ``(D2,D3,D4,A,B,S)``."""
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
    """Evaluate the exact harmonic formula for ``q1/4-q2`` at rank two."""
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
    """Return the canonical heavy-coordinate interpolation ``(eps,p,q,index)``."""
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
    """Construct the canonical rank-two projection from ``eps,p,q``."""
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


def rank_two_heavy_bernstein_coefficients(
    p: np.ndarray,
    q: np.ndarray,
    *,
    atol: float = 1e-10,
) -> np.ndarray:
    """Return the degree-four Bernstein coefficients of the heavy-path defect.

    Let ``x=2 eps`` and ``Delta(x)=q1(P(x/2))/4-q2(P(x/2))``.  Then

        Delta(x) = sum_{k=0}^4 b_k C(4,k) x^k (1-x)^(4-k).

    This helper returns ``(b0,...,b4)``.  The first interior coefficient
    ``b1`` is proved nonnegative for every orthonormal ``p,q``; the global
    positivity of ``b2,b3`` remains the final interpolation target.
    """
    p = np.asarray(p, dtype=np.float64)
    q = np.asarray(q, dtype=np.float64)
    if p.shape != (5,) or q.shape != (5,):
        raise ValueError("p and q must have shape (5,)")
    if not np.isclose(p @ p, 1.0, atol=atol, rtol=atol):
        raise ValueError("p must be a unit vector")
    if not np.isclose(q @ q, 1.0, atol=atol, rtol=atol):
        raise ValueError("q must be a unit vector")
    if abs(float(p @ q)) > atol:
        raise ValueError("p and q must be orthogonal")

    xs = np.asarray([0.0, 0.25, 0.5, 0.75, 1.0], dtype=np.float64)
    values = np.asarray(
        [
            rank_two_projection_direct_defect(
                rank_two_heavy_coordinate_projection(0.5 * x, p, q, atol=atol)
            )
            for x in xs
        ],
        dtype=np.float64,
    )
    power = np.linalg.solve(np.vander(xs, 5, increasing=True), values)
    bernstein = np.zeros(5, dtype=np.float64)
    for k in range(5):
        for j in range(k + 1):
            bernstein[k] += (
                power[j] * math.comb(k, j) / math.comb(4, j)
            )
    return bernstein


def rank_two_projection_direct_defect(projection: np.ndarray) -> float:
    """Return ``q1/4-q2`` by direct hafnian evaluation for cross-checking."""
    q1, q2 = six_variable_gradient_energies(projection)
    return 0.25 * q1 - q2
