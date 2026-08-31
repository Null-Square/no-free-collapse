"""Diagonal-only dual certificates for the six-variable rank-two problem."""

from __future__ import annotations

from math import comb

import numpy as np


def _validate_rank_two_diagonal(diagonal: np.ndarray, *, atol: float = 1e-10) -> np.ndarray:
    d = np.asarray(diagonal, dtype=np.float64)
    if d.shape != (6,):
        raise ValueError("diagonal must have shape (6,)")
    if np.min(d) < -atol or np.max(d) > 1.0 + atol:
        raise ValueError("diagonal entries must lie in [0,1]")
    if not np.isclose(np.sum(d), 2.0, atol=atol, rtol=0.0):
        raise ValueError("rank-two projection diagonals must sum to two")
    return np.clip(d, 0.0, 1.0)


def rank_two_quadratic_dual_certificate(
    diagonal: np.ndarray,
    *,
    atol: float = 1e-10,
) -> float:
    """Return the explicit quadratic-dual lower certificate ``H(d)``.

    Let ``d`` be any vector in ``[0,1]^6`` with ``sum d_i=2``.  For an actual
    rank-two projection with this diagonal, let

        Delta = q1/4 - q2.

    The exact Pluecker defect, chord linearization on

        max(0,d_i+d_j-1) <= m_ij <= d_i d_j,

    and the node potential

        alpha(t) = -t^2 + (3/2)t - 1/4

    give the rigorous lower bound

        8 Delta >= H(d).

    Writing ``r_i=d_i(1-d_i)``, ``R=sum r_i`` and ``E=sum r_i^2``, the
    certificate has the compact form

        H = R(1-R) + 2E - C_low - C_high,

    where

        C_low  = sum_{d_i+d_j<1/2}
                 4 d_i d_j (1-d_i-d_j) (1-2d_i-2d_j),

        C_high = sum_{d_i+d_j>3/2}
                 4 (1-d_i)(1-d_j)(d_i+d_j-1)(2d_i+2d_j-3).

    This routine evaluates the proved dual certificate.  Global
    nonnegativity is currently proved on the high-pair regime
    ``max_{i<j}(d_i+d_j) >= 3/2``; the complementary regime is the remaining
    diagonal-only target.
    """
    d = _validate_rank_two_diagonal(diagonal, atol=atol)
    r = d * (1.0 - d)
    R = float(np.sum(r))
    E = float(r @ r)

    correction = 0.0
    for i in range(6):
        for j in range(i + 1, 6):
            x = float(d[i])
            y = float(d[j])
            s = x + y
            if s < 0.5:
                correction += 4.0 * x * y * (1.0 - s) * (1.0 - 2.0 * s)
            elif s > 1.5:
                correction += (
                    4.0
                    * (1.0 - x)
                    * (1.0 - y)
                    * (s - 1.0)
                    * (2.0 * s - 3.0)
                )
    return R * (1.0 - R) + 2.0 * E - correction


def rank_two_high_pair_bernstein_coefficients(
    diagonal: np.ndarray,
    *,
    atol: float = 1e-10,
) -> tuple[float, np.ndarray]:
    """Return ``(x,beta)`` for the high-pair dual certificate.

    Sort ``d`` decreasingly and assume ``d_0+d_1>=3/2``.  Set

        a = 1-d_0,  b = 1-d_1,
        T = a+b = d_2+...+d_5 <= 1/2,
        x = 2T.

    For ``T>0`` normalize ``u=(a,b)/T`` and ``v=(d_2,...,d_5)/T``.  Then

        H(d) = sum_{k=0}^4 beta_k C(4,k) x^k (1-x)^(4-k),

    with the exact coefficients used in the analytic high-pair proof.
    When ``T=0`` the diagonal is a permutation of ``(1,1,0,0,0,0)`` and the
    returned coefficient vector is ``(0,1/4,25/96,1/32,127/512)``; only
    ``beta_0`` contributes at ``x=0``.
    """
    d = np.sort(_validate_rank_two_diagonal(diagonal, atol=atol))[::-1]
    if d[0] + d[1] < 1.5 - atol:
        raise ValueError("largest pair must have sum at least 3/2")

    a = 1.0 - float(d[0])
    b = 1.0 - float(d[1])
    T = a + b
    x = 2.0 * T
    if T <= atol:
        return 0.0, np.asarray([0.0, 0.25, 25.0 / 96.0, 1.0 / 32.0, 127.0 / 512.0])

    u = np.asarray([a, b], dtype=np.float64) / T
    v = np.asarray(d[2:], dtype=np.float64) / T
    if not np.isclose(np.sum(u), 1.0, atol=100 * atol, rtol=0.0):
        raise ValueError("internal high-pair two-block normalization failed")
    if not np.isclose(np.sum(v), 1.0, atol=100 * atol, rtol=0.0):
        raise ValueError("internal high-pair four-block normalization failed")

    A2 = float(np.sum(u**2))
    A3 = float(np.sum(u**3))
    A4 = float(np.sum(u**4))
    Z2 = float(np.sum(v**2))
    Z3 = float(np.sum(v**3))
    Z4 = float(np.sum(v**4))
    S2 = A2 + Z2
    S3 = A3 + Z3
    S4 = A4 + Z4

    beta = np.zeros(5, dtype=np.float64)
    beta[0] = 0.0
    beta[1] = 0.25
    beta[2] = 1.0 / 6.0 + S2 / 8.0
    beta[3] = -0.25 + 7.0 * S2 / 8.0 - S3 / 2.0
    beta[4] = (
        -1.0
        + 11.0 * S2 / 4.0
        - 5.0 * S3 / 2.0
        + 9.0 * S4 / 8.0
        - 9.0 * S2 * S2 / 16.0
        + A2 * Z2
    )
    return x, beta


def evaluate_degree_four_bernstein(x: float, coefficients: np.ndarray) -> float:
    """Evaluate a degree-four Bernstein expansion at ``x``."""
    coefficients = np.asarray(coefficients, dtype=np.float64)
    if coefficients.shape != (5,):
        raise ValueError("coefficients must have shape (5,)")
    return float(
        sum(
            coefficients[k] * comb(4, k) * x**k * (1.0 - x) ** (4 - k)
            for k in range(5)
        )
    )
