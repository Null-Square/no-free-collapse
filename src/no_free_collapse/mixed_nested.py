"""Exact mixed nested-projection formulas for the contraction proof program."""

from __future__ import annotations

import numpy as np


def _validate_orthonormal_pair(
    u: np.ndarray,
    v: np.ndarray,
    *,
    atol: float = 1e-10,
) -> tuple[np.ndarray, np.ndarray]:
    u = np.asarray(u, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    if u.shape != (6,) or v.shape != (6,):
        raise ValueError("u and v must have shape (6,)")
    if not np.isclose(u @ u, 1.0, atol=atol, rtol=atol):
        raise ValueError("u must be a unit vector")
    if not np.isclose(v @ v, 1.0, atol=atol, rtol=atol):
        raise ValueError("v must be a unit vector")
    if abs(float(u @ v)) > atol:
        raise ValueError("u and v must be orthogonal")
    return u, v


def rank_one_rank_five_phi(t: np.ndarray | float, mu: float) -> np.ndarray | float:
    """Return the diagonal polynomial in the `(1,5,5,5)` mixed coefficient.

    The polynomial is

        phi_mu(t) = 1-mu + (6 mu-5)t + 24 t^2 - 36 t^3.
    """
    values = np.asarray(t, dtype=np.float64)
    result = 1.0 - float(mu) + (6.0 * float(mu) - 5.0) * values + 24.0 * values**2 - 36.0 * values**3
    if np.ndim(values) == 0:
        return float(result)
    return result


def rank_one_rank_five_mixed_formula(
    u: np.ndarray,
    v: np.ndarray,
    *,
    atol: float = 1e-10,
) -> float:
    """Evaluate the exact nested `(1,5,5,5)` homogenized defect coefficient.

    Let

        P = u u^T,
        Q = I - v v^T,
        u perpendicular to v.

    For the symmetric homogenized spectral kernel `D` used in the contraction
    reduction, the exact identity is

        16 D(P,Q,Q,Q)
          = sum_i phi_mu(t_i) u_i^2
            + 12 (sum_i u_i v_i^3)^2,

    where `t_i=v_i^2`, `mu=sum_i t_i^2`, and `phi_mu` is returned by
    :func:`rank_one_rank_five_phi`.

    The accompanying proof note establishes that this quantity is always
    nonnegative.  Simultaneous complementation gives the `(1,1,1,5)` case.
    """
    u, v = _validate_orthonormal_pair(u, v, atol=atol)
    t = v * v
    mu = float(t @ t)
    phi = np.asarray(rank_one_rank_five_phi(t, mu), dtype=np.float64)
    return float((phi @ (u * u) + 12.0 * float(u @ (v**3)) ** 2) / 16.0)
