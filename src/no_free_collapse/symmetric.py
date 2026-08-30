"""Exactly solvable permutation-symmetric order-1 constructions."""

from __future__ import annotations

import math

import numpy as np


def _validate_even(n: int, kappa: float) -> None:
    if n < 2 or n % 2:
        raise ValueError("n must be an even integer at least two")
    if kappa < 1:
        raise ValueError("condition number must be at least one")


def symmetric_branch_gram(
    n: int,
    kappa: float,
    constant_fraction: float = 0.0,
    *,
    increasing: bool = True,
) -> np.ndarray:
    """Permutation- and global-sign-invariant order-1 Gram family.

    The feature ordering is [1,x1,...,xn].  ``constant_fraction`` is
    lambda=a/A, where A=q(0) for the associated quadratic norm profile.

    Increasing branch:
        q(x)=A[1+(kappa-1)(sum_i x_i)^2/n^2], 0<=lambda<=1.

    Decreasing branch:
        q(x)=A[1-(kappa-1)(sum_i x_i)^2/(kappa n^2)],
        0<=lambda<=1/kappa.

    Overall scale A is fixed to one because Born probabilities and the cube
    condition number are invariant under positive rescaling of Q.
    """
    _validate_even(n, kappa)
    lam = float(constant_fraction)
    if increasing:
        if not (0.0 <= lam <= 1.0):
            raise ValueError("increasing branch requires 0<=constant_fraction<=1")
        h = (1.0 - lam) / n
        g = (kappa - lam) / n
    else:
        if kappa == 1:
            if lam != 0.0:
                raise ValueError("at kappa=1 only constant_fraction=0 is allowed")
            h = g = 1.0 / n
        else:
            if not (0.0 <= lam <= 1.0 / kappa):
                raise ValueError("decreasing branch requires 0<=constant_fraction<=1/kappa")
            h = (1.0 - lam) / n
            g = (1.0 / kappa - lam) / n

    Q = np.zeros((n + 1, n + 1), dtype=np.float64)
    Q[0, 0] = lam
    J = np.ones((n, n), dtype=np.float64) / n
    Q[1:, 1:] = h * (np.eye(n) - J) + g * J
    return Q


def mean_field_gram(n: int, kappa: float) -> np.ndarray:
    """The increasing symmetric optimum: lambda=0."""
    return symmetric_branch_gram(n, kappa, 0.0, increasing=True)


def mean_field_full_parity_log_capacity(n: int, kappa: float) -> float:
    """Log exact optimal full-parity coefficient for ``mean_field_gram``."""
    _validate_even(n, kappa)
    if kappa == 1:
        return -math.inf
    beta_sq = n * n / (4.0 * (kappa - 1.0))
    log_denominator = sum(math.log(m * m + beta_sq) for m in range(1, n // 2 + 1))
    return (
        math.log(kappa / (kappa - 1.0))
        + math.lgamma(n + 1.0)
        - n * math.log(2.0)
        - log_denominator
    )


def mean_field_full_parity_capacity(n: int, kappa: float) -> float:
    """Exact optimal full-parity coefficient for the mean-field construction."""
    log_value = mean_field_full_parity_log_capacity(n, kappa)
    return 0.0 if log_value == -math.inf else math.exp(log_value)


def increasing_branch_capacity(n: int, kappa: float, constant_fraction: float) -> float:
    """Exact capacity along the increasing symmetric branch.

    For lambda=constant_fraction,
        C(lambda)=C_mf * (kappa-lambda)/kappa.
    """
    _validate_even(n, kappa)
    lam = float(constant_fraction)
    if not (0.0 <= lam <= 1.0):
        raise ValueError("constant_fraction must lie in [0,1]")
    if kappa == 1:
        return 0.0
    return mean_field_full_parity_capacity(n, kappa) * (kappa - lam) / kappa


def decreasing_branch_full_parity_log_capacity(n: int, kappa: float) -> float:
    """Log capacity at lambda=0 on the decreasing symmetric branch."""
    _validate_even(n, kappa)
    if kappa == 1:
        return -math.inf
    beta_sq = kappa * n * n / (4.0 * (kappa - 1.0))
    log_denominator = sum(math.log(beta_sq - m * m) for m in range(1, n // 2 + 1))
    return (
        math.lgamma(n + 1.0)
        - n * math.log(2.0)
        - math.log(kappa - 1.0)
        - log_denominator
    )


def decreasing_branch_full_parity_capacity(n: int, kappa: float) -> float:
    """Maximum capacity on the decreasing symmetric branch (lambda=0)."""
    log_value = decreasing_branch_full_parity_log_capacity(n, kappa)
    return 0.0 if log_value == -math.inf else math.exp(log_value)


def decreasing_branch_capacity(n: int, kappa: float, constant_fraction: float) -> float:
    """Exact capacity along the decreasing symmetric branch.

    For lambda=constant_fraction,
        C(lambda)=C_minus(0)*(1-lambda).
    """
    _validate_even(n, kappa)
    lam = float(constant_fraction)
    if kappa == 1:
        return 0.0
    if not (0.0 <= lam <= 1.0 / kappa):
        raise ValueError("constant_fraction must lie in [0,1/kappa]")
    return decreasing_branch_full_parity_capacity(n, kappa) * (1.0 - lam)


def mean_field_asymptotic_rate(kappa: float) -> float:
    """Rate I(kappa) with C_n=exp(-n I(kappa)+O(log n)) for fixed kappa."""
    if kappa <= 1:
        if kappa == 1:
            return math.inf
        raise ValueError("condition number must be at least one")
    root = math.sqrt(kappa - 1.0)
    return 0.5 * math.log(kappa / (kappa - 1.0)) + math.atan(root) / root
