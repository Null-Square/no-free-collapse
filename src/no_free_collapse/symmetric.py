"""Exactly solvable permutation-symmetric order-1 constructions."""

from __future__ import annotations

import math

import numpy as np


def mean_field_gram(n: int, kappa: float) -> np.ndarray:
    """Return an order-1 Gram matrix with q(x)=1+t(sum_i x_i)^2."""
    if n < 2 or n % 2:
        raise ValueError("n must be an even integer at least two")
    if kappa < 1:
        raise ValueError("condition number must be at least one")
    t = (kappa - 1.0) / (n * n)
    Q = np.zeros((n + 1, n + 1), dtype=np.float64)
    Q[1:, 1:] = np.eye(n) / n + t * np.ones((n, n))
    return Q


def mean_field_full_parity_log_capacity(n: int, kappa: float) -> float:
    """Log of the exact optimal full-parity coefficient for mean_field_gram."""
    if n < 2 or n % 2:
        raise ValueError("n must be an even integer at least two")
    if kappa < 1:
        raise ValueError("condition number must be at least one")
    if kappa == 1:
        return -math.inf
    beta_sq = n * n / (4.0 * (kappa - 1.0))
    log_denominator = sum(math.log(m * m + beta_sq) for m in range(1, n // 2 + 1))
    return math.log(kappa / (kappa - 1.0)) + math.lgamma(n + 1.0) - n * math.log(2.0) - log_denominator


def mean_field_full_parity_capacity(n: int, kappa: float) -> float:
    """Exact optimal full-parity coefficient for the mean-field construction."""
    log_value = mean_field_full_parity_log_capacity(n, kappa)
    return 0.0 if log_value == -math.inf else math.exp(log_value)


def mean_field_asymptotic_rate(kappa: float) -> float:
    """Rate I(kappa) with C_n=exp(-n I(kappa)+O(log n)) for fixed kappa."""
    if kappa <= 1:
        if kappa == 1:
            return math.inf
        raise ValueError("condition number must be at least one")
    root = math.sqrt(kappa - 1.0)
    return 0.5 * math.log(kappa / (kappa - 1.0)) + math.atan(root) / root
