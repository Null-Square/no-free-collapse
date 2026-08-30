"""Exact symmetry-breaking constructions from matched variable pairs."""

from __future__ import annotations

import math

import numpy as np


def _validate(n: int, kappa: float) -> int:
    if n < 4 or n % 2:
        raise ValueError("n must be an even integer at least four")
    if kappa < 1:
        raise ValueError("condition number must be at least one")
    return n // 2


def matched_pair_gram(n: int, kappa: float) -> np.ndarray:
    """Optimal Gram matrix in the matched pair-block class.

    Feature order is [1,x1,...,xn]. One isotropic readout pair supplies the
    baseline norm. Each remaining pair is rank one and carries an equal share
    of the normalization variation. The cube norm condition number is kappa.
    """
    L = _validate(n, kappa)
    Q = np.zeros((n + 1, n + 1), dtype=np.float64)
    Q[1, 1] = Q[2, 2] = 1.0
    if kappa == 1:
        return Q
    d = (kappa - 1.0) / (L - 1)
    for pair in range(1, L):
        i = 1 + 2 * pair
        Q[i, i] = Q[i + 1, i + 1] = 0.5 * d
        Q[i, i + 1] = Q[i + 1, i] = -0.5 * d
    return Q


def matched_pair_log_capacity(n: int, kappa: float) -> float:
    """Log exact full-parity capacity of the optimal matched-pair family."""
    L = _validate(n, kappa)
    if kappa == 1:
        return -math.inf
    alpha = (L - 1.0) / (kappa - 1.0)
    return -L * math.log(2.0) + math.lgamma(1.0 + alpha) + math.lgamma(L) - math.lgamma(L + alpha)


def matched_pair_capacity(n: int, kappa: float) -> float:
    """Exact optimal capacity within the matched pair-block class."""
    log_value = matched_pair_log_capacity(n, kappa)
    return 0.0 if log_value == -math.inf else math.exp(log_value)


def matched_pair_asymptotic_rate(kappa: float) -> float:
    """Fixed-kappa exponent I_pair with C_pair=exp(-n I_pair+O(log n))."""
    if kappa < 1:
        raise ValueError("condition number must be at least one")
    if kappa == 1:
        return math.inf
    c = 1.0 / (kappa - 1.0)
    return 0.5 * (math.log(2.0) + (1.0 + c) * math.log(1.0 + c) - c * math.log(c))
