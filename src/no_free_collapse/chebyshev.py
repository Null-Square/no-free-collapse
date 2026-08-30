"""Chebyshev-minimax bounds for conditioned Born interaction capacity."""

from __future__ import annotations

import math


def chebyshev_rho(kappa: float) -> float:
    """Return rho=(sqrt(kappa)-1)/(sqrt(kappa)+1) in [0,1)."""
    if kappa < 1:
        raise ValueError("condition number must be at least one")
    if kappa == 1:
        return 0.0
    root = math.sqrt(kappa)
    return (root - 1.0) / (root + 1.0)


def chebyshev_residual_bound(steps: int, kappa: float) -> float:
    """Minimax residual for polynomials R with deg(R)<=steps and R(0)=1."""
    if steps < 0:
        raise ValueError("steps must be non-negative")
    if kappa < 1:
        raise ValueError("condition number must be at least one")
    if steps == 0:
        return 1.0
    if kappa == 1:
        return 0.0
    rho = chebyshev_rho(kappa)
    power = rho**steps
    return 2.0 * power / (1.0 + power * power)


def cancellation_steps(interaction_order: int, preparation_order: int) -> int:
    """Largest m for which 2*r*m is strictly below the target order."""
    if interaction_order < 0:
        raise ValueError("interaction_order must be non-negative")
    if preparation_order <= 0:
        raise ValueError("preparation_order must be positive")
    return (interaction_order - 1) // (2 * preparation_order)


def scalar_chebyshev_interaction_bound(interaction_order: int, preparation_order: int, kappa: float) -> float:
    """Uniform approximation bound for a normalized Born probability coefficient."""
    steps = cancellation_steps(interaction_order, preparation_order)
    if steps <= 0:
        return 1.0
    return chebyshev_residual_bound(steps, kappa)


def born_chebyshev_interaction_bound(interaction_order: int, preparation_order: int, kappa: float) -> float:
    """Physical Born bound using Helstrom trace distance."""
    if interaction_order == 0:
        return 1.0
    steps = cancellation_steps(interaction_order, preparation_order)
    if steps <= 0:
        return 0.5
    return 0.5 * chebyshev_residual_bound(steps, kappa)


def asymptotic_chebyshev_rate(kappa: float, preparation_order: int = 1) -> float:
    """Exponential rate per interaction variable of the large-order bound."""
    if preparation_order <= 0:
        raise ValueError("preparation_order must be positive")
    rho = chebyshev_rho(kappa)
    if rho == 0:
        return math.inf
    return -math.log(rho) / (2.0 * preparation_order)
