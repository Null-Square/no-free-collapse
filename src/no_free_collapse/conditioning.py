"""Conditioning-controlled bounds for input-normalized Born readout.

For an unnormalized polynomial amplitude preparation of Boolean order r,
normalization produces a rational function p=a/q.  This module quantifies how
far that rational function can depart from low-degree polynomial behavior when
the squared norm q is well-conditioned over the Boolean hypercube.
"""

from __future__ import annotations


def condition_number(q_min: float, q_max: float) -> float:
    """Return kappa=q_max/q_min for strictly positive squared norms."""
    if q_min <= 0:
        raise ValueError("q_min must be strictly positive")
    if q_max < q_min:
        raise ValueError("q_max must be at least q_min")
    return q_max / q_min


def centered_relative_variation(q_min: float, q_max: float) -> float:
    """Return the optimal centered relative variation delta in [0,1).

    With q_c=(q_max+q_min)/2 and h=q/q_c-1, every q in [q_min,q_max]
    satisfies |h|<=delta, where

        delta=(q_max-q_min)/(q_max+q_min).
    """
    if q_min <= 0:
        raise ValueError("q_min must be strictly positive")
    if q_max < q_min:
        raise ValueError("q_max must be at least q_min")
    return (q_max - q_min) / (q_max + q_min)


def delta_from_condition_number(kappa: float) -> float:
    """Convert kappa>=1 to delta=(kappa-1)/(kappa+1)."""
    if kappa < 1:
        raise ValueError("condition number must be at least one")
    return (kappa - 1.0) / (kappa + 1.0)


def truncation_bound(
    preparation_order: int,
    steps: int,
    kappa: float,
) -> tuple[int, float]:
    """Degree/error pair for geometric de-normalization.

    If an order-r amplitude preparation has squared-norm condition number
    kappa and 0<=M<=I is a fixed Born effect, then its normalized probability
    p can be uniformly approximated by a polynomial of degree at most 2*r*m
    with error at most delta**m, where m=steps and
    delta=(kappa-1)/(kappa+1).
    """
    if preparation_order < 0:
        raise ValueError("preparation_order must be non-negative")
    if steps < 0:
        raise ValueError("steps must be non-negative")
    delta = delta_from_condition_number(kappa)
    return 2 * preparation_order * steps, delta**steps


def spectral_coefficient_bound(
    interaction_order: int,
    preparation_order: int,
    kappa: float,
) -> float:
    """Bound any Walsh coefficient at the requested interaction order.

    For k>0 and r>0, choose the largest m with 2*r*m<k.  The degree-2*r*m
    truncation has no order-k component, so |f_hat(S)|<=delta**m.
    """
    if interaction_order < 0:
        raise ValueError("interaction_order must be non-negative")
    if preparation_order < 0:
        raise ValueError("preparation_order must be non-negative")
    if preparation_order == 0:
        return 0.0 if interaction_order > 0 else 1.0

    steps = (interaction_order - 1) // (2 * preparation_order)
    if steps <= 0:
        return 1.0
    return delta_from_condition_number(kappa) ** steps


def spectral_tail_energy_bound(
    cutoff_degree: int,
    preparation_order: int,
    kappa: float,
) -> float:
    """Bound Fourier energy strictly above a degree cutoff.

    The guarantee is strongest when cutoff_degree is a multiple of 2*r.  We
    use the largest m such that 2*r*m<=cutoff_degree and obtain

        sum_{|S|>cutoff} |f_hat(S)|^2 <= delta**(2m).
    """
    if cutoff_degree < 0:
        raise ValueError("cutoff_degree must be non-negative")
    if preparation_order < 0:
        raise ValueError("preparation_order must be non-negative")
    if preparation_order == 0:
        return 0.0

    steps = cutoff_degree // (2 * preparation_order)
    return delta_from_condition_number(kappa) ** (2 * steps)


def required_condition_number_for_coefficient(
    coefficient_magnitude: float,
    interaction_order: int,
    preparation_order: int,
) -> float:
    """Invert the spectral bound to lower-bound norm conditioning.

    If an order-k Walsh coefficient has magnitude alpha and
    m=floor((k-1)/(2r))>=1, then delta^m>=alpha and therefore

        kappa >= (1+alpha^(1/m))/(1-alpha^(1/m)).
    """
    if not (0.0 <= coefficient_magnitude < 1.0):
        raise ValueError("coefficient_magnitude must lie in [0,1)")
    if interaction_order < 0:
        raise ValueError("interaction_order must be non-negative")
    if preparation_order <= 0:
        raise ValueError("preparation_order must be positive")

    steps = (interaction_order - 1) // (2 * preparation_order)
    if steps <= 0 or coefficient_magnitude == 0:
        return 1.0

    delta = coefficient_magnitude ** (1.0 / steps)
    return (1.0 + delta) / (1.0 - delta)
