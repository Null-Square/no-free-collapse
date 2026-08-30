"""Small exact rational examples used as theorem witnesses and regression tests."""

from __future__ import annotations

from fractions import Fraction
from typing import Sequence


def normalization_loophole_probability(x: Sequence[int]) -> Fraction:
    """Normalized probability for psi_tilde=(2+x1+x2+x3, 1).

    This amplitude preparation is affine/order-1 before normalization, but its
    normalized Born probability has nonzero order-3 Walsh coefficient -6/65.
    """
    if len(x) != 3 or any(xi not in (-1, 1) for xi in x):
        raise ValueError("x must contain exactly three +/-1 values")
    a = 2 + sum(x)
    return Fraction(a * a, a * a + 1)
