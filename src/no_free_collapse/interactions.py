"""Walsh/Fourier interaction analysis on the Boolean hypercube {-1,+1}^n."""

from __future__ import annotations

from itertools import product
from typing import Callable, Dict, Iterable, Mapping, Sequence, TypeVar

Scalar = TypeVar("Scalar")


def hypercube(n: int) -> Iterable[tuple[int, ...]]:
    """Enumerate {-1,+1}^n in deterministic order."""
    if n < 0:
        raise ValueError("n must be non-negative")
    return product((-1, 1), repeat=n)


def monomial(x: Sequence[int], mask: int):
    """Evaluate chi_mask(x)=prod_{i in mask} x_i."""
    out = 1
    for i, xi in enumerate(x):
        if mask & (1 << i):
            out *= xi
    return out


def evaluate_on_hypercube(n: int, fn: Callable[[tuple[int, ...]], Scalar]) -> Dict[tuple[int, ...], Scalar]:
    return {x: fn(x) for x in hypercube(n)}


def walsh_coefficients(values: Mapping[tuple[int, ...], Scalar], n: int) -> Dict[int, Scalar]:
    """Return exact Walsh coefficients when the scalar type supports exact arithmetic.

    The multilinear expansion is f(x)=sum_A f_hat[A] prod_{i in A} x_i.
    """
    expected = 1 << n
    if len(values) != expected:
        raise ValueError(f"expected {expected} hypercube values, got {len(values)}")

    coeffs: Dict[int, Scalar] = {}
    scale = expected
    for mask in range(expected):
        total = None
        for x, value in values.items():
            term = value * monomial(x, mask)
            total = term if total is None else total + term
        coeffs[mask] = total / scale
    return coeffs


def mask_order(mask: int) -> int:
    return mask.bit_count()


def interaction_degree(coeffs: Mapping[int, Scalar], *, atol: float = 1e-10) -> int:
    """Largest non-negligible Walsh interaction order."""
    degree = 0
    for mask, value in coeffs.items():
        try:
            nonzero = abs(value) > atol
        except TypeError:
            nonzero = value != 0
        if nonzero:
            degree = max(degree, mask_order(mask))
    return degree


def max_abs_by_order(coeffs: Mapping[int, complex]) -> Dict[int, float]:
    out: Dict[int, float] = {}
    for mask, value in coeffs.items():
        order = mask_order(mask)
        out[order] = max(out.get(order, 0.0), float(abs(value)))
    return out
