"""Exact rational checks for the averaged spectral hafnian hierarchy.

These certify supplied instances, not the universal analytic theorem. Both
spectral endpoint assumptions are verified; floats and bools are rejected.
The rational enumerator is capped at order 14 to bound exponential work.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import comb
from typing import Sequence

from .zero_diagonal import _fraction, _psd_exact


def exact_principal_hafnians(matrix: Sequence[Sequence[object]]) -> dict[int, Fraction]:
    """Evaluate all even principal hafnians over Q, with the empty value one."""
    a = [[_fraction(x) for x in row] for row in matrix]
    n = len(a)
    if not 1 <= n <= 14 or any(len(row) != n for row in a):
        raise ValueError("Exact enumeration requires square order 1 through 14")
    if any(a[i][j] != a[j][i] for i in range(n) for j in range(n)):
        raise ValueError("Expected exact symmetry")

    @lru_cache(maxsize=None)
    def haf(mask: int) -> Fraction:
        if not mask:
            return Fraction(1)
        bit = mask & -mask
        i = bit.bit_length()-1
        rest = mask ^ bit
        choices = rest
        result = Fraction(0)
        while choices:
            jbit = choices & -choices
            result += a[i][jbit.bit_length()-1] * haf(rest ^ jbit)
            choices ^= jbit
        return result

    return {mask: haf(mask) for mask in range(1 << n) if mask.bit_count() % 2 == 0}


def exact_hierarchy_certificate(matrix: Sequence[Sequence[object]], r: int,
                                lower: object = 0, upper: object = 1) -> dict[str, object]:
    """Check lo*I <= A <= hi*I and the rational hierarchy with width hi-lo.

    A spectral interval is a certified upper bound for the true diameter.
    The conclusion uses that interval width, not a computed exact eigenvalue.
    """
    a = [[_fraction(x) for x in row] for row in matrix]
    n = len(a)
    if not isinstance(r, int) or isinstance(r, bool) or r < 2 or n < 4*r-2:
        raise ValueError("The hierarchy requires integer r>=2 and n>=4r-2")
    if n > 14 or any(len(row) != n for row in a):
        raise ValueError("Exact certificates are capped at square order 14")
    if any(a[i][j] != a[j][i] for i in range(n) for j in range(n)):
        raise ValueError("Expected exact symmetry")
    lo, hi = _fraction(lower), _fraction(upper)
    if lo > hi:
        raise ValueError("lower must not exceed upper")
    low_slack = [[a[i][j] - (lo if i == j else 0) for j in range(n)] for i in range(n)]
    high_slack = [[(hi if i == j else 0)-a[i][j] for j in range(n)] for i in range(n)]
    if not _psd_exact(low_slack) or not _psd_exact(high_slack):
        raise ValueError("Exact spectral interval hypotheses failed")
    values = exact_principal_hafnians(a)
    h_low = sum((v*v for mask,v in values.items() if mask.bit_count()==2*r-2), Fraction(0))
    h_high = sum((v*v for mask,v in values.items() if mask.bit_count()==2*r), Fraction(0))
    average_low, average_high = h_low/comb(n,2*r-2), h_high/comb(n,2*r)
    right = (hi-lo)**2 * average_low/4
    defect = right-average_high
    if defect < 0:
        raise ArithmeticError("Exact hierarchy conclusion failed")
    return {"dimension": n, "r": r, "lower": lo, "upper": hi,
            "spectral_interval_verified": True, "H_lower": h_low, "H_upper": h_high,
            "average_lower": average_low, "average_upper": average_high,
            "upper_bound": right, "defect": defect, "equality": defect == 0}
