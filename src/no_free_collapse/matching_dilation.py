"""Complementary-subset dilation for weighted Kneser matching operators.

For real symmetric A of order n=2m, the operator on (m-1)-subsets has norm
at most (m+1)/4 times the spectral diameter of A. The analytic theorem has
no dimension cap. These explicitly enumerated diagnostic routines cap n at
10 to prevent accidental exponential allocations. Exact PSD/gradient instance
certification is provided for the six-variable problem.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from typing import Sequence

import numpy as np

from .zero_diagonal import _fraction, _psd_exact


def _order(n: int) -> None:
    if not isinstance(n, int) or isinstance(n, bool) or n < 4 or n % 2:
        raise ValueError("Expected even order at least four")
    if n > 10:
        raise ValueError("Explicit enumeration is capped at order ten; the theorem is not")


def _real_symmetric(matrix: np.ndarray) -> np.ndarray:
    raw = np.asarray(matrix)
    if np.iscomplexobj(raw):
        raise ValueError("The matching-dilation theorem here uses real symmetric matrices")
    a = np.asarray(raw, dtype=float)
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("Expected a square matrix")
    _order(len(a))
    if not np.all(np.isfinite(a)):
        raise ValueError("Matrix entries must be finite")
    if not np.allclose(a, a.T, rtol=0, atol=1e-12):
        raise ValueError("Expected a symmetric matrix")
    return (a + a.T) / 2


@lru_cache(maxsize=4)
def _incidence(n: int) -> tuple[tuple[tuple[int, ...], ...], np.ndarray, np.ndarray]:
    _order(n)
    m = n // 2
    features = tuple(combinations(range(n), m - 1))
    index = {s: j for j, s in enumerate(features)}
    # Exactly one representative of each complementary pair contains vertex 0.
    halves = [r for r in combinations(range(n), m) if 0 in r]
    columns = np.empty((len(halves), n), dtype=np.int64)
    signs = np.empty_like(columns)
    universe = set(range(n))
    for b, r in enumerate(halves):
        rset = set(r)
        complement = universe - rset
        for i in range(n):
            containing = rset if i in rset else complement
            columns[b, i] = index[tuple(sorted(containing - {i}))]
            signs[b, i] = 1 if i in rset else -1
    columns.flags.writeable = False
    signs.flags.writeable = False
    return features, columns, signs


def dilation_incidence(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Return integer sparse row-to-column indices and signs for R_+,R_-.

    Both maps have one nonzero entry per row. R_+ has coefficient +1;
    R_- uses the returned signs. They satisfy R_+^T R_+=R_-^T R_-=(m+1)I.
    Copies are returned to protect the internal immutable cache.
    """
    _, columns, signs = _incidence(n)
    return columns.copy(), signs.copy()


def matching_operator(matrix: np.ndarray) -> np.ndarray:
    """Build T_m(A) directly from disjoint (m-1)-subsets.

    The entry for disjoint E,F is A[i,j], where {i,j} is their complement;
    intersecting subsets have entry zero. Diagonal entries of A are ignored.
    """
    a = _real_symmetric(matrix)
    features, _, _ = _incidence(len(a))
    sets = [set(s) for s in features]
    result = np.zeros((len(features), len(features)))
    universe = set(range(len(a)))
    for e, f in combinations(range(len(features)), 2):
        if sets[e].isdisjoint(sets[f]):
            i, j = sorted(universe - sets[e] - sets[f])
            result[e, f] = result[f, e] = a[i, j]
    return result


def _compress(a: np.ndarray, columns: np.ndarray, signs: np.ndarray, d: int) -> np.ndarray:
    result = np.zeros((d, d), dtype=a.dtype)
    for cols, sg in zip(columns, signs):
        np.add.at(result, (cols[:, None], cols[None, :]), a * sg[:, None] * sg[None, :])
    return result


def dilation_operator(matrix: np.ndarray) -> np.ndarray:
    """Compute T_m(A) through the independent complementary-subset formula."""
    a = _real_symmetric(matrix)
    features, cols, signs = _incidence(len(a))
    plus = _compress(a, cols, np.ones_like(signs), len(features))
    minus = _compress(a, cols, signs, len(features))
    return (plus - minus) / 4


def dilation_psd_slacks(matrix: np.ndarray, bound: float) -> tuple[np.ndarray, np.ndarray]:
    """Return the two PSD decompositions c*L*I-T and c*L*I+T.

    Here c=(m+1)/4. A>=0 and L*I-A>=0 are checked numerically. This is a
    floating-point diagnostic, not a rigorous enclosure; use the exact
    six-variable certificate interface for rational input.
    """
    a = _real_symmetric(matrix)
    ell = float(bound)
    if not np.isfinite(ell) or ell < 0:
        raise ValueError("bound must be finite and nonnegative")
    eigenvalues = np.linalg.eigvalsh(a)
    tol = 1e-11 * max(1.0, ell, np.max(np.abs(eigenvalues)))
    if eigenvalues[0] < -tol or eigenvalues[-1] > ell + tol:
        raise ValueError("Expected 0 <= A <= bound*I")
    features, cols, signs = _incidence(len(a))
    ones = np.ones_like(signs)
    complement = ell * np.eye(len(a)) - a
    minus_slack = (_compress(complement, cols, ones, len(features))
                   + _compress(a, cols, signs, len(features))) / 4
    plus_slack = (_compress(complement, cols, signs, len(features))
                  + _compress(a, cols, ones, len(features))) / 4
    return minus_slack, plus_slack


def _compress_exact(a: list[list[Fraction]], columns: np.ndarray,
                    signs: np.ndarray, d: int) -> list[list[Fraction]]:
    result = [[Fraction(0) for _ in range(d)] for _ in range(d)]
    for cols, sg in zip(columns, signs):
        for i in range(len(a)):
            for j in range(len(a)):
                result[int(cols[i])][int(cols[j])] += a[i][j] * int(sg[i] * sg[j])
    return result


def exact_six_psd_certificate(matrix: Sequence[Sequence[object]],
                              bound: object = 1) -> dict[str, object]:
    """Certify a rational real PSD instance, its dilation and gradient inequality.

    Integers and Fraction values only. Both input spectral hypotheses and both
    operator slacks are checked by exact PSD elimination. The two independent
    operator constructions must agree entry by entry over the rationals.
    """
    a = [[_fraction(x) for x in row] for row in matrix]
    if len(a) != 6 or any(len(row) != 6 for row in a):
        raise ValueError("Exact gradient certificates require order six")
    if any(a[i][j] != a[j][i] for i in range(6) for j in range(6)):
        raise ValueError("Expected exact symmetry")
    ell = _fraction(bound)
    if ell < 0:
        raise ValueError("bound must be nonnegative")
    complement = [[(ell if i == j else Fraction(0)) - a[i][j]
                   for j in range(6)] for i in range(6)]
    if not _psd_exact(a) or not _psd_exact(complement):
        raise ValueError("Exact spectral hypotheses 0 <= A <= bound*I failed")
    features, cols, signs = _incidence(6)
    ones = np.ones_like(signs)
    d = len(features)
    direct = [[Fraction(0) for _ in range(d)] for _ in range(d)]
    for e, f in combinations(range(d), 2):
        if set(features[e]).isdisjoint(features[f]):
            i, j = sorted(set(range(6)) - set(features[e]) - set(features[f]))
            direct[e][f] = direct[f][e] = a[i][j]
    pa = _compress_exact(a, cols, ones, d)
    ma = _compress_exact(a, cols, signs, d)
    pc = _compress_exact(complement, cols, ones, d)
    mc = _compress_exact(complement, cols, signs, d)
    minus = [[(pc[i][j] + ma[i][j]) / 4 for j in range(d)] for i in range(d)]
    plus = [[(mc[i][j] + pa[i][j]) / 4 for j in range(d)] for i in range(d)]
    for i in range(d):
        for j in range(d):
            eye = ell if i == j else Fraction(0)
            if ((pa[i][j] - ma[i][j]) / 4 != direct[i][j]
                    or minus[i][j] != eye - direct[i][j]
                    or plus[i][j] != eye + direct[i][j]):
                raise ArithmeticError("Exact dilation identity failed")
    if not _psd_exact(minus) or not _psd_exact(plus):
        raise ArithmeticError("Exact operator slack verification failed")
    q1 = sum(a[i][j]**2 for i, j in combinations(range(6), 2))
    q2 = Fraction(0)
    for i, j, k, l in combinations(range(6), 4):
        q2 += (a[i][j]*a[k][l] + a[i][k]*a[j][l] + a[i][l]*a[j][k])**2
    haf = Fraction(0)
    for j in range(1, 6):
        k, l, m, n = [i for i in range(1, 6) if i != j]
        haf += a[0][j]*(a[k][l]*a[m][n] + a[k][m]*a[l][n] + a[k][n]*a[l][m])
    gradient_defect = ell**2*q1/4-q2
    hafnian_defect = ell*q1/6-abs(haf)
    if gradient_defect < 0 or hafnian_defect < 0:
        raise ArithmeticError("Exact corollary verification failed")
    return {
        "dimension": 6, "bound": ell, "input_psd_verified": True,
        "upper_spectral_bound_verified": True, "dilation_residual": Fraction(0),
        "operator_minus_slack_psd": True, "operator_plus_slack_psd": True,
        "q1": q1, "q2": q2, "gradient_defect": gradient_defect,
        "hafnian": haf, "hafnian_energy_defect": hafnian_defect,
    }
