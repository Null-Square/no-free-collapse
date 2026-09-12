"""Sharp four-subhafnian energy bounds for zero-diagonal symmetric matrices.

The analytic proof is in docs/zero_diagonal_spectral.md.  Floating-point
routines are diagnostics.  ``exact_certificate`` checks rational input with
exact arithmetic, including both semidefinite spectral-bound conditions.
No routine certifies the unrestricted PSD-contraction conjecture.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from numbers import Integral
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class DefectTerms:
    """Four nonnegative terms when ``bound >= ||matrix||_op``."""
    dimension: int
    bound: float
    edge_energy: float
    subhafnian_energy: float
    occupancy: float
    row_variance: float
    wedge_spreading: float
    spectral_slack: float

    @property
    def defect(self) -> float:
        return (self.dimension - 2) * self.bound**2 * self.edge_energy / 4 - self.subhafnian_energy

    @property
    def term_sum(self) -> float:
        return self.occupancy + self.row_variance + self.wedge_spreading + self.spectral_slack


def _matrix(matrix: np.ndarray, *, atol: float = 1e-12) -> np.ndarray:
    raw = np.asarray(matrix)
    a = np.asarray(raw, dtype=complex if np.iscomplexobj(raw) else float)
    if a.ndim != 2 or a.shape[0] != a.shape[1] or a.shape[0] < 4:
        raise ValueError("Expected a square matrix of order at least four.")
    if not np.all(np.isfinite(a)):
        raise ValueError("Matrix entries must be finite.")
    if not np.isfinite(atol) or atol < 0:
        raise ValueError("atol must be finite and nonnegative.")
    if not np.allclose(a, a.T, rtol=0, atol=atol):
        raise ValueError("Expected a symmetric matrix.")
    if np.any(np.abs(np.diag(a)) > atol):
        raise ValueError("Expected a zero diagonal; removing a diagonal can change the norm.")
    # Explicitly make the accepted floating-point input symmetric and zero-diagonal.
    a = (a + a.T) / 2
    np.fill_diagonal(a, 0.0)
    return a


def subhafnian_energy_direct(matrix: np.ndarray) -> float:
    """Sum haf(B[S])**2 over four-subsets S, by direct O(n**4) enumeration."""
    a = _matrix(matrix)
    return float(sum(abs(a[i,j]*a[k,l] + a[i,k]*a[j,l] + a[i,l]*a[j,k])**2
                     for i,j,k,l in combinations(range(len(a)), 4)))


def subhafnian_energy_trace(matrix: np.ndarray) -> float:
    """O(n**3) trace formula, without enumerating four-subsets.

    Cancellation is possible near zero.  Use direct enumeration or an exact
    certificate when a rigorous sign, rather than an approximate value, matters.
    """
    a = _matrix(matrix)
    squared = np.abs(a)**2
    t = float(squared.sum()/2)
    rows = squared.sum(axis=1)
    a2 = a @ a.conj().T
    return float(t*t/2 + (squared*squared).sum()/2 - rows@rows + np.square(np.abs(a2)).sum()/4)


def defect_terms(matrix: np.ndarray, bound: float | None = None) -> DefectTerms:
    """Evaluate the exact identity in floating point after checking its hypotheses."""
    a = _matrix(matrix)
    norm = float(np.linalg.norm(a, ord=2))
    ell = norm if bound is None else float(bound)
    if not np.isfinite(ell) or ell < 0:
        raise ValueError("The spectral bound must be finite and nonnegative.")
    if norm > ell + 1e-11*max(1.0, norm, ell):
        raise ValueError("The supplied bound is smaller than the numerical operator norm.")
    n = len(a)
    weights = np.abs(a)**2
    rows = weights.sum(axis=1)
    t = float(rows.sum()/2)
    # Sum products directly to avoid cancellation in wedge-spreading energy.
    wedge = sum(weights[i,j]*weights[i,k] for i in range(n)
                for j,k in combinations([v for v in range(n) if v != i], 2))
    a2 = a @ a.conj().T
    return DefectTerms(n, ell, t, subhafnian_energy_direct(a),
                       (n-4)/(2*n)*t*(n*ell*ell/2-t),
                       float(np.square(rows-2*t/n).sum()/2),
                       float(wedge),
                       float((ell*ell*np.trace(a2).real-np.sum(np.abs(a2)**2))/4))


def recover_matching(matrix: np.ndarray) -> tuple[np.ndarray, float]:
    """Recover a unit-modulus matching under the stated quantitative hypotheses.

    Input must be a contraction of order n>=6 with T>=n/4 and defect
    delta <= min((n-4)/64, 1/16).  Returns (matching, proved squared-distance
    upper bound).  Floating-point validation does not replace exact certification.
    """
    a = _matrix(matrix)
    n = len(a)
    if n < 6 or n % 2:
        raise ValueError("Matching recovery requires even order at least six.")
    terms = defect_terms(a, 1.0)
    delta = terms.defect
    if terms.edge_energy < n/4 - 1e-12:
        raise ValueError("The nondegeneracy condition T >= n/4 is required.")
    if delta < -1e-10 or delta > min((n-4)/64, 1/16) + 1e-12:
        raise ValueError("Defect lies outside the proved stability regime.")
    partners = np.argmax(np.abs(a)**2, axis=1)
    matching = np.zeros_like(a)
    for i,j in enumerate(partners):
        if i == j or partners[j] != i or abs(a[i,j])**2 <= .5:
            raise ValueError("Numerical input does not yield the certified mutual matching.")
        matching[i,j] = a[i,j]/abs(a[i,j])
    return matching, (16/(n-4)+16/3)*max(0.0, delta)


def _fraction(x: object) -> Fraction:
    if isinstance(x, Fraction):
        return x
    if isinstance(x, Integral) and not isinstance(x, bool):
        return Fraction(int(x))
    raise TypeError("Exact certificates accept integers and Fraction values, not floats.")


def _psd_exact(a: list[list[Fraction]]) -> bool:
    """Exact symmetric elimination with singular-pivot handling."""
    work = [row[:] for row in a]
    n = len(work)
    for k in range(n):
        pivot = work[k][k]
        if pivot < 0:
            return False
        if pivot == 0:
            if any(work[k][j] != 0 for j in range(k+1,n)):
                return False
            continue
        for i in range(k+1,n):
            for j in range(i,n):
                work[i][j] -= work[i][k]*work[k][j]/pivot
                work[j][i] = work[i][j]
    return True


def exact_certificate(matrix: Sequence[Sequence[object]], bound: object = 1) -> dict[str, object]:
    """Verify the inequality and four-term identity over the rationals.

    Rejects out-of-domain inputs rather than silently projecting them into the
    theorem's domain.  A returned certificate has exact nonnegative terms and
    zero residual; this verifies the supplied instance, not a universal theorem.
    """
    a = [[_fraction(x) for x in row] for row in matrix]
    n = len(a)
    if n < 4 or any(len(row) != n for row in a):
        raise ValueError("Expected a square matrix of order at least four.")
    if any(a[i][j] != a[j][i] for i in range(n) for j in range(n)):
        raise ValueError("Expected an exactly symmetric matrix.")
    if any(a[i][i] for i in range(n)):
        raise ValueError("Expected an exactly zero diagonal.")
    ell = _fraction(bound)
    if ell < 0:
        raise ValueError("The spectral bound must be nonnegative.")
    for sign in (-1,1):
        test = [[(ell if i == j else Fraction(0))+sign*a[i][j]
                 for j in range(n)] for i in range(n)]
        if not _psd_exact(test):
            raise ValueError("The supplied exact spectral bound is invalid.")
    weights = [[x*x for x in row] for row in a]
    rows = [sum(row) for row in weights]
    t = sum(rows)/2
    h = sum((a[i][j]*a[k][l]+a[i][k]*a[j][l]+a[i][l]*a[j][k])**2
            for i,j,k,l in combinations(range(n),4))
    a2 = [[sum(a[i][k]*a[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    tr4 = sum(x*x for row in a2 for x in row)
    terms = {
        "occupancy": Fraction(n-4,2*n)*t*(n*ell*ell/2-t),
        "row_variance": sum((r-2*t/n)**2 for r in rows)/2,
        "wedge_spreading": sum(weights[i][j]*weights[i][k] for i in range(n)
                               for j,k in combinations([v for v in range(n) if v!=i],2)),
        "spectral_slack": (2*ell*ell*t-tr4)/4,
    }
    delta = Fraction(n-2,4)*ell*ell*t-h
    residual = delta-sum(terms.values())
    if residual != 0 or any(v < 0 for v in terms.values()):
        raise ArithmeticError("Exact defect verification failed.")
    return {"dimension": n, "bound": ell, "edge_energy": t,
            "subhafnian_energy": h, "defect": delta, "terms": terms,
            "identity_residual": residual, "spectral_bound_verified": True}


def six_hafnian(matrix: np.ndarray) -> complex | float:
    """Exact polynomial evaluation (floating point) for order six only."""
    a = _matrix(matrix)
    if len(a) != 6:
        raise ValueError("Expected order six.")
    total = 0j if np.iscomplexobj(a) else 0.0
    for j in range(1,6):
        k,l,m,n = [i for i in range(1,6) if i != j]
        total += a[0,j]*(a[k,l]*a[m,n]+a[k,m]*a[l,n]+a[k,n]*a[l,m])
    return total


def recover_six_hafnian_matching(matrix: np.ndarray) -> tuple[np.ndarray, float]:
    """If ||B||<=1 and |haf(B)|>=1-1/48, recover M with ||B-M||_F^2<=34 eps.

    Here eps=1-|haf(B)|; phases, rather than only real signs, are allowed.
    This is a numerical implementation of the analytic stability theorem.
    """
    a = _matrix(matrix)
    if len(a) != 6:
        raise ValueError("Expected order six.")
    defect_terms(a, 1.0)  # Validates the contraction hypothesis numerically.
    eps = 1.0-abs(six_hafnian(a))
    if eps < -1e-10 or eps > 1/48 + 1e-12:
        raise ValueError("Input lies outside the proved near-extremal regime.")
    partners = np.argmax(np.abs(a)**2, axis=1)
    matching = np.zeros_like(a)
    for i,j in enumerate(partners):
        if i == j or partners[j] != i or abs(a[i,j])**2 <= .5:
            raise ValueError("No certified mutual matching was recovered.")
        matching[i,j] = a[i,j]/abs(a[i,j])
    return matching, 34.0*max(0.0, float(eps))
