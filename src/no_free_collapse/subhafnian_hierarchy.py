"""Averaged subhafnian energies and spectral hierarchy diagnostics.

The universal proof is analytic; see docs/subhafnian_hierarchy.md.  This module
provides floating-point combinatorial evaluation for small matrices and the closed-form
bounds used by the reproducibility tests.  It is intentionally not an optimizer.
"""
from __future__ import annotations

from functools import lru_cache
from collections.abc import Mapping, Iterator
from dataclasses import dataclass
from hashlib import sha256
from types import MappingProxyType
from math import comb, factorial

import numpy as np


def _symmetric_matrix(matrix: np.ndarray, *, atol: float = 1e-12) -> np.ndarray:
    raw = np.asarray(matrix)
    a = np.asarray(raw, dtype=complex if np.iscomplexobj(raw) else float)
    if a.ndim != 2 or a.shape[0] != a.shape[1] or len(a) == 0:
        raise ValueError("Expected a square matrix")
    if not np.all(np.isfinite(a)):
        raise ValueError("Matrix entries must be finite")
    if not np.allclose(a, a.T, rtol=0, atol=atol):
        raise ValueError("Expected a symmetric matrix")
    return a / 2 + a.T / 2


MAX_ENUMERATION_ORDER = 20


def _order(r: int, n: int, minimum: int = 0) -> None:
    if not isinstance(r, int) or isinstance(r, bool) or r < minimum or 2*r > n:
        raise ValueError(f"r must be an integer >= {minimum} with 2r <= matrix order")


def _fingerprint(a: np.ndarray) -> str:
    # Hafnians do not use the diagonal; caches remain valid after diagonal shifts.
    b = a.copy()
    np.fill_diagonal(b, 0)
    return sha256(str((b.shape, b.dtype.str)).encode()+b.tobytes()).hexdigest()


@dataclass(frozen=True)
class PrincipalHafnians(Mapping):
    """Read-only, matrix-bound floating-point cache; not an exact certificate."""
    fingerprint: str
    _entries: Mapping

    def __getitem__(self, key: int):
        return self._entries[key]

    def __iter__(self) -> Iterator[int]:
        return iter(self._entries)

    def __len__(self) -> int:
        return len(self._entries)


def _cache(a: np.ndarray, supplied=None) -> PrincipalHafnians:
    if supplied is None:
        return all_principal_hafnians(a)
    if not isinstance(supplied, PrincipalHafnians) or supplied.fingerprint != _fingerprint(a):
        raise ValueError("Use the unmodified cache returned for this matrix by all_principal_hafnians")
    return supplied


def odd_double_factorial(r: int) -> int:
    """Return (2r-1)!!, with (-1)!! convention replaced by 1 at r=0."""
    if not isinstance(r, int) or isinstance(r, bool) or r < 0:
        raise ValueError("r must be a nonnegative integer")
    return factorial(2 * r) // (factorial(r) * 2**r) if r else 1


def all_principal_hafnians(matrix: np.ndarray) -> PrincipalHafnians:
    """Compute every even principal hafnian by a 2^n dynamic program.

    The mask 0 has hafnian 1. Odd masks are omitted.  This is intended for
    theorem diagnostics and tests (n <= 20), not large-scale hafnian
    computation. The enforced cap prevents accidental exponential allocation.
    Values use floating-point arithmetic, not exact arithmetic.
    """
    a = _symmetric_matrix(matrix)
    n = len(a)
    if n > MAX_ENUMERATION_ORDER:
        raise ValueError("Principal-hafnian enumeration is capped at order 20")

    @lru_cache(maxsize=None)
    def haf(mask: int):
        if mask == 0:
            return 1.0 + 0j if np.iscomplexobj(a) else 1.0
        if mask.bit_count() % 2:
            return 0.0 + 0j if np.iscomplexobj(a) else 0.0
        i_bit = mask & -mask
        i = i_bit.bit_length() - 1
        rest = mask ^ i_bit
        total = 0.0 + 0j if np.iscomplexobj(a) else 0.0
        jmask = rest
        while jmask:
            j_bit = jmask & -jmask
            j = j_bit.bit_length() - 1
            total += a[i, j] * haf(rest ^ j_bit)
            jmask ^= j_bit
        if not np.isfinite(total):
            raise FloatingPointError("Hafnian overflow; rescale the input or use rational certification")
        return total

    out: dict[int, complex | float] = {0: haf(0)}
    for mask in range(1, 1 << n):
        if mask.bit_count() % 2 == 0:
            out[mask] = haf(mask)
    return PrincipalHafnians(_fingerprint(a), MappingProxyType(out))


def subhafnian_energy(matrix: np.ndarray, r: int, *, hafnians: PrincipalHafnians | None = None) -> float:
    """Return H_{2r}=sum_{|S|=2r} |haf(A[S])|^2."""
    a = _symmetric_matrix(matrix)
    n = len(a)
    if not isinstance(r, int) or isinstance(r, bool) or r < 0 or 2 * r > n:
        raise ValueError("Expected 0 <= 2r <= matrix order")
    h = _cache(a, hafnians)
    return float(sum(abs(value) ** 2 for mask, value in h.items() if mask.bit_count() == 2 * r))


def averaged_subhafnian_energy(matrix: np.ndarray, r: int, *, hafnians=None) -> float:
    """Average H_{2r}/binom(n,2r)."""
    a = _symmetric_matrix(matrix)
    return subhafnian_energy(a, r, hafnians=hafnians) / comb(len(a), 2 * r)


def roos_G(matrix: np.ndarray, r: int, *, hafnians=None) -> float:
    """Roos's normalized squared-hafnian average G(A,r)."""
    return averaged_subhafnian_energy(matrix, r, hafnians=hafnians) / odd_double_factorial(r) ** 2


def spectral_diameter(matrix: np.ndarray) -> float:
    """Spectral diameter for a real symmetric matrix."""
    a = _symmetric_matrix(matrix)
    if np.iscomplexobj(a):
        raise ValueError("Spectral-diameter hierarchy is stated here for real symmetric matrices")
    eig = np.linalg.eigvalsh(a)
    diameter = float(eig[-1] - eig[0])
    if not np.isfinite(diameter):
        raise FloatingPointError("Spectral diameter overflow; rescale the matrix")
    return diameter


def averaged_step_bound(matrix: np.ndarray, r: int, *, hafnians=None) -> tuple[float, float]:
    """Return (left,right) in the proved averaged H_{2r} recurrence.

    Requires n >= 4r-2 and r>=2:
      H_{2r}/C(n,2r) <= d(A)^2/4 * H_{2r-2}/C(n,2r-2).
    """
    a = _symmetric_matrix(matrix)
    n = len(a)
    if not isinstance(r, int) or isinstance(r, bool) or r < 2 or n < 4 * r - 2:
        raise ValueError("The proved recurrence requires r>=2 and n>=4r-2")
    h = _cache(a, hafnians)
    left = averaged_subhafnian_energy(a, r, hafnians=h)
    right = spectral_diameter(a) ** 2 / 4 * averaged_subhafnian_energy(a, r - 1, hafnians=h)
    return left, right


def normalized_step_coefficients(matrix: np.ndarray, r: int, *, hafnians=None) -> dict[str, float]:
    """Compare the spectral and Roos one-step coefficients for G(A,r).

    Roos gives G_r <= G_1 G_{r-1}.  The new spectral recurrence gives
    G_r <= d^2/[4(2r-1)^2] G_{r-1} when n>=4r-2.  Their minimum is the
    certified hybrid coefficient.
    """
    a = _symmetric_matrix(matrix)
    _order(r, len(a), 2)
    if len(a) < 4 * r - 2:
        raise ValueError("The spectral comparison requires r>=2 and n>=4r-2")
    h = _cache(a, hafnians)
    d = spectral_diameter(a)
    g1 = roos_G(a, 1, hafnians=h)
    spectral = d * d / (4 * (2 * r - 1) ** 2)
    return {
        "roos": g1,
        "spectral": spectral,
        "hybrid": min(g1, spectral),
        "edge_density": 0.0 if d == 0 else 4 * g1 / (d * d),
        "spectral_wins": spectral < g1,
    }


def iterated_hybrid_bound(matrix: np.ndarray, r: int, *, hafnians=None) -> tuple[float, float]:
    """Return (G_r, hybrid upper bound) by iterating the best one-step factor."""
    a = _symmetric_matrix(matrix)
    _order(r, len(a), 1)
    if r >= 2 and len(a) < 4 * r - 2:
        raise ValueError("The hybrid spectral chain requires n>=4r-2")
    h = _cache(a, hafnians)
    g1 = roos_G(a, 1, hafnians=h)
    d = spectral_diameter(a)
    upper = g1
    for s in range(2, r + 1):
        upper *= min(g1, d * d / (4 * (2 * s - 1) ** 2))
    return roos_G(a, r, hafnians=h), upper


def edge_density_upper_bound(n: int) -> float:
    """Universal scale-free bound 4 G_1/d^2 <= 1/(n-1)."""
    if not isinstance(n, int) or isinstance(n, bool) or n < 2:
        raise ValueError("n must be an integer at least two")
    return 1 / (n - 1)


def hybrid_crossover_order(matrix: np.ndarray, r: int, *, hafnians=None) -> int:
    """Largest initial order whose hybrid step uses Roos rather than spectrum.

    The returned k lies in [1,r].  The hybrid chain uses the Roos factor G_1
    through order k and the spectral factor for every step k+1,...,r.  If the
    two factors tie at a step, the tie is assigned to the Roos side.
    """
    a = _symmetric_matrix(matrix)
    if not isinstance(r, int) or isinstance(r, bool) or r < 1:
        raise ValueError("r must be a positive integer")
    if r >= 2 and len(a) < 4 * r - 2:
        raise ValueError("The hybrid spectral chain requires n>=4r-2")
    h = _cache(a, hafnians)
    g1 = roos_G(a, 1, hafnians=h)
    d2 = spectral_diameter(a) ** 2
    k = 1
    for s in range(2, r + 1):
        if g1 <= d2 / (4 * (2 * s - 1) ** 2):
            k = s
        else:
            break
    return k


def closed_form_hybrid_bound(matrix: np.ndarray, r: int, *, hafnians=None) -> tuple[float, float, int]:
    """Return (G_r, one-switch hybrid bound, crossover order k).

    With k from ``hybrid_crossover_order``, the upper bound is
      G_1^k (d^2/4)^(r-k) [ (2k-1)!! / (2r-1)!! ]^2.
    This is algebraically identical to ``iterated_hybrid_bound``.
    """
    a = _symmetric_matrix(matrix)
    if not isinstance(r, int) or isinstance(r, bool) or r < 1:
        raise ValueError("r must be a positive integer")
    if r >= 2 and len(a) < 4 * r - 2:
        raise ValueError("The hybrid spectral chain requires n>=4r-2")
    h = _cache(a, hafnians)
    actual = roos_G(a, r, hafnians=h)
    g1 = roos_G(a, 1, hafnians=h)
    d2 = spectral_diameter(a) ** 2
    k = hybrid_crossover_order(a, r, hafnians=h)
    if g1 == 0:
        return actual, 0.0, k
    ratio = odd_double_factorial(k) / odd_double_factorial(r)
    upper = g1**k * (d2 / 4) ** (r - k) * ratio**2
    return actual, float(upper), k


def sharp_edge_density_upper_bound(n: int) -> float:
    """Sharp fixed-order ceiling, also for odd n; proof in the submission audit.

    It is 1/(n-1) for even n and (n+1)/n**2 for odd n. Equality requires
    a shifted/scaled constant-diagonal projection of rank floor(n/2) or ceil(n/2).
    """
    edge_density_upper_bound(n)  # Validate n with the public legacy function.
    return 4 * (n*n//4) / (n*n*(n-1))
