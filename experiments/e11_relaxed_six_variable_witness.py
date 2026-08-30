"""Deterministic six-variable relaxation witness separating range/degree from PSD."""

from __future__ import annotations

import itertools
import math

import numpy as np


def witness_coefficients() -> np.ndarray:
    """Return pair coefficients for a bounded degree-2 non-PSD witness.

    The polynomial is q(x)=1/2+sum_{i<j} c_ij x_i x_j with c_ij=+/-1/10.
    It maps the six-dimensional Boolean cube into [0,1] and has
    full-parity coefficient of q^3 equal to 3/100, hence the normalized
    hafnian-style objective is 1/200. Its equal-diagonal quadratic matrix is
    not PSD, showing that PSD geometry is essential beyond range/degree.
    """
    signs = [
        1, 1, -1, -1, 1,
        1, 1, -1, -1,
        -1, 1, -1,
        -1, -1,
        -1,
    ]
    return np.asarray(signs, dtype=np.float64) / 10.0


def evaluate_witness() -> tuple[float, float, float, float]:
    coeffs = witness_coefficients()
    xs = list(itertools.product((-1.0, 1.0), repeat=6))
    values = []
    parity_moment = 0.0
    for x in xs:
        q = 0.5
        k = 0
        for i in range(6):
            for j in range(i + 1, 6):
                q += coeffs[k] * x[i] * x[j]
                k += 1
        values.append(q)
        parity_moment += math.prod(x) * q**3
    parity_moment /= len(xs)

    B = np.eye(6) / 12.0
    k = 0
    for i in range(6):
        for j in range(i + 1, 6):
            B[i, j] = B[j, i] = coeffs[k] / 2.0
            k += 1
    return min(values), max(values), parity_moment / math.factorial(3), float(np.linalg.eigvalsh(B)[0])


def main() -> None:
    qmin, qmax, objective, min_eigenvalue = evaluate_witness()
    print(f"q range=[{qmin:.12f}, {qmax:.12f}]")
    print(f"relaxed objective={objective:.12f} (expected 1/200={1/200:.12f})")
    print(f"equal-diagonal quadratic min eigenvalue={min_eigenvalue:.12f}")
    print(f"PSD disjoint-pair value 1/216={1/216:.12f}")


if __name__ == "__main__":
    main()
