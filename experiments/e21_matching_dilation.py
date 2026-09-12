"""Reproduce exact matrix-basis and rational PSD dilation checks."""
from __future__ import annotations
import argparse
from fractions import Fraction
import json
from pathlib import Path
import platform
import numpy as np
from no_free_collapse.matching_dilation import (
    _compress, dilation_incidence, dilation_operator,
    exact_six_psd_certificate, matching_operator,
)


def run(seed: int, samples: int) -> dict:
    if seed < 0 or samples < 1:
        raise ValueError("seed must be nonnegative and samples positive")
    rng = np.random.default_rng(seed)
    basis_count = 0
    maximum_residual = 0.0
    minimum_slack = float("inf")
    for n in (4, 6, 8, 10):
        cols, signs = dilation_incidence(n)
        d = int(cols.max()) + 1
        counts = np.bincount(cols.ravel(), minlength=d)
        if not np.array_equal(counts, np.full(d, n//2+1)):
            raise ArithmeticError("Incidence isometry check failed")
        for i in range(n):
            for j in range(i, n):
                a = np.zeros((n, n), dtype=np.int64)
                a[i, j] = a[j, i] = 1
                difference = (_compress(a, cols, np.ones_like(signs), d)
                              - _compress(a, cols, signs, d))
                if not np.array_equal(difference, 4*matching_operator(a).astype(np.int64)):
                    raise ArithmeticError("Exact basis identity failed")
                basis_count += 1
        for _ in range(samples):
            a = rng.normal(size=(n, n)); a = (a+a.T)/2
            es = np.linalg.eigvalsh(a)
            direct = matching_operator(a)
            residual = float(np.max(np.abs(direct-dilation_operator(a))))
            slack = (n//2+1)*(es[-1]-es[0])/4-np.linalg.norm(direct, 2)
            maximum_residual = max(maximum_residual, residual)
            minimum_slack = min(minimum_slack, float(slack))
            if residual > 1e-10 or slack < -1e-10:
                raise ArithmeticError("Numerical dilation regression failed")
    rational = []
    for _ in range(20):
        z = rng.integers(-3, 4, size=(6, 6))
        q = z@z.T
        bound = int(np.trace(q))
        a = [[Fraction(int(x), bound) for x in row] for row in q]
        rational.append(exact_six_psd_certificate(a, 1))
    return {
        "seed": seed, "samples_per_order": samples,
        "python": platform.python_version(), "numpy": np.__version__,
        "integer_symmetric_basis_matrices_checked": basis_count,
        "random_real_symmetric_matrices_checked": 4*samples,
        "maximum_numerical_dilation_residual": maximum_residual,
        "minimum_numerical_spectral_bound_slack": minimum_slack,
        "exact_rational_psd_certificates": rational,
        "universal_proof": "docs/matching_dilation.md",
        "proved_scope": "sharp all-even-dimensional real spectral operator bound; full six-variable PSD gradient bound",
        "cube_normalized_1_over_216_conjecture_solved": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20260912)
    parser.add_argument("--samples", type=int, default=100)
    parser.add_argument("--output", type=Path, default=Path("results/matching_dilation.json"))
    args = parser.parse_args()
    if args.seed < 0 or args.samples < 1:
        parser.error("seed must be nonnegative and samples positive")
    result = run(args.seed, args.samples)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, default=str)+"\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k != "exact_rational_psd_certificates"}, indent=2))


if __name__ == "__main__":
    main()
