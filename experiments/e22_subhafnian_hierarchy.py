"""Reproduce the all-orders spectral subhafnian hierarchy and crossover tests."""
from __future__ import annotations

import argparse
import json
from math import comb
from pathlib import Path
import platform

import numpy as np

from no_free_collapse.subhafnian_hierarchy import (
    all_principal_hafnians,
    averaged_step_bound,
    averaged_subhafnian_energy,
    edge_density_upper_bound,
    iterated_hybrid_bound,
    roos_G,
    spectral_diameter,
)


def matching_projection(n: int) -> np.ndarray:
    if n % 2:
        raise ValueError("matching projection requires even order")
    a = np.eye(n) / 2
    for i in range(0, n, 2):
        a[i, i + 1] = a[i + 1, i] = 0.5 * (-1 if i % 4 else 1)
    return a


def balanced_rank_one(n: int) -> np.ndarray:
    s = np.where(np.arange(n) % 3 == 0, -1.0, 1.0)
    return np.outer(s, s) / n


def run(seed: int, samples: int) -> dict:
    if seed < 0 or samples < 1:
        raise ValueError("seed must be nonnegative and samples positive")
    rng = np.random.default_rng(seed)
    cases = [(6, 2), (8, 2), (10, 2), (10, 3), (12, 3), (14, 3), (14, 4)]
    random_checks = []
    max_ratio = 0.0
    for n, r in cases:
        local_max = 0.0
        for _ in range(samples):
            a = rng.normal(size=(n, n)); a = (a + a.T) / 2
            h = all_principal_hafnians(a)
            left, right = averaged_step_bound(a, r, hafnians=h)
            ratio = 0.0 if right == 0 else left / right
            if ratio > 1 + 1e-10:
                raise ArithmeticError("spectral hierarchy regression failed")
            local_max = max(local_max, ratio)
            max_ratio = max(max_ratio, ratio)
        random_checks.append({"n": n, "r": r, "samples": samples,
                              "largest_left_over_bound": local_max})

    critical = []
    for r in range(2, 7):
        n = 4 * r - 2
        spectral = 1 / (4 * (2 * r - 1) ** 2)
        rank_one_g1 = 1 / n**2
        matching_g1 = 1 / (4 * (4 * r - 3))
        critical.append({
            "r": r, "critical_order": n,
            "rank_one_top_step_roos_over_spectral": rank_one_g1 / spectral,
            "matching_top_step_roos_over_spectral": matching_g1 / spectral,
            "closed_form_matching_factor": (2 * r - 1) ** 2 / (4 * r - 3),
        })

    # Exact numerical equality checks are kept at sizes where the 2^n diagnostic
    # remains lightweight; larger rows in the table use closed forms above.
    equality = []
    for r in (2, 3, 4):
        n = 4 * r - 2
        for family, a in (("balanced_rank_one", balanced_rank_one(n)),
                          ("matching_projection", matching_projection(n))):
            h = all_principal_hafnians(a)
            left, right = averaged_step_bound(a, r, hafnians=h)
            equality.append({"r": r, "n": n, "family": family,
                             "left": left, "right": right,
                             "absolute_residual": abs(left - right)})

    hybrid = []
    for n, r in ((10, 3), (12, 3), (14, 4)):
        for _ in range(max(2, samples // 10)):
            q = np.linalg.qr(rng.normal(size=(n, n)))[0]
            vals = rng.uniform(0, 1, size=n)
            a = (q * vals) @ q.T
            h = all_principal_hafnians(a)
            actual, upper = iterated_hybrid_bound(a, r, hafnians=h)
            g1 = roos_G(a, 1, hafnians=h)
            d = spectral_diameter(a)
            pure_roos = g1**r
            pure_spectral = g1 * np.prod([d*d/(4*(2*s-1)**2)
                                          for s in range(2, r + 1)])
            if actual > upper * (1 + 1e-10):
                raise ArithmeticError("hybrid regression failed")
            hybrid.append({"n": n, "r": r, "actual_G_r": actual,
                           "hybrid_bound": upper, "pure_roos_bound": pure_roos,
                           "pure_spectral_bound": float(pure_spectral)})

    # Real-interferometer GBS diagnostics: B=O diag(tanh(s_j)) O^T is real PSD
    # and has ||B||<1.  The common collision-free normalization cancels, so the
    # averaged per-pattern probability ratio is the same as the Hbar ratio.
    gbs = []
    for n, r in ((10, 3), (14, 4)):
        for _ in range(max(2, samples // 10)):
            q = np.linalg.qr(rng.normal(size=(n, n)))[0]
            squeezings = rng.uniform(0.0, 1.4, size=n)
            vals = np.tanh(squeezings)
            b = (q * vals) @ q.T
            h = all_principal_hafnians(b)
            current = averaged_subhafnian_energy(b, r, hafnians=h)
            previous = averaged_subhafnian_energy(b, r - 1, hafnians=h)
            d = spectral_diameter(b)
            ratio = 0.0 if previous == 0 else current / previous
            certified = d*d/4
            if ratio > certified * (1 + 1e-10):
                raise ArithmeticError("GBS-sector hierarchy regression failed")
            gbs.append({"modes": n, "pair_order": r,
                        "average_weight_ratio": ratio,
                        "spectral_certificate": certified,
                        "operator_norm": float(np.linalg.norm(b, 2))})

    return {
        "seed": seed, "samples_per_random_case": samples,
        "python": platform.python_version(), "numpy": np.__version__,
        "random_hierarchy_checks": random_checks,
        "largest_random_left_over_bound": max_ratio,
        "critical_sharpness_and_roos_crossover": critical,
        "explicit_critical_equalities": equality,
        "hybrid_psd_checks": hybrid,
        "real_gbs_like_checks": gbs,
        "universal_edge_density_bound": "theta <= 1/(n-1)",
        "analytic_proof": "docs/subhafnian_hierarchy.md",
        "cube_normalized_1_over_216_conjecture_solved": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20260912)
    parser.add_argument("--samples", type=int, default=20)
    parser.add_argument("--output", type=Path,
                        default=Path("results/subhafnian_hierarchy.json"))
    args = parser.parse_args()
    if args.seed < 0 or args.samples < 1:
        parser.error("seed must be nonnegative and samples positive")
    result = run(args.seed, args.samples)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(args.output),
        "random_cases": len(result["random_hierarchy_checks"]),
        "largest_random_left_over_bound": result["largest_random_left_over_bound"],
        "critical_rows": len(result["critical_sharpness_and_roos_crossover"]),
        "gbs_checks": len(result["real_gbs_like_checks"]),
    }, indent=2))


if __name__ == "__main__":
    main()
