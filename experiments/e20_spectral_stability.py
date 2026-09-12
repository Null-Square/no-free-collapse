"""Reproduce spectral-stability checks using only NumPy and standard Python."""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import platform
from time import perf_counter

import numpy as np

from no_free_collapse.zero_diagonal import (
    defect_terms, exact_certificate, recover_six_hafnian_matching,
    six_hafnian, subhafnian_energy_direct, subhafnian_energy_trace,
)


def run(seed: int, samples: int) -> dict:
    """Return reproducibility data; raise on a violated numerical regression."""
    if samples < 1 or seed < 0:
        raise ValueError("samples must be positive and seed must be nonnegative")
    rng = np.random.default_rng(seed)
    max_identity_error = 0.0
    min_defect = float("inf")
    count = 0
    for n in range(4, 13):
        for complex_entries in (False, True):
            for _ in range(samples):
                a = rng.normal(size=(n, n))
                if complex_entries:
                    a = a + 1j * rng.normal(size=(n, n))
                a = (a + a.T) / 2
                np.fill_diagonal(a, 0)
                a /= np.linalg.norm(a, 2)
                terms = defect_terms(a, 1)
                error = abs(terms.defect - terms.term_sum)
                max_identity_error = max(max_identity_error, error)
                min_defect = min(min_defect, terms.defect)
                count += 1
                if terms.defect < -1e-10 or error > 1e-10:
                    raise ArithmeticError("Random regression failure")

    exact = []
    for n in (4, 6, 8, 10, 12):
        a = [[Fraction(0) for _ in range(n)] for _ in range(n)]
        for i in range(0, n, 2):
            a[i][i + 1] = a[i + 1][i] = Fraction(1)
        exact.append(exact_certificate(a, 1))

    stability = []
    for complex_entries in (False, True):
        m = np.zeros((6, 6), complex if complex_entries else float)
        for i in range(0, 6, 2):
            phase = np.exp(1j * (i + 1)) if complex_entries else (-1.0)**(i // 2)
            m[i, i + 1] = m[i + 1, i] = phase
        for strength in (0, 1e-5, 1e-4, 1e-3):
            noise = rng.normal(size=(6, 6))
            if complex_entries:
                noise = noise + 1j * rng.normal(size=(6, 6))
            noise = (noise + noise.T) / 2
            np.fill_diagonal(noise, 0)
            a = m + strength * noise
            a /= max(1.0, np.linalg.norm(a, 2))
            recovered, upper = recover_six_hafnian_matching(a)
            distance = float(np.linalg.norm(a - recovered, "fro")**2)
            if distance > upper + 1e-10:
                raise ArithmeticError("Stability regression failure")
            stability.append({
                "complex": complex_entries, "perturbation": strength,
                "hafnian_modulus": float(abs(six_hafnian(a))),
                "squared_distance": distance, "proved_upper_bound": upper,
            })

    # Single-run timing diagnostics, not a controlled performance benchmark.
    timings = []
    for n in (16, 32, 64):
        a = rng.normal(size=(n, n))
        a = (a + a.T) / 2
        np.fill_diagonal(a, 0)
        start = perf_counter()
        direct = subhafnian_energy_direct(a)
        direct_seconds = perf_counter() - start
        start = perf_counter()
        trace = subhafnian_energy_trace(a)
        trace_seconds = perf_counter() - start
        if not np.isclose(direct, trace, rtol=1e-10):
            raise ArithmeticError("Trace diagnostic disagreement")
        timings.append({
            "n": n, "direct_seconds": direct_seconds,
            "trace_seconds": trace_seconds,
            "relative_difference": abs(direct - trace) / max(1.0, abs(direct)),
        })

    return {
        "seed": seed, "samples_per_dimension_and_field": samples,
        "python": platform.python_version(), "numpy": np.__version__,
        "random_matrices_checked": count,
        "max_identity_absolute_error": max_identity_error,
        "minimum_random_defect": min_defect,
        "exact_matching_certificates": exact,
        "six_hafnian_stability": stability, "timing_diagnostics": timings,
        "proof_scope": "zero-diagonal complex symmetric spectral constraint",
        "full_real_psd_gradient_proof_is_separate": "docs/matching_dilation.md",
        "cube_normalized_1_over_216_conjecture_solved": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20260912)
    parser.add_argument("--samples", type=int, default=250)
    parser.add_argument("--output", type=Path, default=Path("results/spectral_stability.json"))
    args = parser.parse_args()
    if args.samples < 1 or args.seed < 0:
        parser.error("--samples must be positive; --seed must be nonnegative")
    result = run(args.seed, args.samples)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, default=str) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(args.output),
        "random_matrices_checked": result["random_matrices_checked"],
        "max_identity_absolute_error": result["max_identity_absolute_error"],
    }, indent=2))


if __name__ == "__main__":
    main()
