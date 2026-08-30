"""Verify conditioning-controlled decay of normalization-induced interactions."""

from __future__ import annotations

import numpy as np

from no_free_collapse.born import normalized_born_probability, polynomial_state
from no_free_collapse.conditioning import (
    centered_relative_variation,
    spectral_coefficient_bound,
    spectral_tail_energy_bound,
)
from no_free_collapse.interactions import hypercube, walsh_coefficients


def main() -> None:
    n, r, dim = 10, 1, 4
    rng = np.random.default_rng(23)
    masks = [mask for mask in range(1 << n) if mask.bit_count() <= r]
    coefficients = {
        mask: 0.05 * (rng.normal(size=dim) + 1j * rng.normal(size=dim))
        for mask in masks
    }
    coefficients[0] += np.array([3.0, 0.0, 0.0, 0.0])

    raw = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
    effect = raw.conj().T @ raw
    effect /= np.linalg.eigvalsh(effect).max()

    probabilities = {}
    norms = {}
    for x in hypercube(n):
        state = polynomial_state(coefficients, x)
        norms[x] = float(np.vdot(state, state).real)
        probabilities[x] = normalized_born_probability(state, effect)

    q_min, q_max = min(norms.values()), max(norms.values())
    kappa = q_max / q_min
    delta = centered_relative_variation(q_min, q_max)
    coeffs = walsh_coefficients(probabilities, n)

    print(f"n={n} r={r} kappa={kappa:.6f} delta={delta:.6f}")
    print("m cutoff observed_tail_energy theorem_bound")
    for m in range(1, 5):
        cutoff = 2 * r * m
        observed = sum(
            abs(value) ** 2
            for mask, value in coeffs.items()
            if mask.bit_count() > cutoff
        )
        bound = spectral_tail_energy_bound(cutoff, r, kappa)
        print(f"{m} {cutoff:>6} {observed:.6e} {bound:.6e}")

    full_mask = (1 << n) - 1
    observed_full = abs(coeffs[full_mask])
    full_bound = spectral_coefficient_bound(n, r, kappa)
    print(f"full-order coefficient: observed={observed_full:.6e} bound={full_bound:.6e}")


if __name__ == "__main__":
    main()
