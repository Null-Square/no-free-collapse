"""Numerically verify the 2r raw Born-score interaction ceiling."""

import numpy as np

from no_free_collapse.born import born_score, polynomial_state
from no_free_collapse.interactions import evaluate_on_hypercube, max_abs_by_order, walsh_coefficients


def run(n=8, r=2, dim=6, seed=20260830):
    rng = np.random.default_rng(seed)
    coefficients = {
        mask: rng.normal(size=dim) + 1j * rng.normal(size=dim)
        for mask in range(1 << n)
        if mask.bit_count() <= r and rng.random() < 0.35
    }
    a = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
    effect = a.conj().T @ a
    values = evaluate_on_hypercube(n, lambda x: born_score(polynomial_state(coefficients, x), effect))
    return max_abs_by_order(walsh_coefficients(values, n))


if __name__ == "__main__":
    for order, magnitude in sorted(run().items()):
        print(f"order={order}: max_abs={magnitude:.6e}")
