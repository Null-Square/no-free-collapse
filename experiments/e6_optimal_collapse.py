"""Solve the best Born measurement exactly for a fixed preparation Gram matrix."""

from __future__ import annotations

import numpy as np

from no_free_collapse.gram import (
    effect_gram,
    feature_vector,
    gram_probability,
    monomial_masks,
    optimal_absolute_linear_readout_for_fixed_Q,
)
from no_free_collapse.interactions import hypercube, monomial


def main() -> None:
    n, r, dim = 6, 1, 7
    rng = np.random.default_rng(31)
    masks = monomial_masks(n, r)
    coefficients = {
        mask: rng.normal(size=dim) + 1j * rng.normal(size=dim)
        for mask in masks
    }

    _, Q, _ = effect_gram(coefficients, np.eye(dim), masks)
    xs = list(hypercube(n))
    features = np.asarray([feature_vector(x, masks) for x in xs])
    full_mask = (1 << n) - 1
    weights = np.asarray([monomial(x, full_mask) / (1 << n) for x in xs])

    optimum, A, sign = optimal_absolute_linear_readout_for_fixed_Q(Q, features, weights)
    attained = sum(
        weights[i] * gram_probability(x, A, Q, masks)
        for i, x in enumerate(xs)
    )

    denominators = np.einsum("bi,ij,bj->b", features, Q, features).real
    kappa = denominators.max() / denominators.min()

    print(f"n={n} r={r} feature_dim={len(masks)} rank(Q)={np.linalg.matrix_rank(Q)}")
    print(f"norm condition number kappa={kappa:.6f}")
    print(f"exact best |full-order Walsh coefficient|={optimum:.12f}")
    print(f"attained coefficient={attained:.12f} sign={sign:+d}")
    print(f"verification error={abs(abs(attained)-optimum):.3e}")


if __name__ == "__main__":
    main()
