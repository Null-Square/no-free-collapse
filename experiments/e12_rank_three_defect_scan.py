"""Deterministic CPU scan for the exact rank-three stability margin.

This experiment is not a proof.  It stress-tests the proved defect identity and
reports the smallest observed value of

    8 W + 4 L + 10 S4 - S2 - S2^2

over a reproducible family of random rank-three projections.
"""

import numpy as np

from no_free_collapse.projection_gradient import rank_three_projection_stability_terms


def main() -> None:
    rng = np.random.default_rng(20260830)
    minimum = float("inf")
    minimum_seed = None
    for seed in range(20000):
        Q, _ = np.linalg.qr(rng.normal(size=(6, 3)))
        P = Q @ Q.T
        W, L, S2, S4 = rank_three_projection_stability_terms(P)
        margin = 8.0 * W + 4.0 * L + 10.0 * S4 - S2 - S2 * S2
        if margin < minimum:
            minimum = margin
            minimum_seed = seed
    print(f"minimum stability margin: {minimum:.12g}")
    print(f"sample index: {minimum_seed}")


if __name__ == "__main__":
    main()
