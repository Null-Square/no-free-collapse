"""Deterministic random search for asymmetric Gram matrices beating the mean-field family.

This is falsification evidence only, not a proof of global optimality.
"""

from __future__ import annotations

import itertools

import numpy as np

from no_free_collapse.gram import (
    feature_vector,
    monomial_masks,
    optimal_absolute_linear_readout_for_fixed_Q,
)
from no_free_collapse.interactions import monomial
from no_free_collapse.symmetric import mean_field_full_parity_capacity


def cube_condition(Q: np.ndarray, features: np.ndarray) -> float:
    q = np.einsum("bi,ij,bj->b", features.conj(), Q, features).real
    if q.min() <= 0:
        return np.inf
    return float(q.max() / q.min())


def move_to_condition_boundary(
    R: np.ndarray,
    features: np.ndarray,
    kappa: float,
) -> np.ndarray | None:
    """Mix a PSD matrix with identity until cube conditioning reaches kappa."""
    dim = R.shape[0]
    q = np.einsum("bi,ij,bj->b", features.conj(), R, features).real
    R = R * (dim / q.mean())
    if cube_condition(R, features) < kappa:
        return None

    identity = np.eye(dim, dtype=np.complex128)
    lo, hi = 0.0, 1.0
    for _ in range(60):
        t = 0.5 * (lo + hi)
        Q = (1.0 - t) * identity + t * R
        if cube_condition(Q, features) < kappa:
            lo = t
        else:
            hi = t
    return (1.0 - hi) * identity + hi * R


def search(n: int, kappa: float, samples: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    masks = monomial_masks(n, 1)
    xs = list(itertools.product((-1, 1), repeat=n))
    features = np.asarray([feature_vector(x, masks) for x in xs])
    weights = np.asarray([
        monomial(x, (1 << n) - 1) / (1 << n)
        for x in xs
    ])

    best = 0.0
    for sample in range(samples):
        complex_model = bool(sample % 2)
        raw = rng.normal(size=(n + 1, n + 1))
        if complex_model:
            raw = raw + 1j * rng.normal(size=(n + 1, n + 1))
        R = raw.conj().T @ raw
        Q = move_to_condition_boundary(R, features, kappa)
        if Q is None:
            continue
        value, _, _ = optimal_absolute_linear_readout_for_fixed_Q(Q, features, weights)
        best = max(best, value)
    return best


def main() -> None:
    samples = 1000
    print("n kappa random_best mean_field ratio")
    for n in (4, 6):
        for kappa in (2.0, 5.0, 10.0):
            best = search(n, kappa, samples=samples, seed=1000 + 10 * n + int(kappa))
            mean_field = mean_field_full_parity_capacity(n, kappa)
            print(f"{n} {kappa:4.1f} {best:.9f} {mean_field:.9f} {best/mean_field:.6f}")


if __name__ == "__main__":
    main()
