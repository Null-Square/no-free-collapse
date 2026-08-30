"""Exact fixed-norm constructions witnessing tight interaction-order bounds."""

from __future__ import annotations

import numpy as np


def parity_state_coefficients(k: int) -> dict[int, np.ndarray]:
    """A 2D fixed-norm state computing k-way parity with preparation order ceil(k/2).

    Split [k] into A,B and use
        psi=((x_A+x_B)/2, (x_A-x_B)/2).
    Measuring |0><0| yields (1 + prod_i x_i)/2 exactly.
    """
    if k < 1:
        raise ValueError("k must be at least 1")

    left_size = (k + 1) // 2
    left_mask = (1 << left_size) - 1
    right_mask = ((1 << k) - 1) ^ left_mask

    coefficients: dict[int, np.ndarray] = {}
    coefficients[left_mask] = np.array([0.5, 0.5], dtype=np.complex128)
    coefficients[right_mask] = coefficients.get(right_mask, np.zeros(2, dtype=np.complex128)) + np.array(
        [0.5, -0.5], dtype=np.complex128
    )
    return coefficients


def projector_zero(dim: int = 2) -> np.ndarray:
    if dim < 1:
        raise ValueError("dim must be positive")
    effect = np.zeros((dim, dim), dtype=np.complex128)
    effect[0, 0] = 1.0
    return effect
