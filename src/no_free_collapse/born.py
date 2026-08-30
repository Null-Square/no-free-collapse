"""Polynomial amplitude preparations and Born-style readout."""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np

from .interactions import monomial

Array = np.ndarray


def preparation_order(coefficients: Mapping[int, Array]) -> int:
    """Maximum Boolean monomial order appearing in an amplitude preparation."""
    return max((mask.bit_count() for mask in coefficients), default=0)


def polynomial_state(coefficients: Mapping[int, Array], x: Sequence[int]) -> Array:
    """Evaluate psi_tilde(x)=sum_A x_A u_A."""
    if not coefficients:
        raise ValueError("at least one coefficient vector is required")
    first = next(iter(coefficients.values()))
    state = np.zeros_like(np.asarray(first), dtype=np.complex128)
    for mask, vector in coefficients.items():
        state = state + monomial(x, mask) * np.asarray(vector, dtype=np.complex128)
    return state


def _validate_effect(effect: Array, dim: int) -> Array:
    effect = np.asarray(effect, dtype=np.complex128)
    if effect.shape != (dim, dim):
        raise ValueError(f"effect must have shape {(dim, dim)}, got {effect.shape}")
    return effect


def born_score(state: Array, effect: Array) -> float:
    """Unnormalized quadratic Born score psi^* M psi.

    For a fixed-norm preparation this is proportional to the physical Born probability.
    """
    state = np.asarray(state, dtype=np.complex128)
    effect = _validate_effect(effect, state.shape[0])
    score = np.vdot(state, effect @ state)
    if abs(score.imag) > 1e-9:
        raise ValueError("Born score is unexpectedly complex; effect may not be Hermitian")
    return float(score.real)


def normalized_born_probability(state: Array, effect: Array) -> float:
    """Born probability after explicit input-dependent normalization.

    Important: normalization is nonlinear. If ||psi_tilde(x)|| depends on x, this
    operation can create interaction orders higher than the quadratic raw-score bound.
    """
    state = np.asarray(state, dtype=np.complex128)
    norm_sq = float(np.vdot(state, state).real)
    if norm_sq <= 0:
        raise ValueError("state has zero norm")
    return born_score(state, effect) / norm_sq


def squared_norm(state: Array) -> float:
    state = np.asarray(state, dtype=np.complex128)
    return float(np.vdot(state, state).real)
