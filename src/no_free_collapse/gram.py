"""Gram-matrix characterization of polynomial Born reasoning models."""

from __future__ import annotations

from typing import Sequence

import numpy as np

from .interactions import monomial

Array = np.ndarray


def monomial_masks(n: int, order: int) -> list[int]:
    """All Boolean monomial masks of degree at most ``order``."""
    if n < 0 or order < 0:
        raise ValueError("n and order must be non-negative")
    order = min(order, n)
    return [mask for mask in range(1 << n) if mask.bit_count() <= order]


def feature_vector(x: Sequence[int], masks: Sequence[int]) -> Array:
    """Evaluate the selected Walsh monomials at x."""
    return np.asarray([monomial(x, mask) for mask in masks], dtype=np.float64)


def coefficient_gram(
    coefficients: dict[int, Array],
    masks: Sequence[int] | None = None,
) -> tuple[Array, list[int]]:
    """Return Q=U^*U for coefficient vectors arranged by monomial mask."""
    if not coefficients:
        raise ValueError("at least one coefficient vector is required")
    if masks is None:
        masks = sorted(coefficients)
    masks = list(masks)
    first = np.asarray(next(iter(coefficients.values())), dtype=np.complex128)
    U = np.column_stack([
        np.asarray(coefficients.get(mask, np.zeros_like(first)), dtype=np.complex128)
        for mask in masks
    ])
    return U.conj().T @ U, masks


def effect_gram(
    coefficients: dict[int, Array],
    effect: Array,
    masks: Sequence[int] | None = None,
) -> tuple[Array, Array, list[int]]:
    """Return (A,Q,masks) with A=U^* M U and Q=U^*U."""
    if not coefficients:
        raise ValueError("at least one coefficient vector is required")
    if masks is None:
        masks = sorted(coefficients)
    masks = list(masks)
    first = np.asarray(next(iter(coefficients.values())), dtype=np.complex128)
    U = np.column_stack([
        np.asarray(coefficients.get(mask, np.zeros_like(first)), dtype=np.complex128)
        for mask in masks
    ])
    effect = np.asarray(effect, dtype=np.complex128)
    if effect.shape != (U.shape[0], U.shape[0]):
        raise ValueError("effect dimension does not match coefficient vectors")
    Q = U.conj().T @ U
    A = U.conj().T @ effect @ U
    return A, Q, masks


def gram_probability(x: Sequence[int], A: Array, Q: Array, masks: Sequence[int]) -> float:
    """Evaluate p(x)=z^* A z / (z^* Q z)."""
    z = feature_vector(x, masks).astype(np.complex128)
    denominator = np.vdot(z, Q @ z)
    if denominator.real <= 0 or abs(denominator.imag) > 1e-9:
        raise ValueError("Q gives a non-positive or complex denominator")
    numerator = np.vdot(z, A @ z)
    if abs(numerator.imag) > 1e-9:
        raise ValueError("A gives a complex numerator")
    return float(numerator.real / denominator.real)


def realize_gram_pair(
    A: Array,
    Q: Array,
    masks: Sequence[int],
    *,
    atol: float = 1e-10,
) -> tuple[dict[int, Array], Array]:
    """Realize 0<=A<=Q using the minimum latent dimension rank(Q)."""
    A = np.asarray(A, dtype=np.complex128)
    Q = np.asarray(Q, dtype=np.complex128)
    if A.shape != Q.shape or Q.shape != (len(masks), len(masks)):
        raise ValueError("A, Q, and masks have incompatible dimensions")
    A = 0.5 * (A + A.conj().T)
    Q = 0.5 * (Q + Q.conj().T)
    if np.linalg.eigvalsh(A).min(initial=0.0) < -atol:
        raise ValueError("A must be positive semidefinite")
    if np.linalg.eigvalsh(Q - A).min(initial=0.0) < -atol:
        raise ValueError("A must satisfy A <= Q")

    eigenvalues, eigenvectors = np.linalg.eigh(Q)
    if eigenvalues.min(initial=0.0) < -atol:
        raise ValueError("Q must be positive semidefinite")
    keep = eigenvalues > atol
    if not np.any(keep):
        raise ValueError("Q must have nonzero rank")
    lam = eigenvalues[keep]
    V = eigenvectors[:, keep]
    U = np.sqrt(lam)[:, None] * V.conj().T

    reduced_A = V.conj().T @ A @ V
    effect = reduced_A / np.sqrt(np.outer(lam, lam))
    effect = 0.5 * (effect + effect.conj().T)
    if np.linalg.eigvalsh(effect).min(initial=0.0) < -1e-8:
        raise ValueError("constructed effect is not positive semidefinite")
    if np.linalg.eigvalsh(np.eye(effect.shape[0]) - effect).min(initial=0.0) < -1e-8:
        raise ValueError("constructed effect exceeds the identity")

    coefficients = {mask: U[:, j].copy() for j, mask in enumerate(masks)}
    return coefficients, effect


def _psd_sqrt(matrix: Array, *, atol: float = 1e-10) -> Array:
    matrix = np.asarray(matrix, dtype=np.complex128)
    matrix = 0.5 * (matrix + matrix.conj().T)
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    if eigenvalues.min(initial=0.0) < -atol:
        raise ValueError("matrix must be positive semidefinite")
    eigenvalues = np.clip(eigenvalues, 0.0, None)
    return (eigenvectors * np.sqrt(eigenvalues)) @ eigenvectors.conj().T


def optimal_linear_readout_for_fixed_Q(
    Q: Array,
    features: Array,
    weights: Array,
    *,
    atol: float = 1e-10,
) -> tuple[float, Array]:
    """Maximize a weighted output functional over all Born effects for fixed Q."""
    Q = np.asarray(Q, dtype=np.complex128)
    features = np.asarray(features, dtype=np.complex128)
    weights = np.asarray(weights, dtype=np.float64)
    if features.ndim != 2 or features.shape[0] != weights.shape[0]:
        raise ValueError("features and weights have incompatible shapes")
    if Q.shape != (features.shape[1], features.shape[1]):
        raise ValueError("Q dimension does not match feature width")

    Qh = _psd_sqrt(Q, atol=atol)
    denominators = np.einsum("bi,ij,bj->b", features.conj(), Q, features).real
    if np.min(denominators) <= atol:
        raise ValueError("Q must give a strictly positive denominator on every input")

    B = np.einsum(
        "b,bi,bj->ij",
        weights / denominators,
        features.conj(),
        features,
    )
    H = Qh @ B @ Qh
    H = 0.5 * (H + H.conj().T)
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    positive = eigenvalues > atol
    projector = (
        eigenvectors[:, positive] @ eigenvectors[:, positive].conj().T
        if np.any(positive)
        else np.zeros_like(Q)
    )
    A = Qh @ projector @ Qh
    optimum = float(eigenvalues[eigenvalues > 0].sum())
    return optimum, A


def optimal_absolute_linear_readout_for_fixed_Q(
    Q: Array,
    features: Array,
    weights: Array,
    *,
    atol: float = 1e-10,
) -> tuple[float, Array, int]:
    """Maximize the absolute weighted output functional for fixed Q."""
    pos_value, pos_A = optimal_linear_readout_for_fixed_Q(Q, features, weights, atol=atol)
    neg_value, neg_A = optimal_linear_readout_for_fixed_Q(
        Q,
        features,
        -np.asarray(weights),
        atol=atol,
    )
    if pos_value >= neg_value:
        return pos_value, pos_A, 1
    return neg_value, neg_A, -1
