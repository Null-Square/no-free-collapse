import numpy as np

from no_free_collapse.born import normalized_born_probability, polynomial_state
from no_free_collapse.gram import (
    effect_gram,
    feature_vector,
    gram_probability,
    monomial_masks,
    optimal_absolute_linear_readout_for_fixed_Q,
    realize_gram_pair,
)
from no_free_collapse.interactions import hypercube, monomial


def _model(n=5, r=1, dim=4, seed=4):
    rng = np.random.default_rng(seed)
    masks = monomial_masks(n, r)
    coefficients = {
        mask: rng.normal(size=dim) + 1j * rng.normal(size=dim)
        for mask in masks
    }
    raw = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
    effect = raw.conj().T @ raw
    effect /= np.linalg.eigvalsh(effect).max()
    return coefficients, effect, masks


def test_gram_probability_matches_amplitude_model():
    n = 5
    coefficients, effect, masks = _model(n=n)
    A, Q, masks = effect_gram(coefficients, effect, masks)
    for x in hypercube(n):
        state = polynomial_state(coefficients, x)
        direct = normalized_born_probability(state, effect)
        gram = gram_probability(x, A, Q, masks)
        assert abs(direct - gram) < 1e-10


def test_fixed_Q_readout_optimization_is_attained_and_dominates_random_effects():
    n, r = 5, 1
    coefficients, _, masks = _model(n=n, r=r, dim=6, seed=9)
    _, Q, _ = effect_gram(coefficients, np.eye(6), masks)
    xs = list(hypercube(n))
    features = np.asarray([feature_vector(x, masks) for x in xs])
    weights = np.asarray([monomial(x, (1 << n) - 1) / (1 << n) for x in xs])

    optimum, A, sign = optimal_absolute_linear_readout_for_fixed_Q(Q, features, weights)
    attained = sum(
        weights[i] * gram_probability(x, A, Q, masks)
        for i, x in enumerate(xs)
    )
    assert abs(abs(attained) - optimum) < 1e-10
    assert np.sign(attained) == sign or abs(attained) < 1e-12

    rng = np.random.default_rng(123)
    evals, evecs = np.linalg.eigh(0.5 * (Q + Q.conj().T))
    Qh = (evecs * np.sqrt(np.clip(evals, 0.0, None))) @ evecs.conj().T
    for _ in range(50):
        U, _ = np.linalg.qr(rng.normal(size=Q.shape) + 1j * rng.normal(size=Q.shape))
        diag = rng.uniform(0.0, 1.0, size=Q.shape[0])
        C = U @ np.diag(diag) @ U.conj().T
        candidate_A = Qh @ C @ Qh
        value = sum(
            weights[i] * gram_probability(x, candidate_A, Q, masks)
            for i, x in enumerate(xs)
        )
        assert abs(value) <= optimum + 1e-10


def test_gram_pair_has_rank_minimal_realization():
    n = 4
    coefficients, effect, masks = _model(n=n, r=1, dim=5, seed=21)
    A, Q, masks = effect_gram(coefficients, effect, masks)
    rebuilt_coefficients, rebuilt_effect = realize_gram_pair(A, Q, masks)

    latent_dim = len(next(iter(rebuilt_coefficients.values())))
    assert latent_dim == np.linalg.matrix_rank(Q, tol=1e-9)

    for x in hypercube(n):
        original = gram_probability(x, A, Q, masks)
        state = polynomial_state(rebuilt_coefficients, x)
        rebuilt = normalized_born_probability(state, rebuilt_effect)
        assert abs(original - rebuilt) < 1e-9
