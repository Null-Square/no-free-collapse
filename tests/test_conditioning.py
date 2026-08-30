import numpy as np

from no_free_collapse.conditioning import (
    centered_relative_variation,
    delta_from_condition_number,
    required_condition_number_for_coefficient,
    spectral_coefficient_bound,
    spectral_tail_energy_bound,
    truncation_bound,
)
from no_free_collapse.born import born_score, normalized_born_probability, polynomial_state
from no_free_collapse.interactions import hypercube, walsh_coefficients


def _random_model(n=8, r=2, dim=4, seed=7):
    rng = np.random.default_rng(seed)
    masks = [mask for mask in range(1 << n) if mask.bit_count() <= r]
    coefficients = {
        mask: 0.04 * (rng.normal(size=dim) + 1j * rng.normal(size=dim))
        for mask in masks
    }
    coefficients[0] += np.array([4.0, 0.0, 0.0, 0.0])
    raw = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
    effect = raw.conj().T @ raw
    effect /= np.linalg.eigvalsh(effect).max()
    return coefficients, effect


def test_condition_number_delta_identity():
    for kappa in (1.0, 1.5, 2.0, 10.0):
        delta = delta_from_condition_number(kappa)
        assert np.isclose((1 + delta) / (1 - delta), kappa)


def test_geometric_denormalization_bound_and_degree():
    n, r = 8, 2
    coefficients, effect = _random_model(n=n, r=r)
    a = {}
    q = {}
    p = {}
    for x in hypercube(n):
        state = polynomial_state(coefficients, x)
        q[x] = float(np.vdot(state, state).real)
        a[x] = born_score(state, effect)
        p[x] = normalized_born_probability(state, effect)

    q_min, q_max = min(q.values()), max(q.values())
    q_center = 0.5 * (q_min + q_max)
    delta = centered_relative_variation(q_min, q_max)
    kappa = q_max / q_min

    for steps in (1, 2):
        approximant = {
            x: (a[x] / q_center)
            * sum((-(q[x] / q_center - 1.0)) ** t for t in range(steps))
            for x in p
        }
        error = max(abs(p[x] - approximant[x]) for x in p)
        degree_bound, error_bound = truncation_bound(r, steps, kappa)
        assert degree_bound == 2 * r * steps
        assert error <= error_bound + 1e-12

        coeffs = walsh_coefficients(approximant, n)
        above = [abs(value) for mask, value in coeffs.items() if mask.bit_count() > degree_bound]
        assert max(above, default=0.0) < 1e-9


def test_normalized_born_spectral_tail_obeys_conditioning_bound():
    n, r = 8, 2
    coefficients, effect = _random_model(n=n, r=r, seed=11)
    p, q = {}, {}
    for x in hypercube(n):
        state = polynomial_state(coefficients, x)
        q[x] = float(np.vdot(state, state).real)
        p[x] = normalized_born_probability(state, effect)

    kappa = max(q.values()) / min(q.values())
    coeffs = walsh_coefficients(p, n)

    for steps in (1, 2):
        cutoff = 2 * r * steps
        tail_energy = sum(
            abs(value) ** 2
            for mask, value in coeffs.items()
            if mask.bit_count() > cutoff
        )
        assert tail_energy <= spectral_tail_energy_bound(cutoff, r, kappa) + 1e-12

    for mask, value in coeffs.items():
        order = mask.bit_count()
        assert abs(value) <= spectral_coefficient_bound(order, r, kappa) + 1e-10


def test_inverse_conditioning_bound():
    alpha = 0.2
    order = 9
    r = 1
    kappa = required_condition_number_for_coefficient(alpha, order, r)
    assert spectral_coefficient_bound(order, r, kappa) >= alpha - 1e-12
    assert spectral_coefficient_bound(order, r, kappa * 0.999) < alpha
