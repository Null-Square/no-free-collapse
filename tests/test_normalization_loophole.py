from fractions import Fraction

from no_free_collapse.exact_examples import normalization_explosion_probability, normalization_loophole_probability
from no_free_collapse.interactions import evaluate_on_hypercube, interaction_degree, walsh_coefficient, walsh_coefficients


def test_input_dependent_normalization_creates_exact_cubic_interaction():
    values = evaluate_on_hypercube(3, normalization_loophole_probability)
    coeffs = walsh_coefficients(values, 3)

    assert coeffs[0b111] == Fraction(-6, 65)
    assert interaction_degree(coeffs, atol=0) == 3


def test_normalization_can_amplify_order_one_to_full_order():
    for n in range(1, 13):
        values = evaluate_on_hypercube(n, normalization_explosion_probability)
        full = walsh_coefficient(values, n, (1 << n) - 1)
        assert full != 0
        assert (full > 0) == (n % 2 == 1)
