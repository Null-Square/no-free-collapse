from fractions import Fraction

from no_free_collapse.exact_examples import normalization_loophole_probability
from no_free_collapse.interactions import evaluate_on_hypercube, interaction_degree, walsh_coefficients


def test_input_dependent_normalization_creates_exact_cubic_interaction():
    values = evaluate_on_hypercube(3, normalization_loophole_probability)
    coeffs = walsh_coefficients(values, 3)

    assert coeffs[0b111] == Fraction(-6, 65)
    assert interaction_degree(coeffs, atol=0) == 3
