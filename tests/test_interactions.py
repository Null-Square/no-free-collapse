from fractions import Fraction

from no_free_collapse.interactions import evaluate_on_hypercube, interaction_degree, walsh_coefficients


def test_exact_parity_has_only_top_order_coefficient():
    n = 5
    values = evaluate_on_hypercube(n, lambda x: Fraction(1 + __import__("math").prod(x), 2))
    coeffs = walsh_coefficients(values, n)

    assert coeffs[0] == Fraction(1, 2)
    assert coeffs[(1 << n) - 1] == Fraction(1, 2)
    assert all(v == 0 for mask, v in coeffs.items() if mask not in (0, (1 << n) - 1))
    assert interaction_degree(coeffs, atol=0) == n
