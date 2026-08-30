"""Verify the ceil(k/2) fixed-norm construction for k-way parity."""

from itertools import product

import numpy as np

from no_free_collapse.born import born_score, polynomial_state, preparation_order, squared_norm
from no_free_collapse.constructions import parity_state_coefficients, projector_zero
from no_free_collapse.interactions import interaction_degree, walsh_coefficients


if __name__ == "__main__":
    effect = projector_zero()
    for k in range(1, 11):
        coefficients = parity_state_coefficients(k)
        values = {}
        max_norm_error = 0.0
        max_probability_error = 0.0
        for x in product((-1, 1), repeat=k):
            state = polynomial_state(coefficients, x)
            p = born_score(state, effect)
            target = (1 + np.prod(x)) / 2
            values[x] = p
            max_norm_error = max(max_norm_error, abs(squared_norm(state) - 1.0))
            max_probability_error = max(max_probability_error, abs(p - target))
        degree = interaction_degree(walsh_coefficients(values, k), atol=1e-10)
        print(
            f"k={k:2d} prep_order={preparation_order(coefficients):2d} "
            f"output_degree={degree:2d} norm_err={max_norm_error:.1e} p_err={max_probability_error:.1e}"
        )
