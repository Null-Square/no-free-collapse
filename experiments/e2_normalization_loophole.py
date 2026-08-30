"""Exact counterexample: normalization can create higher-order interaction."""

from no_free_collapse.exact_examples import normalization_loophole_probability
from no_free_collapse.interactions import evaluate_on_hypercube, walsh_coefficients


if __name__ == "__main__":
    coeffs = walsh_coefficients(evaluate_on_hypercube(3, normalization_loophole_probability), 3)
    for mask, value in sorted(coeffs.items(), key=lambda item: (item[0].bit_count(), item[0])):
        print(f"mask={mask:03b} order={mask.bit_count()} coefficient={value}")
