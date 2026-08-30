"""Exact family: order-1 unnormalized amplitudes -> full n-way interaction after normalization."""

from no_free_collapse.exact_examples import normalization_explosion_probability
from no_free_collapse.interactions import evaluate_on_hypercube, walsh_coefficient


if __name__ == "__main__":
    for n in range(1, 13):
        values = evaluate_on_hypercube(n, normalization_explosion_probability)
        full = walsh_coefficient(values, n, (1 << n) - 1)
        print(f"n={n:2d} c={2*n+1:2d} full_coefficient={full} approx={float(full):+.6e}")
