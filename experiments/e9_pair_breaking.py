"""Compare symmetry-breaking pairs, symmetric mean-field, and the universal bound."""

from no_free_collapse.chebyshev import born_chebyshev_interaction_bound
from no_free_collapse.paired import matched_pair_asymptotic_rate, matched_pair_capacity
from no_free_collapse.symmetric import mean_field_asymptotic_rate, mean_field_full_parity_capacity


def main() -> None:
    print("n kappa pair mean_field chebyshev_upper best_lower")
    for kappa in (1.5, 2.0, 5.0, 10.0):
        for n in (4, 6, 8, 10, 20):
            pair = matched_pair_capacity(n, kappa)
            mean = mean_field_full_parity_capacity(n, kappa)
            upper = born_chebyshev_interaction_bound(n, 1, kappa)
            print(f"{n:2d} {kappa:5.1f} {pair:.9e} {mean:.9e} {upper:.9e} {max(pair,mean):.9e}")
        print(
            f"rates: pair={matched_pair_asymptotic_rate(kappa):.8f} "
            f"mean_field={mean_field_asymptotic_rate(kappa):.8f}"
        )
        print()


if __name__ == "__main__":
    main()
