"""Compare the universal Chebyshev upper bound to an exact symmetric construction."""

from __future__ import annotations

from no_free_collapse.chebyshev import born_chebyshev_interaction_bound
from no_free_collapse.symmetric import mean_field_full_parity_capacity


def main() -> None:
    print("n kappa construction chebyshev_upper gap_ratio")
    for kappa in (2.0, 5.0, 10.0):
        for n in (4, 6, 8, 10, 20):
            construction = mean_field_full_parity_capacity(n, kappa)
            upper = born_chebyshev_interaction_bound(n, 1, kappa)
            ratio = construction / upper
            print(f"{n:2d} {kappa:5.1f} {construction:.9e} {upper:.9e} {ratio:.6f}")
        print()


if __name__ == "__main__":
    main()
