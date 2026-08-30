"""Verify the universal quadratic-form hafnian bound and the exact m=4 case."""

from __future__ import annotations

from no_free_collapse.hafnian_bounds import (
    bounded_quadratic_hafnian_bound,
    disjoint_pair_hafnian_value,
)


def main() -> None:
    print("m pair_witness universal_upper ratio")
    for m in (4, 6, 8, 10, 12):
        lower = disjoint_pair_hafnian_value(m)
        upper = bounded_quadratic_hafnian_bound(m)
        print(f"{m:2d} {lower:.12e} {upper:.12e} {lower/upper:.9f}")


if __name__ == "__main__":
    main()
