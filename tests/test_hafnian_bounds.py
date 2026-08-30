import math

import numpy as np

from no_free_collapse.hafnian_bounds import (
    bounded_quadratic_hafnian_bound,
    disjoint_pair_hafnian_value,
    disjoint_pair_quadratic,
    full_parity_power_coefficient,
    hafnian,
    normalized_cube_max,
    twice_offdiag,
)


def test_full_parity_coefficient_is_factorial_times_hafnian():
    rng = np.random.default_rng(123)
    for size in (4, 6):
        raw = rng.normal(size=(size, size))
        B = raw @ raw.T
        B /= normalized_cube_max(B)
        L = size // 2
        lhs = full_parity_power_coefficient(B)
        rhs = math.factorial(L) * hafnian(twice_offdiag(B))
        assert abs(lhs - rhs) < 1e-10


def test_random_psd_examples_obey_universal_bound():
    rng = np.random.default_rng(77)
    for size in (4, 6, 8):
        bound = bounded_quadratic_hafnian_bound(size)
        for _ in range(12):
            raw = rng.normal(size=(size, size))
            B = raw @ raw.T
            B /= normalized_cube_max(B)
            value = abs(hafnian(twice_offdiag(B)))
            assert value <= bound + 1e-10


def test_four_variable_bound_is_exact():
    B = disjoint_pair_quadratic(4)
    assert abs(normalized_cube_max(B) - 1.0) < 1e-12
    value = abs(hafnian(twice_offdiag(B)))
    assert math.isclose(value, 1.0 / 16.0, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(value, bounded_quadratic_hafnian_bound(4), rel_tol=0.0, abs_tol=1e-12)


def test_pair_witness_matches_closed_form_for_larger_sizes():
    for size in (4, 6, 8, 10):
        B = disjoint_pair_quadratic(size)
        value = abs(hafnian(twice_offdiag(B)))
        assert math.isclose(value, disjoint_pair_hafnian_value(size), rel_tol=1e-12, abs_tol=1e-15)
        assert value <= bounded_quadratic_hafnian_bound(size) + 1e-15
