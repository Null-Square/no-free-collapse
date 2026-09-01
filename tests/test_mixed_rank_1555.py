from fractions import Fraction
from math import comb

import numpy as np

from no_free_collapse.mixed_nested import rank_one_rank_five_mixed_formula
from no_free_collapse.spectral_polarization import nested_projection_polarized_defect


def _bernstein_coefficients_2d(coeffs, degree_x, degree_y):
    values = []
    for k in range(degree_x + 1):
        row = []
        for ell in range(degree_y + 1):
            total = Fraction(0)
            for i in range(k + 1):
                for j in range(ell + 1):
                    a = coeffs.get((i, j), Fraction(0))
                    total += (
                        a
                        * Fraction(comb(k, i), comb(degree_x, i))
                        * Fraction(comb(ell, j), comb(degree_y, j))
                    )
            row.append(total)
        values.append(row)
    return values


def _random_orthonormal_pair(seed):
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.normal(size=(6, 2)))
    return Q[:, 0], Q[:, 1]


def test_rank_one_rank_five_formula_matches_kernel_and_is_nonnegative():
    for seed in range(200):
        u, v = _random_orthonormal_pair(seed)
        P = np.outer(u, u)
        Q5 = np.eye(6) - np.outer(v, v)
        direct = nested_projection_polarized_defect([P, Q5, Q5, Q5])
        formula = rank_one_rank_five_mixed_formula(u, v)
        assert np.isclose(direct, formula, atol=3e-12, rtol=3e-12)
        assert formula >= -3e-12


def test_rank_one_rank_five_complement_pattern_matches():
    for seed in range(50):
        u, v = _random_orthonormal_pair(1000 + seed)
        P1 = np.outer(u, u)
        Q5 = np.eye(6) - np.outer(v, v)
        left = nested_projection_polarized_defect([P1, Q5, Q5, Q5])
        complement = nested_projection_polarized_defect(
            [
                np.eye(6) - Q5,
                np.eye(6) - Q5,
                np.eye(6) - Q5,
                np.eye(6) - P1,
            ]
        )
        assert np.isclose(left, complement, atol=3e-12, rtol=3e-12)


def test_light_coordinate_bernstein_certificate_is_exactly_nonnegative():
    # After x=1/6+p/3 and r=x+(1/2-x)q, the lower polynomial G has
    # this exact power-basis representation in (p,q).
    coeffs = {
        (3, 2): Fraction(5, 18),
        (3, 1): Fraction(-4, 9),
        (3, 0): Fraction(-8, 9),
        (2, 2): Fraction(-5, 9),
        (2, 1): Fraction(14, 9),
        (2, 0): Fraction(-2, 9),
        (1, 2): Fraction(5, 18),
        (1, 1): Fraction(-10, 9),
        (1, 0): Fraction(7, 9),
        (0, 0): Fraction(1, 3),
    }
    bernstein = _bernstein_coefficients_2d(coeffs, 3, 2)
    expected = [
        [Fraction(1, 3), Fraction(1, 3), Fraction(1, 3)],
        [Fraction(16, 27), Fraction(11, 27), Fraction(17, 54)],
        [Fraction(7, 9), Fraction(2, 3), Fraction(5, 9)],
        [Fraction(0), Fraction(0), Fraction(0)],
    ]
    assert bernstein == expected
    assert min(value for row in bernstein for value in row) >= 0


def test_heavy_scalar_bernstein_certificate_is_exactly_negative():
    # After r=X/2 and z=(1+4Y)/5, the polynomial B in the proof has
    # this exact power-basis representation in (X,Y).
    coeffs = {
        (4, 2): Fraction(36, 25),
        (4, 1): Fraction(-207, 25),
        (4, 0): Fraction(-54, 25),
        (3, 2): Fraction(-12, 5),
        (3, 1): Fraction(39, 1),
        (3, 0): Fraction(99, 10),
        (2, 1): Fraction(-288, 5),
        (2, 0): Fraction(-189, 10),
        (1, 1): Fraction(122, 5),
        (1, 0): Fraction(108, 5),
        (0, 0): Fraction(-13, 1),
    }
    bernstein = _bernstein_coefficients_2d(coeffs, 4, 2)
    expected = [
        [Fraction(-13), Fraction(-13), Fraction(-13)],
        [Fraction(-38, 5), Fraction(-91, 20), Fraction(-3, 2)],
        [Fraction(-107, 20), Fraction(-81, 20), Fraction(-11, 4)],
        [Fraction(-151, 40), Fraction(-83, 20), Fraction(-41, 8)],
        [Fraction(-64, 25), Fraction(-19, 5), Fraction(-6)],
    ]
    assert bernstein == expected
    assert max(value for row in bernstein for value in row) < 0


def test_coordinate_rank_five_boundary_has_zero_mixed_defect():
    v = np.zeros(6)
    v[0] = 1.0
    u = np.zeros(6)
    u[1] = 1.0
    assert np.isclose(rank_one_rank_five_mixed_formula(u, v), 0.0, atol=1e-14)
