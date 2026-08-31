from fractions import Fraction
from math import comb

import numpy as np

from no_free_collapse.rank_two_projection import (
    rank_two_heavy_bernstein_coefficients,
    rank_two_projection_direct_defect,
    rank_two_stability_decomposition,
)


def _random_rank_two_projection(seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.normal(size=(6, 2)))
    return Q @ Q.T


def _random_orthonormal_pair(seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.normal(size=(5, 2)))
    return Q[:, 0], Q[:, 1]


def test_rank_two_stability_decomposition_is_exact():
    for seed in range(50):
        P = _random_rank_two_projection(seed)
        F, commutator_energy, angle_mixing = rank_two_stability_decomposition(P)
        reconstructed = (F + 8.0 * commutator_energy + 8.0 * angle_mixing) / 8.0
        assert commutator_energy >= -1e-12
        assert angle_mixing >= -1e-12
        assert np.isclose(
            reconstructed,
            rank_two_projection_direct_defect(P),
            atol=3e-11,
            rtol=3e-11,
        )


def test_first_heavy_bernstein_coefficient_is_nonnegative_randomly():
    for seed in range(100):
        p, q = _random_orthonormal_pair(seed)
        coeffs = rank_two_heavy_bernstein_coefficients(p, q)
        assert coeffs[1] >= -2e-11


def test_first_heavy_bernstein_coefficient_coordinate_equality_family():
    # q is a coordinate vector and p is any unit vector on the orthogonal
    # four-coordinate support.  This is the sharp b1=1/32 family seen in the
    # analytic formula before the deliberately looser positivity certificate.
    rng = np.random.default_rng(2026)
    for _ in range(10):
        p = np.zeros(5)
        p[:4] = rng.normal(size=4)
        p /= np.linalg.norm(p)
        q = np.asarray([0.0, 0.0, 0.0, 0.0, 1.0])
        coeffs = rank_two_heavy_bernstein_coefficients(p, q)
        assert np.isclose(coeffs[1], 1.0 / 32.0, atol=2e-12, rtol=2e-12)


# Exact rational certificate for the final quartic simplex lemma in the b1
# proof.  After the Timofte--Riener half-degree reduction, the four remaining
# probability coordinates have multiplicities 2+2 or 1+3.  With x in [0,1]
# and a normalized split parameter t in [0,1], the two quartics below are the
# exact upper-bound polynomials that must be <=0.
P22 = {
    (4, 4): Fraction(6), (3, 4): Fraction(-24), (2, 4): Fraction(36),
    (1, 4): Fraction(-24), (0, 4): Fraction(6),
    (4, 3): Fraction(-12), (3, 3): Fraction(48), (2, 3): Fraction(-72),
    (1, 3): Fraction(48), (0, 3): Fraction(-12),
    (4, 2): Fraction(42), (3, 2): Fraction(-147), (2, 2): Fraction(173),
    (1, 2): Fraction(-73), (0, 2): Fraction(5),
    (4, 1): Fraction(-36), (3, 1): Fraction(123), (2, 1): Fraction(-137),
    (1, 1): Fraction(49), (0, 1): Fraction(1),
    (4, 0): Fraction(-27), (3, 0): Fraction(18), (2, 0): Fraction(12),
    (1, 0): Fraction(-9), (0, 0): Fraction(-2),
}

_P13_BRACKET = {
    (4, 4): Fraction(96), (3, 4): Fraction(-384), (2, 4): Fraction(576),
    (1, 4): Fraction(-384), (0, 4): Fraction(96),
    (4, 3): Fraction(-480), (3, 3): Fraction(1632), (2, 3): Fraction(-2016),
    (1, 3): Fraction(1056), (0, 3): Fraction(-192),
    (4, 2): Fraction(624), (3, 2): Fraction(-1932), (2, 2): Fraction(2116),
    (1, 2): Fraction(-932), (0, 2): Fraction(124),
    (4, 1): Fraction(-288), (3, 1): Fraction(792), (2, 1): Fraction(-744),
    (1, 1): Fraction(264), (0, 1): Fraction(-24),
    (4, 0): Fraction(144), (3, 0): Fraction(-216), (2, 0): Fraction(96),
    (1, 0): Fraction(-3), (0, 0): Fraction(3),
}
P13 = {key: -value / 3 for key, value in _P13_BRACKET.items()}


def _affine_pullback(poly, x0, x1, t0, t1):
    """Substitute x=x0+(x1-x0)X and t=t0+(t1-t0)T exactly."""
    out = {}
    dx = x1 - x0
    dt = t1 - t0
    for (a, b), coefficient in poly.items():
        for i in range(a + 1):
            cx = Fraction(comb(a, i)) * (x0 ** (a - i)) * (dx ** i)
            for j in range(b + 1):
                ct = Fraction(comb(b, j)) * (t0 ** (b - j)) * (dt ** j)
                out[(i, j)] = out.get((i, j), Fraction(0)) + coefficient * cx * ct
    return out


def _bernstein_coefficients(poly, nx=4, nt=4):
    """Convert bivariate power coefficients to the tensor Bernstein basis."""
    values = []
    for i in range(nx + 1):
        for j in range(nt + 1):
            value = Fraction(0)
            for a in range(i + 1):
                for b in range(j + 1):
                    value += (
                        poly.get((a, b), Fraction(0))
                        * Fraction(comb(i, a), comb(nx, a))
                        * Fraction(comb(j, b), comb(nt, b))
                    )
            values.append(value)
    return values


def _cell_max_bernstein(poly, xr, tr):
    pulled = _affine_pullback(poly, *xr, *tr)
    return max(_bernstein_coefficients(pulled))


def test_exact_bernstein_certificate_for_two_plus_two_pattern():
    assert _cell_max_bernstein(
        P22, (Fraction(0), Fraction(1, 2)), (Fraction(0), Fraction(1))
    ) == Fraction(-1, 4)
    assert _cell_max_bernstein(
        P22, (Fraction(1, 2), Fraction(1)), (Fraction(0), Fraction(1))
    ) == Fraction(-17, 12)


def test_exact_bernstein_certificate_for_one_plus_three_pattern():
    assert _cell_max_bernstein(
        P13, (Fraction(1, 2), Fraction(1)), (Fraction(0), Fraction(1))
    ) == Fraction(-7, 6)
    assert _cell_max_bernstein(
        P13, (Fraction(0), Fraction(1, 2)), (Fraction(1, 2), Fraction(1))
    ) == Fraction(-7, 8)
    assert _cell_max_bernstein(
        P13, (Fraction(0), Fraction(1, 2)), (Fraction(0), Fraction(1, 2))
    ) == Fraction(0)
