"""Exact witnesses show that linear squared-distance stability has optimal order."""
from fractions import Fraction as F

import numpy as np
import pytest

from no_free_collapse.zero_diagonal import exact_certificate, six_hafnian


@pytest.mark.parametrize("t", [F(1, 1000), F(1, 100), F(1, 10)])
def test_rational_rotation_sharpness(t):
    c = (1 - t*t) / (1 + t*t)
    s = 2*t / (1 + t*t)
    u = [[c, -s, F(0)], [s, c, F(0)], [F(0), F(0), F(1)]]
    b = [[F(0) for _ in range(6)] for _ in range(6)]
    for i in range(3):
        for j in range(3):
            b[i][j+3] = b[j+3][i] = u[i][j]
    cert = exact_certificate(b, 1)
    h = c*c - s*s
    eps = 1 - h
    distance_squared = 8*(1-c)
    assert cert["edge_energy"] == 3
    assert cert["subhafnian_energy"] == 2+h*h
    assert cert["defect"] == 1-h*h
    assert distance_squared / eps == 4/(1+c)
    assert distance_squared / cert["defect"] == 2/(c*c*(1+c))
    assert F(2) < distance_squared / eps < F(21, 10)
    assert c*c > F(1, 2)  # The cross-diagonal matching is the unique nearest one.
    assert np.isclose(six_hafnian(np.asarray(b, dtype=float)), float(h))
