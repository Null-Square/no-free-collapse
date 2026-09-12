"""Exact linear-map checks, rational PSD certificates and sharpness examples."""
from fractions import Fraction as F
from itertools import combinations

import numpy as np
import pytest

from no_free_collapse.matching_dilation import (
    _compress, dilation_incidence, dilation_operator, dilation_psd_slacks,
    exact_six_psd_certificate, matching_operator,
)
from no_free_collapse.zero_diagonal import subhafnian_energy_direct


@pytest.mark.parametrize("n", [4, 6, 8, 10])
def test_exact_integer_linear_map_on_every_symmetric_basis(n):
    cols, signs = dilation_incidence(n)
    d = int(cols.max())+1
    rows = len(cols)*n
    rp = np.zeros((rows, d), dtype=np.int64)
    rm = rp.copy()
    rp[np.arange(rows), cols.ravel()] = 1
    rm[np.arange(rows), cols.ravel()] = signs.ravel()
    assert np.array_equal(rp.T@rp, (n//2+1)*np.eye(d, dtype=np.int64))
    assert np.array_equal(rm.T@rm, rp.T@rp)
    for i in range(n):
        for j in range(i, n):
            a = np.zeros((n, n), dtype=np.int64)
            a[i, j] = a[j, i] = 1
            # Integer coefficient identity: no sampled entries or tolerance.
            difference = (_compress(a, cols, np.ones_like(signs), d)
                          - _compress(a, cols, signs, d))
            direct = matching_operator(a).astype(np.int64)
            assert np.array_equal(difference, 4*direct)


@pytest.mark.parametrize("n", [4, 6, 8, 10])
def test_random_spectral_diameter_and_psd_slacks(n):
    rng = np.random.default_rng(123+n)
    for _ in range(12):
        a = rng.normal(size=(n, n)); a = (a+a.T)/2
        es = np.linalg.eigvalsh(a)
        t = matching_operator(a)
        assert np.allclose(t, dilation_operator(a), atol=1e-12)
        coefficient = (n//2+1)/4
        assert np.linalg.norm(t, 2) <= coefficient*(es[-1]-es[0])+1e-10
        b = (a-es[0]*np.eye(n))/(es[-1]-es[0])
        tb = matching_operator(b)
        minus, plus = dilation_psd_slacks(b, 1)
        assert np.allclose(minus, coefficient*np.eye(len(tb))-tb, atol=1e-12)
        assert np.allclose(plus, coefficient*np.eye(len(tb))+tb, atol=1e-12)
        assert np.linalg.eigvalsh(minus)[0] >= -1e-10
        assert np.linalg.eigvalsh(plus)[0] >= -1e-10


@pytest.mark.parametrize("n", [4, 6, 8, 10])
def test_sharp_general_dimension_constant(n):
    p = np.ones((n, n))/n
    assert np.isclose(np.linalg.norm(matching_operator(p), 2), (n//2+1)/4)


@pytest.mark.parametrize("family", ["rank_one", "rank_five", "pairs", "diagonal"])
def test_exact_six_equality_families(family):
    sign = [1, -1, 1, 1, -1, 1]
    if family in ("rank_one", "rank_five"):
        a = [[F(sign[i]*sign[j], 6) for j in range(6)] for i in range(6)]
        if family == "rank_five":
            a = [[F(i == j)-a[i][j] for j in range(6)] for i in range(6)]
    elif family == "pairs":
        a = [[F(1, 2) if i == j else F(0) for j in range(6)] for i in range(6)]
        for i in range(0, 6, 2):
            a[i][i+1] = a[i+1][i] = F(-1 if i == 2 else 1, 2)
    else:
        a = [[F(i, 5) if i == j else F(0) for j in range(6)] for i in range(6)]
    cert = exact_six_psd_certificate(a, 1)
    assert cert["gradient_defect"] == 0
    assert cert["hafnian_energy_defect"] == 0


def test_exact_nonprojective_rational_psd_instances():
    rng = np.random.default_rng(392)
    for _ in range(12):
        z = rng.integers(-3, 4, size=(6, 6))
        q = z@z.T
        ell = int(np.trace(q))
        a = [[F(int(x), ell) for x in row] for row in q]
        cert = exact_six_psd_certificate(a, 1)
        assert cert["gradient_defect"] >= 0
        assert cert["dilation_residual"] == 0


def test_global_gradient_on_all_projection_ranks_and_nonprojective_inputs():
    rng = np.random.default_rng(3029)
    for rank in range(1, 7):
        for _ in range(30):
            u = np.linalg.qr(rng.normal(size=(6, 6)))[0]
            eig = np.zeros(6); eig[:rank] = rng.uniform(.1, 1, size=rank); eig[0] = 1
            for values in (eig, (eig > 0).astype(float)):
                a = (u*values)@u.T
                off = a-np.diag(np.diag(a))
                q1 = float(np.sum(off*off)/2)
                q2 = subhafnian_energy_direct(off)
                assert q2 <= q1/4+1e-11


def test_conference_scaling_does_not_refute_psd_operator_bound():
    # Symmetric Paley conference matrix of order six, constructed modulo five.
    s = np.zeros((6, 6), dtype=int)
    s[0, 1:] = 1; s[1:, 0] = 1
    for i in range(5):
        for j in range(5):
            if i != j:
                s[i+1, j+1] = 1 if (i-j)%5 in (1, 4) else -1
    assert np.array_equal(s@s, 5*np.eye(6, dtype=int))
    k = s/np.sqrt(5)
    p = (np.eye(6)+k)/2
    assert np.isclose(np.linalg.norm(matching_operator(k), 2), 3/np.sqrt(5))
    assert np.linalg.norm(matching_operator(k), 2) > 1
    assert np.isclose(np.linalg.norm(matching_operator(p), 2), 3/(2*np.sqrt(5)))
    assert np.linalg.norm(matching_operator(p), 2) < 1


def test_invalid_input_and_exact_type_boundaries():
    for a in (np.zeros((5, 5)), np.zeros((6, 7)), np.zeros((12, 12)),
              np.eye(6, dtype=complex), np.full((6, 6), np.nan)):
        with pytest.raises(ValueError): matching_operator(a)
    a = np.eye(6); a[0, 1] = 1
    with pytest.raises(ValueError): matching_operator(a)
    with pytest.raises(ValueError): dilation_psd_slacks(-np.eye(6), 1)
    with pytest.raises(ValueError): dilation_psd_slacks(np.eye(6), .5)
    with pytest.raises(ValueError): dilation_psd_slacks(np.eye(6), np.nan)
    with pytest.raises(TypeError): exact_six_psd_certificate(np.eye(6).tolist())
    with pytest.raises(ValueError): exact_six_psd_certificate((-np.eye(6,dtype=int)).tolist())
    with pytest.raises(ValueError): exact_six_psd_certificate(np.eye(6,dtype=int).tolist(), 0)
    with pytest.raises(ValueError): exact_six_psd_certificate([[0]*4]*4)
