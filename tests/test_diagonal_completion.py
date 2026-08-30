import math

import numpy as np

from no_free_collapse.completion import (
    completion_dual_value,
    completion_matrix,
    completion_primal_trace,
    offdiag_cube_extrema,
    offdiag_cube_max,
    six_variable_range_hafnian_bound,
)
from no_free_collapse.hafnian_bounds import (
    disjoint_pair_quadratic,
    hafnian,
    twice_offdiag,
)


def _pair_dual(size: int) -> np.ndarray:
    Y = np.zeros((size, size), dtype=np.float64)
    for i in range(0, size, 2):
        Y[i : i + 2, i : i + 2] = np.asarray([[1.0, -1.0], [-1.0, 1.0]])
    return Y


def _conference_witness() -> np.ndarray:
    signs = [
        1, 1, -1, -1, 1,
        1, 1, -1, -1,
        -1, 1, -1,
        -1, -1,
        -1,
    ]
    S = np.zeros((6, 6), dtype=np.float64)
    k = 0
    for i in range(6):
        for j in range(i + 1, 6):
            S[i, j] = S[j, i] = signs[k]
            k += 1
    return S


def test_pair_point_has_exact_minimum_completion_trace_one_half():
    B = disjoint_pair_quadratic(6)
    C = twice_offdiag(B)
    primal = completion_primal_trace(C, np.diag(B))
    dual = completion_dual_value(C, _pair_dual(6))
    assert math.isclose(primal, 0.5, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(dual, primal, rel_tol=0.0, abs_tol=1e-12)

    minimum, maximum = offdiag_cube_extrema(C)
    assert math.isclose(minimum, -0.5, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(maximum, 0.5, rel_tol=0.0, abs_tol=1e-12)
    range_bound = six_variable_range_hafnian_bound(C)
    value = abs(hafnian(C))
    assert math.isclose(range_bound, 1.0 / 192.0, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(value / range_bound, 8.0 / 9.0, rel_tol=0.0, abs_tol=1e-12)


def test_rank_one_point_has_exact_minimum_completion_trace_one_sixth():
    C = np.full((6, 6), 1.0 / 18.0, dtype=np.float64)
    np.fill_diagonal(C, 0.0)
    diagonal = np.full(6, 1.0 / 36.0)
    J = np.ones((6, 6), dtype=np.float64)
    Y = (6.0 / 5.0) * (np.eye(6) - J / 6.0)
    primal = completion_primal_trace(C, diagonal)
    dual = completion_dual_value(C, Y)
    assert math.isclose(primal, 1.0 / 6.0, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(dual, primal, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(offdiag_cube_max(C), 5.0 / 6.0, rel_tol=0.0, abs_tol=1e-12)

    minimum, maximum = offdiag_cube_extrema(C)
    assert math.isclose(minimum, -1.0 / 6.0, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(maximum, 5.0 / 6.0, rel_tol=0.0, abs_tol=1e-12)
    range_bound = six_variable_range_hafnian_bound(C)
    value = abs(hafnian(C))
    assert math.isclose(range_bound, 5.0 / 1728.0, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(value / range_bound, 8.0 / 9.0, rel_tol=0.0, abs_tol=1e-12)


def test_non_psd_one_over_200_witness_needs_larger_trace_completion():
    S = _conference_witness()
    assert np.allclose(S @ S, 5.0 * np.eye(6), atol=1e-12, rtol=0.0)

    C = S / 10.0
    d = math.sqrt(5.0) / 20.0
    diagonal = np.full(6, d)
    A = 0.5 * C
    Y = np.eye(6) - A / d

    primal = completion_primal_trace(C, diagonal)
    dual = completion_dual_value(C, Y)
    expected_tau = 3.0 * math.sqrt(5.0) / 10.0
    assert math.isclose(primal, expected_tau, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(dual, expected_tau, rel_tol=0.0, abs_tol=1e-12)

    # The equal-diagonal trace-1/2 representation used by the bounded
    # degree-two relaxation is not PSD.
    relaxed = completion_matrix(C, np.full(6, 1.0 / 12.0))
    assert float(np.linalg.eigvalsh(relaxed)[0]) < -1e-3

    minimum, s = offdiag_cube_extrema(C)
    value = abs(hafnian(C))
    assert math.isclose(minimum, -0.5, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(s, 0.5, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(value, 1.0 / 200.0, rel_tol=0.0, abs_tol=1e-12)
    range_bound = six_variable_range_hafnian_bound(C)
    assert math.isclose(range_bound, 1.0 / 192.0, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(value / range_bound, 24.0 / 25.0, rel_tol=0.0, abs_tol=1e-12)

    # Once the exact PSD completion cost tau(C) is inserted, the witness is
    # comfortably below the proposed six-variable inequality.
    ratio = 54.0 * value / (primal * s * (primal + s))
    assert ratio < 0.69

    # The generic range-only theorem is already sufficient once the exact
    # completion penalty is included.
    assert 54.0 * range_bound < primal * s * (primal + s)
