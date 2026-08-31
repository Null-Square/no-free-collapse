import numpy as np

from no_free_collapse.rank_two_dual import (
    evaluate_degree_four_bernstein,
    rank_two_high_pair_bernstein_coefficients,
    rank_two_quadratic_dual_certificate,
)
from no_free_collapse.rank_two_projection import rank_two_projection_direct_defect


def _random_rank_two_projection(seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.normal(size=(6, 2)))
    return Q @ Q.T


def _random_high_pair_diagonal(seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    T = rng.uniform(1e-5, 0.5)
    u = rng.dirichlet(np.ones(2))
    v = rng.dirichlet(np.ones(4))
    d = np.asarray([1.0 - T * u[0], 1.0 - T * u[1], *(T * v)])
    rng.shuffle(d)
    return d


def test_quadratic_dual_certificate_lower_bounds_actual_projection_defect():
    for seed in range(100):
        P = _random_rank_two_projection(seed)
        H = rank_two_quadratic_dual_certificate(np.diag(P))
        assert 8.0 * rank_two_projection_direct_defect(P) >= H - 5e-10


def test_high_pair_certificate_reconstructs_exact_bernstein_polynomial():
    for seed in range(100):
        d = _random_high_pair_diagonal(seed)
        x, beta = rank_two_high_pair_bernstein_coefficients(d)
        H = rank_two_quadratic_dual_certificate(d)
        assert np.isclose(
            evaluate_degree_four_bernstein(x, beta),
            H,
            atol=2e-11,
            rtol=2e-11,
        )


def test_high_pair_bernstein_coefficients_have_proved_positive_bounds():
    for seed in range(100):
        d = _random_high_pair_diagonal(seed)
        _, beta = rank_two_high_pair_bernstein_coefficients(d)
        assert np.isclose(beta[0], 0.0, atol=1e-14, rtol=0.0)
        assert np.isclose(beta[1], 1.0 / 4.0, atol=1e-14, rtol=0.0)
        assert beta[2] >= 25.0 / 96.0 - 2e-12
        assert beta[3] >= 1.0 / 32.0 - 2e-12
        assert beta[4] >= 127.0 / 512.0 - 2e-12


def test_beta4_maclaurin_equality_family():
    # The endpoint bound beta4 >= 127/512 is sharp when both normalized
    # groups are uniform: u=(1/2,1/2), v=(1/4,...,1/4).
    T = 0.37
    d = np.asarray([1.0 - T / 2.0, 1.0 - T / 2.0, T / 4.0, T / 4.0, T / 4.0, T / 4.0])
    _, beta = rank_two_high_pair_bernstein_coefficients(d)
    assert np.isclose(beta[2], 25.0 / 96.0, atol=2e-14, rtol=0.0)
    assert np.isclose(beta[4], 127.0 / 512.0, atol=2e-14, rtol=0.0)


def test_high_pair_certificate_equality_only_at_coordinate_endpoint_numerically():
    coordinate = np.asarray([1.0, 1.0, 0.0, 0.0, 0.0, 0.0])
    assert np.isclose(rank_two_quadratic_dual_certificate(coordinate), 0.0, atol=1e-14)

    for seed in range(100):
        d = _random_high_pair_diagonal(seed)
        assert rank_two_quadratic_dual_certificate(d) > 0.0
