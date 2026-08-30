import numpy as np

from no_free_collapse.projection_gradient import six_variable_gradient_energies


def _block_projection(parts):
    p = np.zeros((6, 6), dtype=np.float64)
    offset = 0
    for size in parts:
        p[offset : offset + size, offset : offset + size] = 1.0 / size
        offset += size
    return p


def _random_two_direction_projection(parts, seed):
    rng = np.random.default_rng(seed)
    p = np.zeros((6, 6), dtype=np.float64)
    offset = 0
    for size in parts:
        weights = rng.random(size)
        weights /= np.linalg.norm(weights)
        block = np.outer(weights, weights)
        p[offset : offset + size, offset : offset + size] = block
        offset += size
    return p


def test_exact_partition_constants_at_equal_weights():
    expected = {
        (1, 5): 9.0 / 50.0,
        (2, 4): 33.0 / 160.0,
        (3, 3): 1.0 / 6.0,
    }
    for parts, ratio in expected.items():
        p = _block_projection(parts)
        q1, q2 = six_variable_gradient_energies(p)
        assert np.isclose(q2 / q1, ratio, atol=1e-12, rtol=0.0)


def test_two_plus_four_extremizer_has_exact_energies():
    p = _block_projection((2, 4))
    q1, q2 = six_variable_gradient_energies(p)
    assert np.isclose(q1, 5.0 / 8.0, atol=1e-12, rtol=0.0)
    assert np.isclose(q2, 33.0 / 256.0, atol=1e-12, rtol=0.0)
    assert np.isclose(q2, (33.0 / 160.0) * q1, atol=1e-12, rtol=0.0)


def test_random_two_direction_projections_respect_sharp_global_constant():
    bound = 33.0 / 160.0
    for parts in ((1, 5), (2, 4), (3, 3)):
        for seed in range(100):
            p = _random_two_direction_projection(parts, seed)
            q1, q2 = six_variable_gradient_energies(p)
            if q1 > 1e-14:
                assert q2 <= bound * q1 + 2e-12
