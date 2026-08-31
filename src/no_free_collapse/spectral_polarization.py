"""Spectral-chain polarization for the six-variable hafnian-gradient defect.

This module records an exact reduction from a normalized PSD contraction to a
four-linear defect evaluated on nested spectral projections.  The mixed
nested-projection inequality itself remains conjectural; the identities here
are exact and CPU-verifiable.
"""

from __future__ import annotations

import itertools

import numpy as np

from .projection_gradient import SIX_VARIABLE_EDGES, six_variable_gradient_energies


FOUR_SUBSETS: tuple[tuple[int, int, int, int], ...] = tuple(
    itertools.combinations(range(6), 4)
)


def _validate_six_symmetric(matrix: np.ndarray, *, atol: float = 1e-10) -> np.ndarray:
    A = np.asarray(matrix, dtype=np.float64)
    if A.shape != (6, 6):
        raise ValueError("matrix must have shape (6, 6)")
    if not np.allclose(A, A.T, atol=atol, rtol=atol):
        raise ValueError("matrix must be symmetric")
    return A


def _validate_projection(projection: np.ndarray, *, atol: float = 1e-10) -> np.ndarray:
    P = _validate_six_symmetric(projection, atol=atol)
    if not np.allclose(P @ P, P, atol=atol, rtol=atol):
        raise ValueError("matrix must be an orthogonal projection")
    return P


def zeon_degree_four_product(
    left: np.ndarray,
    right: np.ndarray,
    *,
    atol: float = 1e-10,
) -> np.ndarray:
    """Return the 15 degree-four coefficients of ``Omega_left Omega_right``.

    With commuting square-free generators ``z_i^2=0``, put

        Omega_A = sum_{i<j} A_ij z_i z_j.

    For a four-set ``{a,b,c,d}``, the returned coefficient is the sum over the
    six ordered ways to choose one matching edge from ``left`` and the
    complementary edge from ``right``.  In particular,

        ||zeon_degree_four_product(A,A)||^2 = 4 q2(A).
    """
    A = _validate_six_symmetric(left, atol=atol)
    B = _validate_six_symmetric(right, atol=atol)
    values: list[float] = []
    for a, b, c, d in FOUR_SUBSETS:
        values.append(
            float(
                A[a, b] * B[c, d]
                + A[c, d] * B[a, b]
                + A[a, c] * B[b, d]
                + A[b, d] * B[a, c]
                + A[a, d] * B[b, c]
                + A[b, c] * B[a, d]
            )
        )
    return np.asarray(values, dtype=np.float64)


def _offdiag_inner(left: np.ndarray, right: np.ndarray) -> float:
    return float(sum(left[i, j] * right[i, j] for i, j in SIX_VARIABLE_EDGES))


def nested_projection_polarized_defect(
    projections: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | list[np.ndarray],
    *,
    atol: float = 1e-10,
    check_nested: bool = True,
) -> float:
    """Evaluate the symmetric four-linear mixed projection defect.

    For four projections ``P1,...,P4`` define

        D(P1,P2,P3,P4)
          = (1/24) sum_{a<b} <off(Pa),off(Pb)>
            - (1/12) sum_pairings
                <Omega_Pa Omega_Pb, Omega_Pc Omega_Pd>.

    Setting all four arguments equal gives exactly ``q1(P)/4-q2(P)``.

    When ``check_nested=True`` this routine also verifies

        P1 <= P2 <= P3 <= P4

    in projection order.  Global nonnegativity on all such nested quadruples
    is the remaining mixed inequality that would imply the PSD-contraction
    theorem; this function evaluates the exact quantity but does not assert
    that conjecture.
    """
    if len(projections) != 4:
        raise ValueError("exactly four projections are required")
    Ps = [_validate_projection(P, atol=atol) for P in projections]

    if check_nested:
        for i in range(4):
            for j in range(i + 1, 4):
                if not np.allclose(Ps[i] @ Ps[j], Ps[i], atol=20 * atol, rtol=20 * atol):
                    raise ValueError("projections must be nested in nondecreasing order")

    q1_part = sum(_offdiag_inner(Ps[i], Ps[j]) for i, j in itertools.combinations(range(4), 2)) / 24.0

    products = {
        (i, j): zeon_degree_four_product(Ps[i], Ps[j], atol=atol)
        for i in range(4)
        for j in range(i, 4)
    }
    pairings = (
        ((0, 1), (2, 3)),
        ((0, 2), (1, 3)),
        ((0, 3), (1, 2)),
    )
    q2_part = 0.0
    for (i, j), (k, l) in pairings:
        key1 = (min(i, j), max(i, j))
        key2 = (min(k, l), max(k, l))
        q2_part += float(products[key1] @ products[key2])
    q2_part /= 12.0
    return float(q1_part - q2_part)


def normalized_spectral_chain(
    contraction: np.ndarray,
    *,
    atol: float = 1e-10,
) -> tuple[float, np.ndarray, tuple[np.ndarray, ...]]:
    """Return ``(lambda_max, delta, projections)`` for a PSD contraction.

    If ``A`` is nonzero with ``0 <= A <= I``, normalize ``B=A/lambda_max``.
    Sorting the eigenvalues of ``B`` decreasingly gives

        B = sum_{k=1}^6 delta_k P_k,
        delta_k >= 0,
        sum delta_k = 1,

    where ``P_k`` is the projection onto the first ``k`` eigenvectors and the
    ``P_k`` form a nested chain.  Repeated eigenvalues merely create zero
    ``delta_k``, so the basis choice inside a repeated eigenspace is harmless.
    """
    A = _validate_six_symmetric(contraction, atol=atol)
    eigenvalues, eigenvectors = np.linalg.eigh(A)
    if float(np.min(eigenvalues)) < -20.0 * atol:
        raise ValueError("matrix must be positive semidefinite")
    if float(np.max(eigenvalues)) > 1.0 + 20.0 * atol:
        raise ValueError("matrix must satisfy A <= I")

    eigenvalues = np.clip(eigenvalues, 0.0, 1.0)
    lambda_max = float(np.max(eigenvalues))
    if lambda_max <= atol:
        return 0.0, np.zeros(6, dtype=np.float64), tuple(np.zeros((6, 6)) for _ in range(6))

    order = np.argsort(eigenvalues)[::-1]
    lam = eigenvalues[order] / lambda_max
    U = eigenvectors[:, order]
    extended = np.concatenate([lam, np.asarray([0.0])])
    delta = extended[:-1] - extended[1:]
    delta[np.abs(delta) < 100.0 * atol] = 0.0
    if float(np.min(delta)) < -100.0 * atol:
        raise ValueError("internal spectral weights are not nonnegative")

    projections: list[np.ndarray] = []
    for k in range(1, 7):
        P = U[:, :k] @ U[:, :k].T
        projections.append(P)

    if not np.isclose(float(np.sum(delta)), 1.0, atol=100.0 * atol, rtol=0.0):
        raise ValueError("internal spectral weights must sum to one")
    return lambda_max, np.asarray(delta, dtype=np.float64), tuple(projections)


def spectral_chain_polarized_average(
    contraction: np.ndarray,
    *,
    atol: float = 1e-10,
) -> tuple[float, float]:
    """Return ``(direct_defect, polarized_average)`` for normalized ``A``.

    For nonzero ``A`` this routine first normalizes ``B=A/lambda_max`` and
    writes ``B=sum delta_k P_k`` using :func:`normalized_spectral_chain`.
    The exact polarization identity is

        q1(B)/4-q2(B)
          = sum_{i,j,k,l} delta_i delta_j delta_k delta_l
              D(P_i,P_j,P_k,P_l).

    The second returned number is the right-hand side evaluated directly.
    Equality of the two outputs is an identity; nonnegativity of every mixed
    nested term ``D`` remains conjectural.
    """
    lambda_max, delta, projections = normalized_spectral_chain(contraction, atol=atol)
    if lambda_max == 0.0:
        return 0.0, 0.0

    B = np.asarray(contraction, dtype=np.float64) / lambda_max
    q1, q2 = six_variable_gradient_energies(B)
    direct = 0.25 * q1 - q2

    total = 0.0
    active = [i for i, weight in enumerate(delta) if weight > atol]
    cache: dict[tuple[int, int, int, int], float] = {}
    for i in active:
        for j in active:
            for k in active:
                for l in active:
                    indices = tuple(sorted((i, j, k, l)))
                    if indices not in cache:
                        Ps = [projections[index] for index in indices]
                        cache[indices] = nested_projection_polarized_defect(Ps, atol=atol)
                    weight = float(delta[i] * delta[j] * delta[k] * delta[l])
                    total += weight * cache[indices]
    return float(direct), float(total)
