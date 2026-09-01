# Results ledger

This file is the reviewer-facing index of mathematical claims in the repository. It is designed to answer four questions quickly:

1. What exactly is proved?
2. Where is the proof written?
3. Which tests protect the supporting identities/certificates?
4. Which nearby statements are still conjectural?

## Status legend

- **PROVED** — a complete analytic proof is recorded in `docs/`.
- **SOLVED EXACTLY** — an extremal class is optimized with a matching construction and upper bound.
- **PROVED REDUCTION** — an exact equivalence/reduction is proved, but the reduced target may remain open.
- **NUMERICAL EVIDENCE** — diagnostic only; never used as a theorem.
- **OPEN** — explicitly not claimed.

Tests verify algebraic identities, constructions, exact constants, rational certificates, and representative random instances. Tests are not presented as substitutes for analytic proofs.

## A. Interaction order, normalization, and conditioning

| ID | Claim | Status | Proof note | Main tests |
| --- | --- | --- | --- | --- |
| A1 | Fixed-norm order-`r` polynomial preparation has Born interaction degree at most `2r`. | **PROVED** | [`math.md`](math.md) | [`test_degree_bound.py`](../tests/test_degree_bound.py), [`test_interactions.py`](../tests/test_interactions.py) |
| A2 | The factor-of-two degree ceiling is tight via a fixed-norm parity construction. | **PROVED** | [`math.md`](math.md) | [`test_degree_bound.py`](../tests/test_degree_bound.py) |
| A3 | Input-dependent normalization can create arbitrarily high exact interaction degree from low-order amplitudes. | **PROVED** | [`math.md`](math.md) | [`test_normalization_loophole.py`](../tests/test_normalization_loophole.py) |
| A4 | Conditioning controls high-order Walsh leakage exponentially, with `rho=(sqrt(kappa)-1)/(sqrt(kappa)+1)`. | **PROVED** | [`conditioning.md`](conditioning.md), [`chebyshev.md`](chebyshev.md) | [`test_conditioning.py`](../tests/test_conditioning.py), [`test_chebyshev.py`](../tests/test_chebyshev.py) |
| A5 | Every order-`r` model has an exact Gram form; the minimum latent dimension for fixed `Q` is `rank(Q)` and the optimal linear collapse is spectral. | **PROVED** | [`gram.md`](gram.md) | [`test_gram.py`](../tests/test_gram.py) |

## B. Exact solvable classes

| ID | Claim | Status | Proof note | Main tests |
| --- | --- | --- | --- | --- |
| B1 | The full real order-1 class invariant under coordinate permutations and global sign reversal is optimized exactly by the mean-field normalizer. | **SOLVED EXACTLY** | [`symmetric_optimality.md`](symmetric_optimality.md) | [`test_symmetric.py`](../tests/test_symmetric.py) |
| B2 | The symmetric optimum is not globally optimal; a matched-pair symmetry-breaking construction beats it. | **PROVED** | [`paired.md`](paired.md) | [`test_paired.py`](../tests/test_paired.py) |
| B3 | The full matched pair-block class is optimized exactly, including its rank and Gamma-function capacity formula. | **SOLVED EXACTLY** | [`paired.md`](paired.md) | [`test_paired.py`](../tests/test_paired.py) |

## C. Hafnian / low-conditioning regime

| ID | Claim | Status | Proof note | Main tests |
| --- | --- | --- | --- | --- |
| C1 | The first full-order quadratic-normalizer term is a hafnian and obeys the universal minimax/Chebyshev bound. | **PROVED** | [`hafnian_bound.md`](hafnian_bound.md) | [`test_hafnian_bounds.py`](../tests/test_hafnian_bounds.py) |
| C2 | The four-variable extremal hafnian problem is solved globally with optimum `1/16`. | **SOLVED EXACTLY** | [`hafnian_bound.md`](hafnian_bound.md) | [`test_hafnian_bounds.py`](../tests/test_hafnian_bounds.py) |
| C3 | The equal disjoint-pair six-variable point is an unrestricted first-order PSD local optimum. | **PROVED** | [`pair_local_optimality.md`](pair_local_optimality.md) | [`test_pair_local_optimality.py`](../tests/test_pair_local_optimality.py) |
| C4 | Every nonzero zero-slope tangent at that six-variable point has strictly negative second variation; zero-slope feasible line segments are globally non-improving. | **PROVED** | [`six_variable_second_order.md`](six_variable_second_order.md) | [`test_six_variable_second_order.py`](../tests/test_six_variable_second_order.py) |
| C5 | The candidate six-variable inequality `54|haf(C)| <= t s (t+s)` holds on the entire rank-one PSD stratum. | **PROVED** | [`six_variable_rank_one.md`](six_variable_rank_one.md) | [`test_six_variable_rank_one.py`](../tests/test_six_variable_rank_one.py) |
| C6 | The six diagonal variables can be eliminated exactly by the minimum-trace PSD completion `tau(C)`, with elliptope dual. | **PROVED REDUCTION** | [`diagonal_completion.md`](diagonal_completion.md) | [`test_diagonal_completion.py`](../tests/test_diagonal_completion.py) |
| C7 | The range-only bound `|haf(C)| <= a s(a+s)/48` holds without a PSD assumption and reduces any PSD counterexample to a thin `tau≈a` shell. | **PROVED** | [`six_variable_range_bound.md`](six_variable_range_bound.md) | see theorem-specific tests in `tests/` and exact moment identities in the proof note |

## D. Six-variable projection-gradient theorem

Define

\[
q_1(A)=\sum_{i<j}A_{ij}^2,
\qquad
q_2(A)=\sum_{i<j}\operatorname{haf}(A_{\widehat i,\widehat j})^2.
\]

The completed projection theorem is

\[
\boxed{q_2(P)\le\frac14q_1(P)}
\]

for every real `6 x 6` orthogonal projection `P`.

| ID | Rank / component | Status | Proof note | Main tests |
| --- | --- | --- | --- | --- |
| D1 | Perfect-matching operator and zeon reformulation; rank-one theorem. | **PROVED** | [`projection_gradient.md`](projection_gradient.md) | [`test_projection_gradient.py`](../tests/test_projection_gradient.py) |
| D2 | Rank five by projection complementation. | **PROVED** | [`projection_gradient.md`](projection_gradient.md) | [`test_projection_gradient.py`](../tests/test_projection_gradient.py) |
| D3 | Rank-three equal-diagonal fourth-power identity and equality geometry. | **PROVED** | [`rank_three_equal_diagonal.md`](rank_three_equal_diagonal.md) | [`test_rank_three_equal_diagonal.py`](../tests/test_rank_three_equal_diagonal.py) |
| D4 | Exact rank-three defect identity. | **PROVED** | [`rank_three_defect_identity.md`](rank_three_defect_identity.md) | [`test_rank_three_defect_identity.py`](../tests/test_rank_three_defect_identity.py) |
| D5 | Global rank-three projection theorem. | **PROVED** | [`rank_three_global_gradient.md`](rank_three_global_gradient.md) | [`test_rank_three_global_gradient.py`](../tests/test_rank_three_global_gradient.py) |
| D6 | Rank-two Pluecker identity and two-direction sharp theorem. | **PROVED** | [`rank_two_plucker.md`](rank_two_plucker.md), [`rank_two_two_direction.md`](rank_two_two_direction.md) | [`test_rank_two_plucker.py`](../tests/test_rank_two_plucker.py), [`test_rank_two_two_direction.py`](../tests/test_rank_two_two_direction.py) |
| D7 | Arbitrary rank-two balanced-diagonal region `max_i d_i<=1/2`. | **PROVED** | [`rank_two_balanced_diagonal.md`](rank_two_balanced_diagonal.md) | [`test_rank_two_balanced_diagonal.py`](../tests/test_rank_two_balanced_diagonal.py) |
| D8 | Rank-two high-pair region via diagonal quadratic dual. | **PROVED** | [`rank_two_dual_high_pair.md`](rank_two_dual_high_pair.md) | [`test_rank_two_dual_high_pair.py`](../tests/test_rank_two_dual_high_pair.py) |
| D9 | Remaining rank-two middle strip; hence all rank-two projections. | **PROVED** | [`rank_two_global_projection.md`](rank_two_global_projection.md) | [`test_rank_two_global_projection.py`](../tests/test_rank_two_global_projection.py) |
| D10 | Rank four by complementation; ranks zero/six trivial. | **PROVED** | [`rank_two_global_projection.md`](rank_two_global_projection.md), [`projection_gradient.md`](projection_gradient.md) | [`test_rank_two_global_projection.py`](../tests/test_rank_two_global_projection.py), [`test_projection_gradient.py`](../tests/test_projection_gradient.py) |
| D11 | All projection ranks `0,...,6`. | **PROVED** | dependency of D1–D10 | full `pytest` suite |

### Equality / sharpness notes

The projection theorem is a `1/4` contraction theorem. Some rank strata have sharper internal constants. In particular, the two-direction rank-two class has sharp ratio `33/160`, attained by the equal-weight `2+4` block geometry. The repository does **not** currently claim `33/160` as the sharp global rank-two constant.

## E. PSD-contraction frontier

| ID | Claim | Status | Proof note | Main tests |
| --- | --- | --- | --- | --- |
| E1 | A normalized PSD contraction has an exact spectral-chain representation by nested projections with nonnegative barycentric increments. | **PROVED** | [`spectral_polarization.md`](spectral_polarization.md) | [`test_spectral_polarization.py`](../tests/test_spectral_polarization.py) |
| E2 | The mixed four-slot homogenized kernel reconstructs `q1/4-q2` exactly from the nested spectral chain. | **PROVED REDUCTION** | [`spectral_polarization.md`](spectral_polarization.md) | [`test_spectral_polarization.py`](../tests/test_spectral_polarization.py) |
| E3 | Mixed nested rank patterns `(1,5,5,5)` and `(1,1,1,5)` have nonnegative kernel. | **PROVED** | [`mixed_rank_1555.md`](mixed_rank_1555.md) | [`test_mixed_rank_1555.py`](../tests/test_mixed_rank_1555.py) |
| E4 | `q2(A)<=q1(A)/4` for every `0<=A<=I`. | **OPEN** | reduced to nested mixed kernel positivity by E1–E2 | diagnostics only; no theorem claim |
| E5 | `q2(A)<=lambda_max(A)^2 q1(A)/4` for every PSD `A`. | **OPEN** | would follow from E4 by scaling | — |

Important terminology: the spectral kernel is symmetric in four slots after homogenizing the quadratic term with `sum delta_k=1`; it is **not** literally matrix-four-linear. This prevents invalid rank-one termwise expansion arguments.

## F. Open extremal problems

| ID | Statement | Status |
| --- | --- | --- |
| F1 | Full six-variable sharp PSD hafnian inequality `54|haf(C)| <= tau(C)s(C)(tau(C)+s(C))`. | **OPEN globally**; proved on substantial strata and outside the thin completion shell. |
| F2 | General even-`m` sharp hafnian optimum `m^{-m/2}` for the equal pair construction. | **OPEN**. |
| F3 | Sharp global rank-two projection constant equals `33/160`. | **OPEN strengthening**; `1/4` theorem is proved globally. |

## Numerical diagnostics and failed stronger statements

The project keeps failed proof routes visible in the relevant notes when they materially constrain future work. In particular:

- the stronger operator norm claim `||T(A)||_op <= lambda_max(A)` is false;
- separate eigenvalue convexity does not reduce PSD contractions to projections;
- Bernoulli rounding of eigenvectors does not preserve the defect in the required direction;
- atom-by-atom positivity of the mixed spectral kernel is false;
- coefficientwise Bernstein positivity along an arbitrary interpolation is not assumed unless proved.

These failures are part of the robustness story: the main theorem statements above do not rely on them.

## Paper baseline

For the first manuscript, the natural cutoff is:

- all results A1–E3 as completed material;
- E4–E5 and F1–F3 as open problems / future work;
- exploratory numerical observations included only when clearly labeled and useful for motivation.

This makes the paper complete without requiring the unresolved full PSD-contraction theorem.
