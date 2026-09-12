# Current theorem ledger

**Submission-candidate audit, 12 September 2026.** This ledger and the primary manuscript supersede historical statements of unresolved status in the earlier rank-specific progress notes.

## Primary paper

| ID | Precise claim | Status / proof |
|---|---|---|
| E6 | For real symmetric A of order 2m, m>=2, `||T_m(A)|| <= (m+1)d(A)/4`; the coefficient is sharp in every such dimension. | Proved: [dilation](matching_dilation.md) |
| H1 | For r>=2 at critical size N=4r-2, `H_{2r} <= d(A)^2 H_{2r-2}/4`. | Proved; sharp one-step constant: [hierarchy](subhafnian_hierarchy.md) |
| H2 | For r>=2 and N>=4r-2, the same recurrence holds for binomially averaged principal energies. | Proved by interlacing and exact counting |
| H3 | In Roos normalization, `G_r <= G_{r-1} min(G_1,d^2/[4(2r-1)^2])`; the product has at most one switch. | Proved; odd-N justification in [audit](SUBMISSION_AUDIT.md) |
| H4 | Hybrid products are log-subadditive within the admissible order range. | Proved; no claim of global capacity optimality |
| H5 | `theta=4G_1/d^2 <= 4 floor(N^2/4)/(N^2(N-1))` for non-scalar real symmetric A. | Proved; equality iff the normalized matrix is a constant-diagonal projection of rank floor(N/2) or ceil(N/2) |
| E4/E5 | Full real six-variable PSD-gradient theorem: `q2(A)<=lambda_max(A)^2 q1(A)/4`. | **Proved globally**, not only for projections |
| E7 | Nondegenerate PSD equalities are scaled signed balanced rank-one projections, their complements, and signed matching half-rank projections. Diagonal PSD matrices give degenerate scalar equality. | Proved: [dilation](matching_dilation.md), expanded in manuscript |
| H6 | Averaged unconditional per-pattern collision-free GBS probabilities satisfy the hierarchy for real symmetric, pure, zero-displacement sampling matrices. | Corollary; r>=2 and N>=4r-2; no collisionful or threshold claim |

The general hierarchy has sharp witnesses at each critical size. Complete equality classification at r>=3, sharp fixed-N constants for N>4r-2, and sharpness of the iterated hybrid are **not** asserted.

## Companion result

For complex symmetric zero-diagonal B of order N>=4, the four-subhafnian energy is bounded by `(N-2)||B||^2 T(B)/4`. A four-term defect identity proves equality and constructive matching stability. In order six, `|haf(B)|<=||B||^3` and near equality yields a phase matching within squared Frobenius distance `34 epsilon` under the stated hypotheses. The linear stability order is sharp; 34 is not claimed optimal. See [companion note](zero_diagonal_spectral.md). This is not stability around all three PSD equality families.

## Retained earlier results

The following are research background, not equally weighted headline contributions of the primary paper.

| Group | Completed content | Sources |
|---|---|---|
| A | Fixed-norm degree bound; normalization witness; conditioning-controlled leakage; Gram form and fixed-representation rank | [math](math.md), [Chebyshev](chebyshev.md), [Gram](gram.md) |
| B | Exact symmetric and matched-pair model classes | [symmetric optimality](symmetric_optimality.md), [pairs](paired.md) |
| C | Four-variable cube optimum; six-variable local optimality, rank-one and completion/range reductions | [hafnian](hafnian_bound.md), [local](pair_local_optimality.md), [completion](diagonal_completion.md), [range](six_variable_range_bound.md) |
| D | Independent proofs on every six-variable real projection rank | [projection overview](projection_gradient.md), [rank two](rank_two_global_projection.md), [rank three](rank_three_global_gradient.md) |
| E1-E3 | Spectral-chain representation and selected mixed-rank certificates | [polarization](spectral_polarization.md), [mixed ranks](mixed_rank_1555.md) |

Fixed Gram rank is not minimal rank across all functional representations. Classwise capacity optima are not unrestricted optima. Numerical evidence is not a proof.

## Open questions, outside the submission claims

The global cube-normalized six-variable 1/216 optimum and stronger 54-hafnian completion inequality remain open. The general even-variable cube conjecture and sharp global rank-two ratio 33/160 also remain open. The general PSD theorem bypasses, rather than proves, positivity of every mixed spectral coefficient.

## Correction register

The previously rejected PSD operator bound is true; the conference example applied to a differently scaled involution. Coordinate projections belong to the degenerate rank-one equality locus. The old odd-N Roos citation required an additional parity-free argument, now supplied. The old density ceiling 1/(N-1) remains valid but is improved to (N+1)/N^2 when N is odd.

## Verification and publication status

See [reproducibility](REPRODUCIBILITY.md), [audit](SUBMISSION_AUDIT.md), and [submission status](../submission/README.md). Analytic proof, exact instance checks, floating-point diagnostics, and human approval are tracked separately. No journal submission, acceptance, or exhaustive priority review is implied by this ledger.
