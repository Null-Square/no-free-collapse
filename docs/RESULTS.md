# Results ledger

**Updated September 12, 2026.** This is the current authority for theorem status.
Older progress notes retain their historical derivations; statements there that
rank two, rank three, or the full PSD-gradient bound are still open are superseded
by this ledger, [`matching_dilation.md`](matching_dilation.md), and
[`subhafnian_hierarchy.md`](subhafnian_hierarchy.md).

## Main results of the publication branch

For real symmetric `A` of order `2m`, `m>=2`, the complementary-edge weighted
Kneser operator satisfies the sharp bound

\[
\|T_m(A)\|_{op}\le\frac{m+1}{4}
(\lambda_{max}(A)-\lambda_{min}(A)).
\]

The proof is an explicit difference of two complementary-subset compressions,
with two PSD slack certificates. This operator theorem now yields an all-orders
spectral hierarchy for averaged principal subhafnian energies. With

\[
H_{2r}(A)=\sum_{|S|=2r}|haf(A[S])|^2,\qquad
\overline H_{2r}(A)=\frac{H_{2r}(A)}{\binom N{2r}},
\]

one has, for every real symmetric `N x N` matrix with `N>=4r-2`,

\[
\boxed{\overline H_{2r}(A)\le\frac{d(A)^2}{4}\overline H_{2r-2}(A)},
\qquad d(A)=\lambda_{max}(A)-\lambda_{min}(A).
\]

The factor `1/4` is sharp at every order in the uniform admissible family.
At the critical size `N=4r-2`, balanced signed rank-one projections, their
rank-`N-1` complements, and signed perfect-matching half-rank projections all
attain equality. No complete equality classification is claimed for `r>=3`.

In the normalized energy

\[
G_r=\overline H_{2r}/((2r-1)!!)^2,
\]

Roos's published log-subadditivity gives `G_r<=G_{r-1}G_1`, while the new
spectral recurrence gives

\[
G_r\le\frac{d(A)^2}{4(2r-1)^2}G_{r-1}.
\]

Taking the smaller factor at each order yields a **one-switch hybrid bound**.
The scale-free parameter `theta=4G_1/d(A)^2` satisfies the sharp universal
bound `theta<=1/(N-1)`, with equality exactly (after shift and scale) at
rank-`N/2` constant-diagonal `1/2` projections. At critical matching
projections, the new spectral step improves the Roos one-step coefficient by
`(2r-1)^2/(4r-3)=Theta(r)`.

In order six the operator theorem also proves

\[
q_2(A)\le\frac{d(A)^2}{4}q_1(A),\qquad
|haf(A)|\le\frac{d(A)}6q_1(A),
\]

so the full PSD-contraction theorem and its scaled PSD version, formerly E4
and E5, are **PROVED**. All nondegenerate PSD equality cases are classified:
scaled signed balanced rank-one projections, their rank-five complements, and
signed pair-block rank-three projections. Diagonal PSD matrices are the
degenerate scalar-bound equality cases.

For a real symmetric pure Gaussian-boson-sampling matrix, the averaged
collision-free per-pattern probabilities inherit the same spectral hierarchy,
because the common pattern-independent normalization cancels. This is an
unconditional collision-free sector statement, not a claim about threshold
Torontonians or collisionful outcomes.

The separate Boolean-cube-normalized `1/216` optimum remains **OPEN**.

## Status legend

**PROVED** means a complete analytic argument is recorded. **SOLVED EXACTLY**
means a constrained extremal class has a matching upper bound and construction.
**PROVED REDUCTION** does not mean its reduced target is solved.
**NUMERICAL EVIDENCE** is diagnostic only. **OPEN** means no theorem is claimed.
Tests check identities and certificates; they do not replace analytic proofs.

## A. Interaction order, normalization, and conditioning

| ID | Claim | Status | Proof |
| --- | --- | --- | --- |
| A1 | Fixed-norm order-r preparation has degree at most 2r. | PROVED; standard degree mechanism | [math](math.md) |
| A2 | The degree factor two is attained by a parity construction. | PROVED | [math](math.md) |
| A3 | Input-dependent normalization can create arbitrarily high exact degree. | PROVED | [math](math.md) |
| A4 | Conditioning controls high-order Walsh leakage exponentially. | PROVED | [conditioning](conditioning.md), [Chebyshev](chebyshev.md) |
| A5 | Exact Gram form, minimal dimension rank(Q) for a fixed Gram representation, and spectral optimal readout. | PROVED | [Gram](gram.md) |

The corresponding `test_degree_bound`, `test_interactions`, `test_normalization_loophole`,
`test_conditioning`, `test_chebyshev`, and `test_gram` modules protect the implementations.
A fixed representation's rank is not asserted to minimize dimension across all
representations of the same function. High-degree probabilities do not by themselves
establish high-degree threshold decisions or a computational advantage.

## B. Exact solvable classes

| ID | Claim | Status | Proof |
| --- | --- | --- | --- |
| B1 | Exact permutation- and global-sign-invariant real order-one optimum. | SOLVED EXACTLY in that class | [symmetric optimality](symmetric_optimality.md) |
| B2 | A matched-pair construction beats the symmetric optimum. | PROVED | [paired](paired.md) |
| B3 | Exact matched pair-block capacity, rank, and Gamma-function formula. | SOLVED EXACTLY in that class | [paired](paired.md) |

Tests: `test_symmetric.py`, `test_paired.py`. Neither classwise optimum is
asserted to solve the unrestricted normalized-model capacity problem.

## C. Hafnian and low-conditioning regime

| ID | Claim | Status | Proof |
| --- | --- | --- | --- |
| C1 | Leading full-order term is a hafnian with a minimax/Chebyshev upper bound. | PROVED | [hafnian bound](hafnian_bound.md) |
| C2 | Global four-variable PSD cube-normalized optimum is 1/16. | SOLVED EXACTLY | [hafnian bound](hafnian_bound.md) |
| C3 | Six-variable equal-pair point is an unrestricted first-order PSD local optimum. | PROVED | [local optimality](pair_local_optimality.md) |
| C4 | Strict second variation on nonzero zero-slope tangents; non-improving zero-slope feasible lines. | PROVED | [second order](six_variable_second_order.md) |
| C5 | Candidate 54-hafnian inequality holds on the rank-one PSD stratum. | PROVED | [rank one](six_variable_rank_one.md) |
| C6 | Exact minimum-trace diagonal completion with elliptope dual. | PROVED REDUCTION | [completion](diagonal_completion.md) |
| C7 | Range-only hafnian bound and thin-shell reduction. | PROVED | [range bound](six_variable_range_bound.md) |

These results do not prove the six-variable global cube optimum.

## D. Historical projection-gradient proof components

Every real six-dimensional orthogonal projection satisfies `q2<=q1/4`.
This now follows at once from the global dilation theorem, but the independent
rank-specific arguments remain available.

| ID | Component | Status | Proof |
| --- | --- | --- | --- |
| D1 | Matching operator, square-free interpretation, rank-one theorem. | PROVED | [projection gradient](projection_gradient.md) |
| D2 | Rank five by complementation. | PROVED | [projection gradient](projection_gradient.md) |
| D3 | Rank-three equal-diagonal identity and equality geometry. | PROVED | [equal diagonal](rank_three_equal_diagonal.md) |
| D4 | Exact rank-three defect identity. | PROVED | [defect](rank_three_defect_identity.md) |
| D5 | Global rank-three theorem. | PROVED | [rank three](rank_three_global_gradient.md) |
| D6 | Rank-two Pluecker identity and sharp two-direction class. | PROVED | [Pluecker](rank_two_plucker.md), [two directions](rank_two_two_direction.md) |
| D7 | Rank-two balanced-diagonal region. | PROVED | [balanced diagonal](rank_two_balanced_diagonal.md) |
| D8 | Rank-two high-pair region. | PROVED | [high pair](rank_two_dual_high_pair.md) |
| D9 | Remaining rank-two strip and global rank-two theorem. | PROVED | [rank two](rank_two_global_projection.md) |
| D10 | Rank four by complementation; ranks zero/six trivial. | PROVED | [rank two](rank_two_global_projection.md) |
| D11 | All projection ranks. | PROVED | D1-D10; now also [global dilation](matching_dilation.md) |

Rank-one equality includes coordinate projections (`q1=q2=0`) as well as
balanced signed rank-one projections. Nonzero rank alone does not exclude
degenerate equality. The two-direction rank-two constant 33/160 is not asserted
here to be the sharp global rank-two constant.

## E. Spectral operator, PSD gradient, and subhafnian hierarchy

| ID | Claim | Status | Proof / tests |
| --- | --- | --- | --- |
| E1 | Nested-projection spectral-chain representation. | PROVED | [spectral polarization](spectral_polarization.md), `test_spectral_polarization.py` |
| E2 | Homogenized mixed kernel reconstructs the defect. | PROVED REDUCTION | [spectral polarization](spectral_polarization.md) |
| E3 | Specific nested mixed rank patterns have nonnegative kernel. | PROVED | [mixed ranks](mixed_rank_1555.md), `test_mixed_rank_1555.py` |
| E4 | `q2(A)<=q1(A)/4` for every `0<=A<=I`. | **PROVED** | [dilation](matching_dilation.md), `test_matching_dilation.py` |
| E5 | `q2(A)<=lambda_max(A)^2 q1(A)/4` for every real PSD A. | **PROVED** | [dilation](matching_dilation.md), `test_matching_dilation.py` |
| E6 | Sharp all-even-dimensional weighted Kneser spectral-diameter bound and explicit PSD certificates. | **PROVED** | [dilation](matching_dilation.md), exact integer matrix-basis checks |
| E7 | Complete six-variable PSD equality classification and sharp hafnian-energy inequality. | **PROVED** | [dilation](matching_dilation.md), exact rational equality tests |
| E8 | Critical all-orders subhafnian step `H_{2r}<=d^2 H_{2r-2}/4` in order `4r-2`, sharp for every `r`. | **PROVED** | [hierarchy](subhafnian_hierarchy.md), exact recursion tests |
| E9 | Arbitrary-ambient averaged hierarchy `Hbar_{2r}<=d^2 Hbar_{2r-2}/4` for `N>=4r-2`. | **PROVED** | [hierarchy](subhafnian_hierarchy.md), exact double counting + random diagnostics |
| E10 | Roos/spectral one-switch hybrid and closed form. | **PROVED** | [hierarchy](subhafnian_hierarchy.md), `test_subhafnian_hierarchy.py` |
| E11 | Sharp spectral edge-density bound `theta<=1/(N-1)` and equality class. | **PROVED** | [hierarchy](subhafnian_hierarchy.md), equality tests |
| E12 | Collision-free real GBS averaged photon-sector corollary. | **PROVED consequence** | [hierarchy](subhafnian_hierarchy.md), `e22_subhafnian_hierarchy.py` diagnostics |

E4-E12 bypass, rather than prove, coefficientwise positivity of every mixed
spectral kernel. The earlier warning against invalid rank-one termwise arguments
remains valid.

## F. Open problems retained

| ID | Statement | Status |
| --- | --- | --- |
| F1 | Full six-variable cube-normalized optimum 1/216 and the stronger 54-hafnian completion inequality. | OPEN globally |
| F2 | General even-dimensional cube-normalized optimum m^(-m/2). | OPEN |
| F3 | Sharp global rank-two projection ratio 33/160. | OPEN strengthening |
| F4 | Complete equality classification for the critical subhafnian step when `r>=3`. | OPEN strengthening |
| F5 | Best fixed-ambient-dimension constants for `N>4r-2`. | OPEN strengthening |

## G. Complex zero-diagonal spectral bounds and stability

For complex symmetric zero-diagonal B of order n>=4, with T its off-diagonal
squared energy and H4 its summed four-subhafnian energy,
`H4 <= (n-2)||B||_op^2 T/4`. An exact four-term nonnegative defect proves this,
classifies equality, and gives constructive matching recovery. The coefficient
is sharp in every even order. In order six, `|haf(B)|<=||B||_op^3`, and a
unit contraction with `|haf(B)|>=1-epsilon`, `epsilon<=1/48`, lies within
squared Frobenius distance `34 epsilon` of a phase-weighted perfect matching.
A rational rotation family shows that the linear squared-distance order cannot
be improved uniformly to little-o(epsilon). The numerical constant 34 is not
claimed optimal.

Proof: [zero-diagonal spectral note](zero_diagonal_spectral.md). Tests:
`test_zero_diagonal.py`, `test_spectral_stability_sharpness.py`. These complex
stability results are not a quantitative classification of all three PSD
families in E7.

## Corrected failed-route record

The earlier sentence claiming `||T(A)||<=lambda_max(A)` false for PSD A was
incorrect and is withdrawn. The conference example refutes `||T(K)||<=1` for
certain symmetric involutions K; its PSD projection is `(I+K)/2` and its
matching operator is smaller by a factor two. E6 explicitly proves the PSD
bound. The regression test records both scalings.

Other historical failed routes remain distinct: separate eigenvalue convexity,
Bernoulli eigenvector rounding in the needed direction, atomwise mixed-kernel
positivity, and unproved coefficientwise Bernstein assertions are not used by
the new theorem.

## Publication and verification

The focused primary contribution is now E6-E12, with E4-E5 as the motivating
six-variable consequence and G as companion stability work. The main manuscript
has been expanded to the all-orders hierarchy and its Roos/GBS consequences.
The new research overlay contains 124 targeted tests, all passing locally after
the hierarchy additions. Exact algebraic checks, floating-point diagnostics,
and analytic proofs are reported separately. Repository-wide CI is the final
integration gate for the branch revision.

Independent mathematical review, complete priority comparison, author metadata,
and venue-specific submission preparation remain external publication gates.
No journal submission or merge is performed by this research branch.
