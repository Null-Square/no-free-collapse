# Six-variable hafnian gradient structure

**Current status:** the full real PSD-gradient theorem is now proved, not only
its projection cases. See [`matching_dilation.md`](matching_dilation.md) for the
complementary-subset dilation, the sharper spectral-diameter inequality, and
all PSD equality cases. The separate cube-normalized `1/216` problem remains open.

## Definitions

For a real symmetric `6 x 6` matrix A, index edges by `e={i,j}` and put

\[
a_e=A_{ij},\qquad h_e=\operatorname{haf}(A_{\widehat i,\widehat j}),
\quad q_1(A)=\sum_ea_e^2,\quad q_2(A)=\sum_eh_e^2.
\]

The 15-by-15 matching operator T(A) has entry zero for intersecting edges e,f;
otherwise its entry is `a_g`, where g is the edge complementary to e union f.
Direct matching enumeration gives

\[
T(A)a=2h,\qquad4q_2(A)=\|T(A)a\|^2,\qquad3\operatorname{haf}(A)=a^Th.
\]

The new dilation proves

\[
\|T(A)\|_{op}\le\lambda_{max}(A)-\lambda_{min}(A).
\]

Consequently

\[
q_2(A)\le\frac{(\lambda_{max}(A)-\lambda_{min}(A))^2}{4}q_1(A).
\]

For PSD A this implies `q2(A)<=lambda_max(A)^2 q1(A)/4`; for `0<=A<=I`
it implies `q2(A)<=q1(A)/4`. No conjectural mixed-kernel positivity is used.

## Corrected operator warning

Earlier notes rejected `||T(A)||<=lambda_max(A)` for PSD A. That rejection was
incorrect. The conference counterexample concerns an involution K with
`||T(K)||=3/sqrt(5)>1`; its PSD projection is `P=(I+K)/2`, so
`||T(P)||=3/(2 sqrt(5))<1`. The PSD operator statement is now proved by an
explicit pair of PSD compression identities. The regression test in
`test_matching_dilation.py` keeps these two scalings distinct.

## Square-free interpretation

For commuting generators satisfying `z_i^2=0`, set
`Omega_A=sum_{i<j} A_ij z_i z_j`. In the Euclidean coefficient norm,
`||Omega_A||^2=q1(A)` and `||Omega_A^2||^2=4q2(A)`.
This is a reformulation in established zeon algebra machinery, not a claim to
invent zeon algebras or their norm theory.

## Rank-one proof and its degenerate equality cases

Let `P=uu^T`, `sum u_i^2=1`, and `d_i=u_i^2`. Then

\[
q_1(P)=e_2(d),\qquad q_2(P)=9e_4(d).
\]

Maclaurin's inequalities and Cauchy-Schwarz give

\[
e_4(d)\le\frac{e_2(d)^2}{15},\qquad
 e_2(d)=\frac{1-\sum_i d_i^2}{2}\le\frac5{12}.
\]

Therefore `q2<=3 e2^2/5<=e2/4=q1/4`. If `e2=0`, u is a signed coordinate
vector and `q1=q2=0`. If `e2>0`, equality forces `e2=5/12`, hence all
`d_i=1/6`. These are exactly the equality cases. Nonzero rank of P alone does
not exclude the coordinate-projection case.

## Historical rank-specific proofs

Complementation preserves q1 and q2 because off-diagonal entries change sign
and the four-variable hafnian has degree two. Rank five follows from rank one,
and rank four from rank two. Ranks zero and six are trivial. The independent
middle-rank arguments are in
[`rank_two_global_projection.md`](rank_two_global_projection.md) and
[`rank_three_global_gradient.md`](rank_three_global_gradient.md).
Their historical statements of what remained open are superseded by
[`RESULTS.md`](RESULTS.md).

For rank-three P, `K=2P-I` is a trace-zero symmetric involution. Writing
`y_i=P_ii-1/2` and `x_ij=P_ij^2` gives
`sum_i y_i=0` and `sum_{j!=i}x_ij=1/4-y_i^2`. The exact defect involving wedge
spreading and the commutator energy is recorded in
[`rank_three_defect_identity.md`](rank_three_defect_identity.md).

## Equality in the completed global theorem

For nonzero PSD A with L=lambda_max(A), nondegenerate equality occurs exactly
when A/L is `ss^T/6`, `I-ss^T/6`, or `(I+M)/2`, with s a sign vector and M a
signed perfect-matching matrix. All diagonal PSD matrices give the degenerate
scalar equality. The complete necessity and sufficiency proof is in the
[dilation note](matching_dilation.md).

Complex zero-diagonal spectral bounds and quantitative phase-matching stability
are separately proved in [`zero_diagonal_spectral.md`](zero_diagonal_spectral.md).
Those stability statements do not classify near-equality in all three real PSD
families; the hypotheses and conclusions must remain distinct.
