# Six-variable hafnian gradient structure

**Status, September 12, 2026.** The gradient inequality is proved for every
real six-dimensional orthogonal projection: see
[`rank_two_global_projection.md`](rank_two_global_projection.md) and
[`rank_three_global_gradient.md`](rank_three_global_gradient.md).
The extension to arbitrary PSD contractions remains open. The new
[`zero_diagonal_spectral.md`](zero_diagonal_spectral.md) proves it for
half-diagonal PSD contractions, including non-idempotent ones, and gives
complex-symmetric spectral bounds and quantitative matching stability.
These statements have different hypotheses and must not be conflated.

## Definitions and the unrestricted target

For a real symmetric `6 x 6` matrix `A`, index edges by `e={i,j}` and define

\[
a_e=A_{ij},\qquad h_e=\operatorname{haf}(A_{\widehat i,\widehat j}),
\qquad q_1(A)=\sum_ea_e^2,\quad q_2(A)=\sum_eh_e^2.
\]

The unrestricted PSD gradient conjecture is

\[
q_2(A)\le\frac{\lambda_{\max}(A)^2}{4}q_1(A),\qquad A\succeq0.
\tag{G}
\]

Equation (G) is **not proved globally** in this repository.

## Weighted perfect-matching operator

Index rows and columns of `T(A)` by the 15 edges of `K_6`. If edges `e,f`
intersect, put `T(A)_{e,f}=0`; otherwise let `g` be the unique remaining edge
completing a perfect matching and put `T(A)_{e,f}=a_g`.
The matrix is symmetric and

\[
\boxed{T(A)a=2h},\qquad
4q_2(A)=\|T(A)a\|_2^2,\qquad q_1(A)=\|a\|_2^2.
\]

The factor two counts the two orders of each matching on the four vertices
complementary to `e`. Thus (G) is a bound on the specific vector `a` generated
by `A`, not an operator-norm bound on arbitrary edge vectors.
The stronger generic operator-norm routes considered in earlier revisions
have been rejected; see [`RESULTS.md`](RESULTS.md) and the explicit involution
counterexample in [`rank_three_defect_identity.md`](rank_three_defect_identity.md).
They should not be described as available proof routes.

## Square-free interpretation

Let `z_i` commute and satisfy `z_i^2=0`, and set

\[
\Omega_A=\sum_{i<j}A_{ij}z_i z_j.
\]

In the Euclidean coefficient norm,

\[
\|\Omega_A\|^2=q_1(A),\qquad \|\Omega_A^2\|^2=4q_2(A).
\]

The degree-four coefficients of `Omega_A^2` are twice the corresponding
four-variable hafnians. Zeon algebras and their general norm theory are prior
art; this identity is a reformulation, not a claim to invent that machinery.

## Rank-one theorem, including the degenerate equality cases

Let `P=uu^T` with `sum_i u_i^2=1`. Then

\[
\boxed{q_2(P)\le\frac14q_1(P).}
\]

Equality holds exactly in either of the following cases:

1. `u` is a signed coordinate vector, so `P` is a coordinate projection and
   `q1=q2=0`;
2. `u_i^2=1/6` for all six coordinates.

Equivalently, under the additional assumption `q1(P)>0`, equality holds
exactly in the second case. Nonzero rank of `P` alone does not exclude the
first case.

**Proof.** Write `d_i=u_i^2`, so `d_i>=0` and `sum_i d_i=1`. Then

\[
q_1(P)=e_2(d),\qquad q_2(P)=9e_4(d).
\]

Each four-variable hafnian has three equal matching products, giving the
second formula. Maclaurin's inequalities and Cauchy-Schwarz give

\[
e_4(d)\le\frac{e_2(d)^2}{15},\qquad
e_2(d)=\frac{1-\sum_i d_i^2}{2}\le\frac5{12},
\]

hence

\[
q_2(P)\le\frac35e_2(d)^2\le\frac14e_2(d)=\frac14q_1(P).
\]

If `e2=0`, exactly one `d_i` equals one, and equality is degenerate. If `e2>0`,
equality in the final inequality forces `e2=5/12`; equality in Cauchy then
forces all six `d_i=1/6`. These conditions also suffice.

## Complementation and the completed projection theorem

For an orthogonal projection, `P^c=I-P` has opposite off-diagonal entries.
Therefore `q1(P^c)=q1(P)` and `q2(P^c)=q2(P)` because the four-variable hafnian
has degree two. Rank five follows from rank one; rank four follows from rank
two. Ranks zero and six are trivial. The rank-two and rank-three proofs linked
above complete the entire real projection locus.

## Rank-three involution and defect coordinates

For rank-three `P`, set `K=2P-I`. Then

\[
K^T=K,\qquad K^2=I,\qquad\operatorname{tr}K=0,
\]

and `q1(P)=q1(K)/4`, `q2(P)=q2(K)/16`. The rank-three inequality is therefore
`q2(K)<=q1(K)` on this involution locus.
Writing `y_i=P_ii-1/2` and `x_ij=P_ij^2`, one has

\[
\sum_i y_i=0,\qquad \sum_{j\ne i}x_{ij}=\frac14-y_i^2.
\]

The nonnegative energies

\[
W=\sum_i\sum_{j<k\atop j,k\ne i}x_{ij}x_{ik},\qquad
L=\sum_{i<j}x_{ij}(y_i-y_j)^2
=\frac12\|[\operatorname{diag}(y),P]\|_F^2
\]

enter the exact identity proved in
[`rank_three_defect_identity.md`](rank_three_defect_identity.md). Its resulting
stability inequality was subsequently proved in
[`rank_three_global_gradient.md`](rank_three_global_gradient.md); it is no
longer an open rank-three target.

## Prior-art boundary

Hafnians, generalized Laplace expansions, averaged squared subhafnians, and
zeon norm inequalities are established topics. Relevant sources include
B. Roos, *New inequalities for permanents and hafnians and some generalizations*,
Linear and Multilinear Algebra 73 (2025), and T. Lindell and G. S. Staples,
*Norm Inequalities in Zeon Algebras*, Advances in Applied Clifford Algebras 29
(2019). The spectral proof note gives a proposition-level comparison with the
displayed bound in Roos's Theorem 2.3 and explicitly states the remaining
novelty-review limitations.
