# First-order local optimality of the disjoint-pair hafnian construction

This note records a first-order theorem in the **unrestricted real PSD problem**. It is stronger than the matched pair-block optimality theorem: the perturbation direction is an arbitrary symmetric matrix, subject only to being the derivative of a differentiable feasible PSD path.

Let `n=2L` and define

\[
B_\star=\frac1{2n}\operatorname{diag}(J_2,\ldots,J_2),
\qquad
\Phi(B)=\operatorname{haf}(2\operatorname{offdiag}B).
\]

The matrix `B_star` is PSD, has rank `L`, satisfies

\[
\max_{x\in\{\pm1\}^n}x^TB_\star x=1,
\]

and attains

\[
\Phi(B_\star)=n^{-L}=n^{-n/2}.
\]

## Theorem

Let `B(epsilon)` be any differentiable one-sided feasible path with

\[
B(0)=B_\star,\qquad B(\epsilon)\succeq0,
\qquad \max_x x^TB(\epsilon)x\le1
\]

for all sufficiently small `epsilon>=0`. Write

\[
H=B'(0^+).
\]

Then

\[
\boxed{D\Phi(B_\star)[H]\le0.}
\]

Because `Phi(B_star)>0`, the same first-order conclusion holds for `|Phi|` in a neighborhood of `B_star`.

This is a first-order stationarity/local-optimality result. It does **not** by itself prove a strict local maximum, a second-order maximum, or unrestricted global optimality.

## Pair coordinates

For `j=1,...,L`, define

\[
f_j=e_{2j-1}+e_{2j},\qquad
g_j=e_{2j-1}-e_{2j}.
\]

The `f_j` span the range directions of `B_star`, while the `g_j` span its nullspace.

## Active Boolean constraints

Every Boolean vector whose two entries agree inside each matched pair is active at `B_star`. Such a vector can be written

\[
x(s)=\sum_{j=1}^L s_j f_j,\qquad s_j\in\{\pm1\},
\]

and satisfies `x(s)^T B_star x(s)=1`.

Feasibility of the path therefore gives, for every sign vector `s`,

\[
x(s)^T Hx(s)\le0.
\]

Average this inequality uniformly over all `s`. The cross terms vanish, leaving

\[
\boxed{\sum_{j=1}^L f_j^THf_j\le0.}
\]

## PSD tangent constraint on the nullspace

Since `g_j` lies in the nullspace of `B_star`, PSD feasibility gives

\[
0\le g_j^TB(\epsilon)g_j
=\epsilon\,g_j^THg_j+o(\epsilon).
\]

Hence

\[
g_j^THg_j\ge0
\]

for every `j`, and therefore

\[
\boxed{\sum_{j=1}^L g_j^THg_j\ge0.}
\]

## Matched-edge consequence

For each pair,

\[
4H_{2j-1,2j}=f_j^THf_j-g_j^THg_j.
\]

Summing and using the two inequalities above yields

\[
\boxed{\sum_{j=1}^L H_{2j-1,2j}\le0.}
\]

This is the key PSD-specific first-order restriction. The Boolean active-set condition alone does not supply the nonnegative nullspace term.

## Hafnian derivative

At `B_star`, the scaled off-diagonal matrix

\[
C_\star=2\operatorname{offdiag}(B_\star)
\]

has value `1/n` on each matched edge and zero on every unmatched edge. In the first derivative of the hafnian, only perturbations of the existing matched edges contribute. An unmatched edge requires at least one additional unmatched edge to complete a perfect matching, so it first appears at second order.

Therefore

\[
D\Phi(B_\star)[H]
=2n^{1-L}\sum_{j=1}^L H_{2j-1,2j}.
\]

Combining with the matched-edge inequality gives

\[
\boxed{D\Phi(B_\star)[H]\le0.}
\]

Equivalently, if

\[
A(H)=\sum_j f_j^THf_j,\qquad
N(H)=\sum_j g_j^THg_j,
\]

then

\[
D\Phi(B_\star)[H]
=\frac12 n^{1-L}\bigl(A(H)-N(H)\bigr),
\]

with `A(H)<=0` and `N(H)>=0` for every feasible tangent.

## CPU verification

The helper functions

- `disjoint_pair_tangent_invariants`,
- `disjoint_pair_hafnian_directional_derivative`

in `src/no_free_collapse/hafnian_bounds.py` encode the exact identities used above. The test file `tests/test_pair_local_optimality.py` checks:

1. the pair-sum/pair-difference identity for arbitrary symmetric directions;
2. the analytic derivative against direct finite differences of the hafnian;
3. explicit PSD- and cube-feasible paths whose derivatives are non-positive.

## Relevance to the six-variable problem

For `n=6`, this proves that the `1/216` disjoint-pair point cannot be improved to first order by **any** differentiable feasible PSD perturbation. A counterexample to global optimality, if one exists, must therefore arise through genuinely nonlinear/finite displacement or higher-order behavior rather than a simple improving tangent direction.

The next target is a second-order certificate or a global six-variable inequality, ideally

\[
54|\operatorname{haf}(C)|\le ts(t+s),
\]

with `t=Tr(B)` and `s=max_x sum_{i<j} C_ij x_i x_j`.
