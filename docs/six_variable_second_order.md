# Six-variable second-order certificate at the disjoint-pair point

This note strengthens the first-order theorem specifically for the six-variable hafnian problem.

Let

\[
B_\star=\frac1{12}\operatorname{diag}(J_2,J_2,J_2),
\qquad
\Phi(B)=\operatorname{haf}(2\operatorname{offdiag}B).
\]

Then

\[
\Phi(B_\star)=\frac1{216}.
\]

The feasible set is

\[
\mathcal F=\left\{B\succeq0:\ \max_{x\in\{\pm1\}^6}x^TBx\le1\right\}.
\]

The previous first-order theorem proves

\[
D\Phi(B_\star)[H]\le0
\]

for every feasible tangent direction `H`. Here we characterize the directions where equality occurs and prove that every nonzero critical tangent has strictly negative second variation along any twice-differentiable feasible path.

## Orthonormal pair coordinates

For `j=1,2,3`, define

\[
F_j=\frac{e_{2j-1}+e_{2j}}{\sqrt2},\qquad
G_j=\frac{e_{2j-1}-e_{2j}}{\sqrt2}.
\]

In the ordered basis `(F_1,F_2,F_3,G_1,G_2,G_3)`,

\[
B_\star=
\begin{pmatrix}
\frac16 I_3&0\\
0&0
\end{pmatrix}.
\]

Write a symmetric tangent as

\[
H=
\begin{pmatrix}
H_{FF}&H_{FG}\\
H_{FG}^T&H_{GG}
\end{pmatrix}.
\]

## Critical first-order directions

Assume `H` is feasible to first order and

\[
D\Phi(B_\star)[H]=0.
\]

From PSD tangent feasibility on the nullspace,

\[
H_{GG}\succeq0.
\]

From the active Boolean vertices, for every `s in {+/-1}^3`,

\[
s^T H_{FF}s\le0.
\]

The first-order proof expresses the hafnian derivative as the difference between the average active curvature and the nullspace trace. Equality of the derivative therefore forces both nonnegative losses to vanish. Consequently

\[
H_{GG}=0,
\]

and every active first derivative is exactly zero:

\[
s^T H_{FF}s=0
\qquad\forall s\in\{\pm1\}^3.
\]

A quadratic form that vanishes on every sign vector has zero off-diagonal entries and zero trace. Hence every critical tangent has the form

\[
\boxed{
H=
\begin{pmatrix}
D&X\\
X^T&0
\end{pmatrix},
\qquad
D=\operatorname{diag}(d_1,d_2,d_3),
\qquad d_1+d_2+d_3=0.
}
\]

The matrix `X` is arbitrary at first order.

## Exact pure-H hafnian curvature

A direct expansion of the 15 perfect matchings in the six-variable hafnian gives

\[
\Phi(B_\star+\varepsilon H)
=\frac1{216}
+\varepsilon^2 Q(H)
+O(\varepsilon^3),
\]

where

\[
\boxed{
Q(H)
=-\frac{d_1^2+d_2^2+d_3^2}{12}
-\frac13\sum_{j\ne k}X_{jk}^2.
}
\]

The diagonal entries `X_jj` do not appear in this direct quadratic hafnian term. Geometrically, those directions mix the range and nullspace inside the same original pair without changing an off-diagonal hafnian edge to second order by themselves.

## PSD curvature supplies the missing penalty

Now take an arbitrary twice-differentiable feasible path

\[
B(\varepsilon)
=B_\star+\varepsilon H+\frac{\varepsilon^2}{2}K+o(\varepsilon^2),
\]

with a nonzero critical `H` as above.

### PSD second-order constraint

Using the Schur complement in `(F,G)` coordinates,

\[
\frac{\varepsilon^2}{2}K_{GG}
-\varepsilon^2 X^T(\tfrac16 I_3)^{-1}X
+o(\varepsilon^2)\succeq0.
\]

Therefore

\[
\boxed{K_{GG}\succeq12X^TX.}
\]

In particular,

\[
\operatorname{tr}K_{GG}\ge12\|X\|_F^2.
\]

### Active-cube second-order constraint

Every pair-aligned Boolean vertex was active at `B_star`, and its first derivative vanishes for a critical `H`. The second derivative must therefore be non-positive. Averaging over the eight active sign vectors gives

\[
\boxed{\operatorname{tr}K_{FF}\le0.}
\]

### Effect on the matched hafnian edges

For each original pair,

\[
K_{2j-1,2j}
=\frac12\left((K_{FF})_{jj}-(K_{GG})_{jj}\right).
\]

Hence

\[
\sum_{j=1}^3K_{2j-1,2j}
=\frac12\left(\operatorname{tr}K_{FF}-\operatorname{tr}K_{GG}\right)
\le-6\|X\|_F^2.
\]

At `n=6`, the first derivative formula is

\[
D\Phi(B_\star)[K]
=\frac1{18}\sum_j K_{2j-1,2j}.
\]

Because the path contains `(epsilon^2/2)K`, the `K` contribution to the coefficient of `epsilon^2` is at most

\[
-\frac16\|X\|_F^2.
\]

## Strict second-order theorem

Combining the direct hafnian curvature with the forced PSD curvature gives

\[
\boxed{
[\varepsilon^2]\,\Phi(B(\varepsilon))
\le
-\frac{\|d\|_2^2}{12}
-\frac{\|\operatorname{offdiag}X\|_F^2}{3}
-\frac{\|X\|_F^2}{6}.
}
\]

Therefore, for every nonzero critical tangent `H`,

\[
\boxed{[\varepsilon^2]\,\Phi(B(\varepsilon))<0.}
\]

So a smooth feasible path cannot leave the pair point with zero first-order loss and then improve at second order.

This is stronger than first-order stationarity, but it still does **not** prove unrestricted global optimality.

## Exact zero-slope segment corollary

There is a useful stronger statement for straight feasible segments.

Let `B in F` and set `H=B-B_star`. Since the feasible set is convex, the entire segment `B_star+epsilon H` is feasible for `0<=epsilon<=1`. If

\[
D\Phi(B_\star)[H]=0,
\]

then the linear path has no second-order correction `K`. The PSD curvature condition above therefore forces `X=0`. Thus

\[
H_{FF}=\operatorname{diag}(d_1,d_2,d_3),
\qquad
\sum_jd_j=0,
\]

and no other blocks are present.

In original coordinates, the only nonzero hafnian edges are the three matched edges with weights

\[
\frac16+d_1,\qquad \frac16+d_2,\qquad \frac16+d_3.
\]

PSD feasibility makes these weights nonnegative, and their sum is `1/2`. Therefore AM-GM gives

\[
\Phi(B)
=\prod_{j=1}^3\left(\frac16+d_j\right)
\le\left(\frac16\right)^3
=\frac1{216},
\]

with equality only at `B=B_star`.

Hence any global counterexample connected to `B_star` by a straight segment must begin with a **strictly negative** first derivative and later rebound; it cannot lie on a zero-slope face.

## CPU verification

The implementation in `src/no_free_collapse/hafnian_bounds.py` provides:

- `disjoint_pair_basis_blocks`;
- `six_variable_critical_hafnian_quadratic_coefficient`;
- `six_variable_critical_path_second_order_bound`.

The tests in `tests/test_six_variable_second_order.py` verify:

1. the critical block decomposition;
2. the exact quadratic coefficient by evaluating the cubic hafnian polynomial at `epsilon=+1` and `epsilon=-1`;
3. strict negativity for random nonzero critical tangents;
4. an exact Schur-complement PSD path that attains the second-order curvature bound and remains cube-feasible for small positive and negative perturbations.

## Research implication

The six-variable search is now constrained in three layers:

1. the pair point is first-order optimal in the unrestricted PSD problem;
2. every nonzero zero-slope tangent has strictly negative second variation;
3. every zero-slope straight feasible segment is globally non-improving by AM-GM.

The remaining global problem is to rule out a finite-displacement configuration whose segment from `B_star` initially decreases and later exceeds `1/216`, or to find such a configuration.
