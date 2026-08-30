# Hafnian upper bound from Boolean minimax approximation

This note isolates the first-order, low-conditioning object that appears when a quadratic normalization begins to generate a full parity interaction.

Let `m=2L` and let `B` be a real positive-semidefinite matrix. Define

\[
q(x)=x^T Bx,\qquad x\in\{\pm1\}^{m},
\]

and normalize so that

\[
0\le q(x)\le1
\]

on the Boolean cube. Put

\[
C=2\,\operatorname{offdiag}(B).
\]

## Coefficient identity

Because

\[
q(x)=\operatorname{tr}B+\sum_{i<j} C_{ij}x_i x_j,
\]

the full `m`-way Fourier coefficient of `q^L` can only arise by selecting exactly `L` quadratic monomials whose edges cover every vertex exactly once. Those selections are precisely perfect matchings. Every perfect matching can be ordered in `L!` ways, hence

\[
\boxed{
\widehat{q^L}([m])=L!\,\operatorname{haf}(C).
}
\]

This is the rigorous form of the observation that the first full-order term created by a quadratic normalizer is controlled by perfect-matchings combinatorics.

## Universal upper bound

Let `P` be any polynomial of degree at most `L-1`. Since `q` has Boolean degree at most two,

\[
\deg P(q)\le2L-2<m,
\]

so `P(q)` has zero full-parity Fourier coefficient. Therefore

\[
|\widehat{q^L}([m])|
=|\widehat{q^L-P(q)}([m])|
\le\|q^L-P(q)\|_\infty.
\]

The best degree-`L-1` uniform approximation to `t^L` on `[0,1]` is equivalent to the minimum sup norm of a monic degree-`L` polynomial. By the classical Chebyshev minimax theorem,

\[
\inf_{\deg P\le L-1}\|t^L-P(t)\|_{L_\infty[0,1]}
=2^{1-2L}.
\]

Combining with the coefficient identity yields

\[
\boxed{
|\operatorname{haf}(2\operatorname{offdiag}B)|
\le
\frac{2^{1-2L}}{L!}.
}
\]

The bound uses PSD only to ensure `q>=0`; the approximation step itself only needs `q(x) in [0,1]` and Boolean degree two.

## Exact four-variable theorem

For `L=2`,

\[
|\operatorname{haf}(2\operatorname{offdiag}B)|\le\frac1{16}.
\]

The disjoint-pair matrix with two rank-one blocks

\[
\frac18\begin{pmatrix}1&1\\1&1\end{pmatrix}
\]

has `max_x x^T Bx=1` and attains `1/16`. Hence the extremal problem is solved globally for four variables.

## Six-variable gap

For `m=6`, the theorem gives

\[
|\operatorname{haf}(2\operatorname{offdiag}B)|\le\frac1{192},
\]

while the equal disjoint-pair construction gives

\[
6^{-3}=\frac1{216}.
\]

So only a 12.5 percent gap remains. Numerical optimization over the larger class of globally-sign-symmetric degree-two functions `q:{+/-1}^6->[0,1]` can beat `1/216`, but the corresponding quadratic form is not PSD. Thus PSD geometry, rather than degree or range alone, is the remaining obstruction in the six-variable problem.

## Prior-art boundary

The ingredients are classical or established separately:

- hafnians are sums over perfect matchings;
- Chebyshev gives the minimax monic polynomial on an interval;
- Fourier coefficients are bounded by uniform approximation error.

The research object here is the constrained composition of these facts for the Boolean-cube normalization problem arising in collapse-based reasoning. The next theorem target is to close the PSD gap for six variables and then determine whether the disjoint-pair value `m^(-m/2)` is globally optimal for all even `m`.
