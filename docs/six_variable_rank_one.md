# Six-variable candidate inequality: exact rank-one theorem

The main six-variable proof target is

\[
54\left|\operatorname{haf}(C)\right|\le t\,s\,(t+s),
\]

where

\[
C=2\operatorname{offdiag}(B),\qquad
t=\operatorname{tr}B,
\]

and

\[
s=\max_{x\in\{\pm1\}^6}\sum_{i<j}C_{ij}x_ix_j.
\]

Since

\[
t+s=\max_x x^TBx,
\]

this inequality would imply the desired global bound `|haf(C)|<=1/216` after normalizing the cube maximum to one.

This note proves the candidate inequality globally on the entire **rank-one PSD stratum**.

## Theorem

Let

\[
B=uu^T\succeq0,
\qquad u\in\mathbb R^6.
\]

Then

\[
\boxed{
54\left|\operatorname{haf}(2\operatorname{offdiag}B)\right|
\le t\,s\,(t+s).
}
\]

For nondegenerate vectors with all six coordinates nonzero, equality occurs exactly when the six absolute coordinate magnitudes are equal. Degenerate zero-slack cases give the trivial equality `0=0`.

## Proof

Write

\[
a_i=|u_i|,
\qquad
A=\sum_{i=1}^6 a_i.
\]

Because `B=uu^T`,

\[
x^TBx=(u^Tx)^2.
\]

The maximizing sign vector aligns every coordinate, so

\[
M:=\max_xx^TBx=A^2.
\]

Also

\[
t=\sum_i a_i^2,
\qquad
s=M-t=A^2-\sum_i a_i^2
=2e_2(a_1,\ldots,a_6).
\]

Here `e_k` denotes the elementary symmetric polynomial of degree `k`.

For the hafnian, every perfect matching contributes the same absolute monomial. There are `15` perfect matchings on six labeled vertices, and each matching contains three scaled edges `2u_i u_j`. Therefore

\[
\left|\operatorname{haf}(C)\right|
=15\cdot2^3\prod_{i=1}^6a_i
=120e_6(a).
\]

Maclaurin's inequalities give

\[
\left(\frac{e_6}{\binom66}\right)^{1/6}
\le
\left(\frac{e_2}{\binom62}\right)^{1/2}.
\]

Hence

\[
e_6
\le
\left(\frac{e_2}{15}\right)^3
=\left(\frac{s}{30}\right)^3.
\]

Thus

\[
54|\operatorname{haf}(C)|
\le
54\cdot120\left(\frac{s}{30}\right)^3
=\frac6{25}s^3.
\]

It remains to compare this with `t s M`.

By Cauchy-Schwarz,

\[
M=A^2\le6\sum_i a_i^2=6t,
\]

so

\[
t\ge\frac M6.
\]

Consequently

\[
s=M-t\le\frac{5M}{6}.
\]

If `s=0`, the claim is trivial. Otherwise,

\[
\frac6{25}s^2
\le
\frac6{25}\left(\frac{5M}{6}\right)^2
=\frac{M^2}{6}
\le tM.
\]

Multiplying by `s` gives

\[
\frac6{25}s^3\le tsM=ts(t+s),
\]

which completes the proof.

## Equality

For positive `a_i`, equality in Maclaurin and Cauchy-Schwarz requires

\[
a_1=\cdots=a_6.
\]

After normalizing `M=1`, this gives

\[
t=\frac16,
\qquad
s=\frac56,
\qquad
|\operatorname{haf}(C)|=\frac5{1944},
\]

which is exactly the six-variable all-parallel stationary family already seen numerically.

## Significance

The candidate inequality is now rigorously known on two important pieces of the geometry:

1. the full rank-one PSD boundary, by the theorem above;
2. the disjoint-pair point, where equality gives `1/216`, together with first- and second-order local certificates in the unrestricted PSD problem.

A global proof must interpolate between these geometrically very different equality structures rather than only controlling one of them.
