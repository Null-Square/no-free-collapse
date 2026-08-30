# Symmetric-class optimality theorem

This note closes the conditioned full-parity problem exactly inside a large, natural symmetry class. It does **not** claim unrestricted global optimality.

Let `n=2L>=4` be even. Use order-1 features `[1,x_1,...,x_n]`, and consider real Gram matrices `Q` invariant under all coordinate permutations and under global sign reversal `x -> -x`. The latter removes constant/linear cross terms. Such a Gram matrix has the representation

\[
Q=a\,|0\rangle\langle 0|+hP_{\rm std}+gP_{\rm triv},
\qquad a,h,g\ge0,
\]

where `P_triv` is the projector onto the all-ones variable direction and `P_std` is its orthogonal complement.

Writing

\[
s=\sum_i x_i,
\]

the squared norm is forced to be

\[
q(s)=a+hn+\frac{g-h}{n}s^2=A+t s^2,
\qquad A=a+hn.
\]

Thus every preparation in this symmetry class has only two possible nonconstant norm profiles: increasing or decreasing in `s^2`.

## Theorem

Among all such preparations with cube condition number at most `kappa>1`, followed by the **optimal possible two-outcome Born measurement**, the maximum full `n`-way parity coefficient is

\[
\boxed{
C^{\rm sym}_{n}(\kappa)=
\frac{\kappa}{\kappa-1}
\frac{n!}{2^n\prod_{m=1}^{n/2}
\left(m^2+\frac{n^2}{4(\kappa-1)}\right)}
}.
\]

It is attained by

\[
Q_{00}=0,\qquad
Q_{\rm var}=\frac1n I+\frac{\kappa-1}{n^2}J,
\]

so that

\[
q(x)=1+\frac{\kappa-1}{n^2}\left(\sum_i x_i\right)^2.
\]

This is exactly `mean_field_gram`.

## Proof: increasing branch

Suppose `t>=0`. Let the actual condition number be `kappa'<=kappa`; after rescaling,

\[
q(s)=A\left[1+(\kappa'-1)\frac{s^2}{n^2}\right].
\]

Write

\[
\lambda=a/A\in[0,1].
\]

PSD then forces

\[
h=\frac{A(1-\lambda)}n,
\qquad
g=\frac{A(\kappa'-\lambda)}n.
\]

Set

\[
R_0=\mathbb E\frac{\chi(x)}{1+(\kappa'-1)s^2/n^2},
\qquad \chi(x)=\prod_i x_i.
\]

Because `E chi=0`,

\[
\mathbb E\frac{\chi s^2}{1+(\kappa'-1)s^2/n^2}
=-\frac{n^2}{\kappa'-1}R_0.
\]

Permutation symmetry splits the signed Helstrom operator into one constant eigenvalue, one all-ones-variable eigenvalue, and an `(n-1)`-fold standard eigenvalue. They are

\[
\lambda R_0,
\qquad
-\frac{\kappa'-\lambda}{\kappa'-1}R_0,
\qquad
\frac{(1-\lambda)\kappa'}{(n-1)(\kappa'-1)}R_0.
\]

Their trace is zero. Therefore the exact optimal absolute Born coefficient is

\[
C_+(\lambda,\kappa')
=\frac{\kappa'-\lambda}{\kappa'-1}|R_0|.
\]

For fixed `kappa'`, this is strictly decreasing in `lambda`, so the optimal allocation is `lambda=0`: no constant feature norm.

The alternating-binomial identity gives

\[
|R_0|=
\frac{n!}{2^n\prod_{m=1}^{n/2}
\left(m^2+\frac{n^2}{4(\kappa'-1)}\right)}.
\]

Hence `C_+(0,kappa')` is the stated mean-field formula.

It is strictly increasing in `kappa'`. Writing `y=kappa'-1` and `L=n/2`, its logarithmic derivative reduces to

\[
\frac1y\left[
L^2\sum_{m=1}^{L}\frac1{L^2+m^2y}
-\frac1{1+y}
\right].
\]

The `m=L` term is exactly `1/(1+y)` and all earlier terms are positive, so the derivative is positive. Thus the full condition budget `kappa'=kappa` is optimal.

## Proof: decreasing branch

Suppose `t<0`. With actual condition number `kappa'`, rescale to

\[
q(s)=A\left[1-\frac{\kappa'-1}{\kappa'}\frac{s^2}{n^2}\right].
\]

Again let `lambda=a/A`. PSD now requires

\[
0\le\lambda\le1/\kappa'.
\]

The same representation split gives

\[
C_-(\lambda,\kappa')=(1-\lambda)C_-(0,\kappa'),
\]

so `lambda=0` is optimal. Its exact endpoint capacity is

\[
C_-(0,\kappa')=
\frac1{\kappa'-1}
\frac{n!}{2^n\prod_{m=1}^{L}
\left(\frac{\kappa'n^2}{4(\kappa'-1)}-m^2\right)}.
\]

The increasing endpoint strictly dominates it for `n>=4`. Put

\[
b=\frac{L^2}{\kappa'-1}.
\]

The ratio simplifies to

\[
\frac{C_+}{C_-}
=
\frac{\prod_{m=0}^{L-1}(b+L^2-m^2)}
{\prod_{j=1}^{L}(b+j^2)}.
\]

Pair numerator index `m` with denominator index `j=L-m`. Then

\[
L^2-m^2-(L-m)^2=2m(L-m)\ge0,
\]

with a strict inequality for some `m` whenever `L>=2`. Hence every paired numerator factor is at least the denominator factor and at least one is strictly larger. Therefore `C_+>C_-`.

Combining both branches proves the theorem.

## Why this matters

The unrestricted Gram objective is neither convex nor concave, so naive group averaging does not prove that an unrestricted optimum is symmetric. This theorem therefore gives a rigorous **symmetry-class optimum**, not the final global result.

The remaining problem is sharply defined:

> Can a symmetry-breaking real or complex Gram matrix beat `C_sym_n(kappa)`?

The repository's asymmetric searches have not found such a counterexample for the tested small cases, but only a proof or an explicit counterexample can close the unrestricted problem.
