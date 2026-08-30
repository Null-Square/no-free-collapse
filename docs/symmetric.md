# Exact symmetric order-1 construction

The universal Chebyshev bound is an upper bound. To study tightness we need explicit preparations whose optimal collapse can be solved exactly.

Assume `n=2L` is even and `kappa>1`. Use order-1 features `[1,x_1,...,x_n]` but set the constant coefficient to zero. Choose the variable Gram block

\[
Q_{\rm var}=\frac1n I+tJ,
\qquad t=\frac{\kappa-1}{n^2}.
\]

Then

\[
q(x)=x^TQ_{\rm var}x
=1+t\left(\sum_i x_i\right)^2.
\]

Since an even-dimensional cube contains inputs with zero sum and inputs with absolute sum `n`,

\[
q_{\min}=1,\qquad q_{\max}=\kappa.
\]

The preparation has Gram rank exactly `n`.

## Exact optimal full-parity collapse

Let

\[
\chi(x)=\prod_{i=1}^n x_i
\]

and

\[
S_n=\mathbb E\frac{\chi(x)}{q(x)}.
\]

Permutation symmetry reduces the signed readout operator to the trivial and standard representations of the symmetric group. Its one trivial eigenvalue has sign opposite to the `n-1` equal standard eigenvalues, and trace zero. The exact Helstrom optimum is

\[
C^{\rm mf}_{n}(\kappa)
=\frac{\kappa}{\kappa-1}|S_n|.
\]

Writing `n=2L` and

\[
\beta=\frac{n}{2\sqrt{\kappa-1}},
\]

the alternating binomial sum is

\[
S_n
=2^{-n}\sum_{j=0}^n
(-1)^j{n\choose j}
\frac{1}{1+\frac{\kappa-1}{n^2}(n-2j)^2}.
\]

Using the partial-fraction decomposition of `1/(u^2+beta^2)` together with

\[
\sum_{j=0}^n(-1)^j{n\choose j}\frac1{j+z}
=\frac{n!}{z(z+1)\cdots(z+n)},
\]

gives

\[
|S_n|
=
\frac{n!}{2^n\prod_{m=1}^{n/2}(m^2+\beta^2)}.
\]

Therefore

\[
\boxed{
C^{\rm mf}_{n}(\kappa)
=
\frac{\kappa}{\kappa-1}
\frac{n!}
{2^n
\prod_{m=1}^{n/2}
\left(m^2+\frac{n^2}{4(\kappa-1)}\right)}
}.
\]

The repository verifies this formula against the general Gram/eigenvalue optimizer on every tested small case.

This is currently a **constructive lower bound** on the unrestricted capacity, not a claim of global optimality.

## Asymptotic rate

For fixed `kappa>1`, Stirling's formula plus a Riemann-sum evaluation gives

\[
C^{\rm mf}_{n}(\kappa)
=\exp[-n I(\kappa)+O(\log n)],
\]

where

\[
\boxed{
I(\kappa)
=
\frac12\log\frac{\kappa}{\kappa-1}
+
\frac{\arctan\sqrt{\kappa-1}}{\sqrt{\kappa-1}}
}.
\]

So even though input-dependent normalization can create a formally nonzero interaction of arbitrary order, this explicit well-conditioned family has exponentially vanishing usable full-order signal.

## Open extremal problem

For order one and full parity, define the unrestricted conditioned capacity

\[
C_{n,1,n}(\kappa,d)
=\sup |\widehat p([n])|.
\]

We now have a rigorous sandwich for `d>=n`:

\[
C^{\rm mf}_{n}(\kappa)
\le C_{n,1,n}(\kappa,d)
\le
\frac{1}{2|T_{\lfloor(n-1)/2\rfloor}((\kappa+1)/(\kappa-1))|}.
\]

Random searches over nonsymmetric positive-semidefinite Gram matrices have not yet beaten the mean-field construction for the small cases tested, but this is evidence only. Proving or disproving global optimality is the next high-value theorem target.
