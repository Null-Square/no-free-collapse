# Mathematical core

## Boolean interaction basis

For `x in {-1,+1}^n`, write `x_A = product_{i in A} x_i`. Every real-valued function on the Boolean hypercube has a unique multilinear Walsh expansion

\[
f(x)=\sum_{A\subseteq[n]} \widehat f(A)x_A.
\]

Its interaction degree is the largest `|A|` with nonzero `f_hat(A)`. Changing from the `{-1,+1}` Walsh basis to the `{0,1}` multilinear/Mobius-Harsanyi basis preserves polynomial degree, so the order cutoff is basis-independent in the sense needed here.

## Fixed-norm Born interaction bound

Let an amplitude preparation have order at most `r`:

\[
|\psi(x)\rangle=\sum_{|A|\le r}x_A|u_A\rangle.
\]

Assume it is already physically normalized for every input, or more generally has an input-independent norm. For a fixed Hermitian effect `M`,

\[
p(x)=\langle\psi(x)|M|\psi(x)\rangle
=\sum_{A,B}x_Ax_B\langle u_A|M|u_B\rangle.
\]

On the Boolean cube, `x_A x_B = x_{A triangle B}` and

\[
|A\triangle B|\le |A|+|B|\le 2r.
\]

Therefore

\[
\deg p\le 2r.
\]

This is the exact ceiling tested in `test_degree_bound.py`.

## Normalization is not free

For an arbitrary polynomial vector `psi_tilde(x)`, a physical probability is

\[
p(x)=\frac{\langle\tilde\psi|M|\tilde\psi\rangle}{\langle\tilde\psi|\tilde\psi\rangle}.
\]

If the denominator depends on `x`, this is a rational function on the hypercube and may contain interactions above order `2r` after multilinear reduction. Thus a theorem stated only in terms of the *unnormalized numerator* would not characterize an implementation that performs input-dependent normalization.

An exact counterexample is

\[
\tilde\psi(x)=(2+x_1+x_2+x_3,\;1),\qquad M=|0\rangle\langle0|.
\]

The amplitudes are order 1, but after normalization

\[
p(x)=\frac{(2+x_1+x_2+x_3)^2}{(2+x_1+x_2+x_3)^2+1}
\]

has cubic Walsh coefficient

\[
\widehat p(\{1,2,3\})=-\frac{6}{65}.
\]

The project therefore treats input-dependent normalization as part of state-preparation computation, not as a free readout operation.

## Tight parity witness

Let `A,B` partition `[k]` with both sizes at most `ceil(k/2)`, and define `u=x_A`, `v=x_B`. The 2D state

\[
|\psi(x)\rangle=\frac{u+v}{2}|0\rangle+\frac{u-v}{2}|1\rangle
\]

has norm 1 for every input and preparation order `ceil(k/2)`. Measuring `|0><0|` gives

\[
P(0)=\frac{(u+v)^2}{4}=\frac{1+uv}{2}=\frac{1+\prod_i x_i}{2}.
\]

Hence the factor-of-two upper bound is attained by the parity family.

## Status

The degree calculation itself is standard polynomial-method mathematics. The research contribution under investigation is a reasoning-specific resource framework: what pre-collapse interaction resources are required by reasoning objectives, how normalization/nonlinearity should be charged, and which stronger approximation or sparsity lower bounds survive comparison with unrestricted classical continuous latent computation.

## Normalization explosion theorem

The three-variable loophole extends to arbitrary order. For any `n>=1`, let

\[
\widetilde\psi_n(x)=\left(c+\sum_{i=1}^n x_i,\;1\right),\qquad c=2n+1,
\]

and measure the first coordinate after normalization. The amplitudes before normalization are affine (interaction order 1), while

\[
p_n(x)=\frac{(c+\sum_i x_i)^2}{(c+\sum_i x_i)^2+1}.
\]

The full Walsh coefficient is

\[
\widehat p_n([n])
=-2^{-n}\sum_{j=0}^n(-1)^j\binom nj
\frac{1}{(c+n-2j)^2+1}.
\]

Set `a=c-n=n+1` and reindex. Using

\[
\frac1{y^2+1}=\frac1{2i}\left(\frac1{y-i}-\frac1{y+i}\right)
\]

and the finite-difference identity

\[
\sum_{k=0}^n(-1)^k\binom nk\frac1{z+k}
=\frac{n!}{z(z+1)\cdots(z+n)},
\]

one obtains a nonzero closed form proportional to the imaginary part of

\[
\frac1{\left((a-i)/2\right)_{n+1}}.
\]

For `a=n+1`, its argument is

\[
\theta=\sum_{k=0}^n\arctan\frac1{n+1+2k},
\]

with `0<theta<1<pi`. Hence the imaginary part is strictly positive and

\[
\widehat p_n([n])\neq0
\]

for every `n`. Its sign alternates as `(-1)^{n+1}`.

Thus input-dependent normalization can amplify an order-1 unnormalized preparation to *full interaction order n*, with no fixed ceiling as `n` grows. Any no-free-collapse resource theory must therefore include normalization/nonlinear renormalization in the charged preparation computation.
