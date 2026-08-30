# Sharp rank-two theorem for two orthogonal latent directions

This note solves a structured rank-two projection class that contains the numerical rank-two maximizer found in the six-variable hafnian-gradient search.

Let `P` be a rank-two orthogonal projection on `R^6` that is, after a coordinate permutation and latent-space rotation, a direct sum of two rank-one blocks:

\[
P=aa^T\oplus bb^T,
\]

where the supports of `a` and `b` are disjoint, `||a||=||b||=1`, and the support sizes are `m` and `n=6-m`.

Define

\[
q_1(P)=\sum_{i<j}P_{ij}^2,
\qquad
q_2(P)=\sum_{i<j}\operatorname{haf}(P_{\widehat i,\widehat j})^2.
\]

Write `p_i=a_i^2` on the first support and `q_j=b_j^2` on the second.  Then `p` and `q` are probability vectors.

## Exact formulas

Because cross-block entries of `P` vanish,

\[
q_1=e_2(p)+e_2(q).
\tag{1}
\]

For a complementary four-set, a nonzero hafnian can occur only in one of three ways:

1. all four vertices lie in the first block, contributing the rank-one value;
2. all four lie in the second block;
3. two lie in each block, where exactly one perfect matching survives.

Therefore

\[
\boxed{
q_2=9e_4(p)+9e_4(q)+e_2(p)e_2(q).
}
\tag{2}
\]

## Maclaurin reduction

For a probability vector of support size `r`, Maclaurin gives

\[
e_4\le \binom r4\left(\frac{e_2}{\binom r2}\right)^2.
\tag{3}
\]

Also

\[
e_2\le\frac{r-1}{2r},
\tag{4}
\]

with equality exactly at equal weights.

Up to swapping the two blocks, only the support partitions `1+5`, `2+4`, and `3+3` occur.

### Partition `1+5`

The size-one block has `e_2=e_4=0`.  For the size-five block, (3) gives

\[
9e_4\le\frac9{20}e_2^2.
\]

Hence

\[
\frac{q_2}{q_1}\le\frac9{20}e_2\le\frac9{50}.
\]

Equality requires the five nonzero weights to be equal.

### Partition `2+4`

Put

\[
A=e_2(p),\qquad B=e_2(q).
\]

Then

\[
0\le A\le\frac14,
\qquad
0\le B\le\frac38,
\]

and the size-four Maclaurin bound gives `9e_4(q)<=B^2/4`.  Thus

\[
\frac{q_2}{q_1}
\le
\frac{AB+B^2/4}{A+B}.
\]

The right-hand side is increasing in both `A` and `B` on the indicated rectangle, so its maximum occurs at

\[
A=\frac14,
\qquad B=\frac38.
\]

Therefore

\[
\boxed{
\frac{q_2}{q_1}\le\frac{33}{160}=0.20625.
}
\tag{5}
\]

Equality requires equal weights in both blocks: two coordinates of squared weight `1/2` and four coordinates of squared weight `1/4`.

### Partition `3+3`

Both fourth elementary symmetric polynomials vanish, hence

\[
\frac{q_2}{q_1}=\frac{AB}{A+B}.
\]

With `A,B<=1/3`, the maximum is

\[
\frac16,
\]

attained at equal weights in both blocks.

## Theorem

Among all rank-two projections supported on two orthogonal latent directions,

\[
\boxed{
q_2(P)\le\frac{33}{160}q_1(P)<\frac14q_1(P).
}
\]

The sharp nondegenerate equality case is, up to coordinate permutations and sign switches,

\[
P=\frac12J_2\oplus\frac14J_4.
\]

Thus the strongest numerical rank-two configuration is not an accidental optimizer: it is the exact extremizer of the full two-direction class.

## Significance

The global rank-two problem is still open, but this theorem removes the most obvious candidate equality geometry and leaves a strict margin below the `1/4` constant needed for the projection-gradient program.  Any rank-two obstruction to the global `1/4` theorem must therefore use at least three genuinely different frame directions in the plane.
