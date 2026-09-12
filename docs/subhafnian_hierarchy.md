# Spectral hierarchy for averaged principal subhafnian energies

**Status: analytic theorem, September 12, 2026.** This note derives an all-orders
consequence of the complementary-subset dilation. It is an upper bound for
averaged **principal** subhafnian energies, not a lower bound for repeated-index
PSD hafnians. All main recurrence statements below have integer r>=2.

## 1. Definitions

Let `A=A^T` be a real symmetric `N x N` matrix and put

\[
d(A)=\lambda_{\max}(A)-\lambda_{\min}(A).
\]

For integers `1<=r<=N/2`, define

\[
H_{2r}(A)=\sum_{\substack{S\subseteq[N]\\|S|=2r}}
|\operatorname{haf}(A[S])|^2,
\qquad
\overline H_{2r}(A)=\frac{H_{2r}(A)}{\binom N{2r}}.
\]

In Roos's normalization,

\[
G_r(A)=\frac{\overline H_{2r}(A)}{((2r-1)!!)^2}.
\]

Thus `G_1=H_2/binom(N,2)` is the average squared off-diagonal entry.

## 2. Critical-dimension theorem

Let integer `r>=2` and let `B` be real symmetric of order `4r-2`. Then

\[
\boxed{H_{2r}(B)\le\frac{d(B)^2}{4}H_{2r-2}(B).}\tag{1}
\]

The constant `1/4` is sharp for every r.

### Proof from the weighted Kneser dilation

Set `m=2r-1`. Index a vector v by the `(2r-2)`-subsets E and put
`v_E=haf(B[E])`. For fixed E, its complement has size 2r, and the hafnian
Laplace expansion gives

\[
(T_m(B)v)_E
=\sum_{\{i,j\}\subseteq E^c}B_{ij}\operatorname{haf}
(B[E^c\setminus\{i,j\}])
=r\,\operatorname{haf}(B[E^c]).\tag{2}
\]

Every perfect matching is counted once for each of its r edges. Therefore

\[
\|v\|^2=H_{2r-2}(B),\qquad
\|T_m(B)v\|^2=r^2H_{2r}(B).
\]

The sharp dilation theorem gives `||T_m(B)||<=r d(B)/2`. Substitution and
cancellation of r squared proves (1).

### Sharp families

At `N=4r-2`, all three PSD contractions `ss^T/N`, `I-ss^T/N`, and `(I+M)/2`
attain equality, where s is a sign vector and M a signed perfect matching.
For the rank-one family,

\[
H_{2q}=\binom N{2q}\frac{((2q-1)!!)^2}{N^{2q}}.
\]

At critical order `binom(N,2r)=binom(N,2r-2)` and `(2r-1)^2/N^2=1/4`.
The complement has the same squared energies. For the matching family,

\[
H_{2q}=\binom{N/2}{q}4^{-q}.
\]

Since `N/2=2r-1`, consecutive binomial coefficients at `q=r-1,r` agree.
These examples prove sharpness. No complete equality classification for r>=3
is claimed.

## 3. Arbitrary ambient dimension

For `N>=4r-2`, apply (1) to each principal `(4r-2)`-subset U. Cauchy
interlacing gives `d(A[U])<=d(A)`. Each 2r-subset is counted
`binom(N-2r,2r-2)` times, whereas each `(2r-2)`-subset is counted
`binom(N-2r+2,2r)` times. Their ratio is

\[
\frac{\binom{N-2r+2}{2r}}{\binom{N-2r}{2r-2}}
=\frac{\binom N{2r}}{\binom N{2r-2}}.
\]

Thus

\[
\boxed{\overline H_{2r}(A)\le\frac{d(A)^2}{4}\overline H_{2r-2}(A).}\tag{3}
\]

Iteration gives

\[
\overline H_{2r}(A)\le
\left(\frac{d(A)^2}{4}\right)^{r-1}\overline H_2(A).\tag{4}
\]

The dimension-free one-step factor is sharp across the admissible `(N,r)`.
This is not a claim of the best fixed-N constant for `N>4r-2`, nor of a sharp
iterated product. For d=0 all positive-order energies vanish.

## 4. Roos comparison, including odd ambient order

Equation (3) gives

\[
G_r\le\frac{d(A)^2}{4(2r-1)^2}G_{r-1}.\tag{5}
\]

Roos's Theorem 2.3 states logarithmic subadditivity in even ambient order.
The following application of his parity-free Corollary 3.1 justifies
`G_{a+b}<=G_a G_b` for arbitrary `N>=2(a+b)`.

Put `F_j(S)=haf(A[S])/(2j-1)!!` and `r=a+b`. Matching enumeration gives

\[
\sum_{T\subset S,\ |T|=2a}\operatorname{haf}(A[T])
\operatorname{haf}(A[S\setminus T])=\binom ra\operatorname{haf}(A[S]).
\]

Since `C(r,a)(2r-1)!!=C(2r,2a)(2a-1)!!(2b-1)!!`, this is the normalized
subset convolution of F_a and F_b. Roos Corollary 3.1 with ambient [N],
subset size 2r and block sizes (2a,2b) bounds its average squared modulus
by `G_a G_b`, with no ambient parity restriction. This is an application of
Roos, not a claim to invent a replacement inequality.

Combining the Roos step and (5),

\[
\boxed{G_r\le G_{r-1}
\min\{G_1,d(A)^2/[4(2r-1)^2]\}.}\tag{6}
\]

Consequently

\[
\boxed{G_r\le G_1\prod_{s=2}^r
\min\{G_1,d(A)^2/[4(2s-1)^2]\}.}\tag{7}
\]

For d>0 the spectral factors strictly decrease, while G_1 is fixed; there is
at most one switch. If k is the last Roos-selected order, assigning ties to
Roos, the product is

\[
G_1^k\left(\frac{d(A)^2}{4}\right)^{r-k}
\left(\frac{(2k-1)!!}{(2r-1)!!}\right)^2.\tag{8}
\]

If G_1=0, every positive-order energy is zero. Otherwise denote (7)'s bound
by U_s and its factors by `m_j=min(G_1,d^2/[4(2j-1)^2])`. For a>=b,

\[
\frac{U_{a+b}}{U_aU_b}
=\frac{\prod_{j=a+1}^{a+b}m_j}{G_1\prod_{j=2}^{b}m_j}\le1.
\]

Here `m_{a+1}<=G_1` and `m_{a+j}<=m_j`. Therefore the hybrid is log-subadditive
within its admissible order range. Larger Roos partitions applied to these
same lower-order bounds do not improve it; no broader capacity optimality
claim is intended.

## 5. Sharp fixed-dimension density crossover

For d(A)>0 let `theta=4G_1/d(A)^2`. The spectral step is strictly smaller
than the Roos step exactly when `theta>1/(2r-1)^2`. The sharp fixed-N ceiling is

\[
\boxed{\theta\le\frac{4\lfloor N^2/4\rfloor}{N^2(N-1)}
=\begin{cases}1/(N-1),&N\text{ even},\\(N+1)/N^2,&N\text{ odd}.\end{cases}}\tag{9}
\]

Normalize `B=(A-lambda_min I)/d`, whose eigenvalues are in [0,1]. Then

\[
2H_2(B)=\operatorname{tr}(B^2)-\sum_i B_{ii}^2
\le\sum_i\lambda_i^2-\frac{(\sum_i\lambda_i)^2}{N}.
\]

The last expression is strictly convex in each individual eigenvalue. Its
maximum occurs only at cube vertices, with value `k(N-k)/N` when k entries
are one. The maximizing ranks are floor(N/2), ceil(N/2). Equality additionally
requires constant diagonal. This gives exactly the stated shifted/scaled
constant-diagonal projection equality locus. Matching projections establish
existence in even dimensions; real Fourier sine/cosine pairs, plus the constant
column when the rank is odd, establish it in odd dimensions.

At critical `N=4r-2`, balanced rank one has `theta=1/(2r-1)^2`, exactly the
crossover. Matching projections have `theta=1/(4r-3)` and attain the spectral
step. The ratio of the particular Roos coefficient to the spectral coefficient
is `(2r-1)^2/(4r-3)`, which grows linearly with r. This is not a comparison
against every possible prior inequality.

## 6. Restricted collision-free Gaussian boson sampling

For a pure, zero-displacement state the collision-free photon-number-resolving
probability is `p_B(S)=C(B)|haf(B[S])|^2`, with common positive C(B). For real
symmetric B, integer r>=2 and N>=4r-2, it cancels from (3), yielding

\[
\frac{\sum_{|S|=2r}p_B(S)}{\binom N{2r}}
\le\frac{d(B)^2}{4}
\frac{\sum_{|S|=2r-2}p_B(S)}{\binom N{2r-2}}.\tag{10}
\]

Physical `||B||<1` gives a factor at most `||B||^2<1`. In the real-interferometer
PSD subclass `B=O diag(tanh s_j) O^T` with nonnegative s_j, it is at most
`||B||^2/4<1/4`. This bounds averaged unconditional probability per pattern,
not necessarily decreasing total sector mass. It makes no statement about
collisionful repetitions, threshold-detector Torontonians, or experimental
computational advantage. Hamilton et al. and Kruse et al. supply the probability
formula; the manuscript includes their bibliographic details.

## 7. Prior-art and verification scope

Primary comparator: B. Roos, *New inequalities for permanents and hafnians and
some generalizations*, Linear and Multilinear Algebra 73 (2025), 1634-1667,
DOI 10.1080/03081087.2024.2436061; arXiv:1906.06176v2, Theorem 2.3 and
Corollary 3.1. The spectral step uses eigenvalue information absent from that
specific entry-energy recurrence. Kummer's ordinary Kneser construction and
Feinsilver's Boolean-lattice representation are relevant background. The full
Lindell--Staples text remains an identified access/priority-review gap; see the
[literature log](../submission/literature_review.md).

The tests verify the integer hafnian recursion, global multiplicities, sharp
families through r=4, odd-order normalized convolutions, exact rational spectral
intervals, hybrid formulas, fixed-N density witnesses and boundary rejection.
The floating enumerator is capped at N=20 and its read-only cache is tied to
matrix entries; exact rational enumeration is capped at N=14. These caps limit
computation, not the analytic theorem. See [the audit](SUBMISSION_AUDIT.md).
