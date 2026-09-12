# Spectral hierarchy for averaged principal subhafnian energies

**Status: analytic theorem, September 12, 2026.** This note derives an all-orders
consequence of the complementary-subset dilation. It is an upper bound for
averaged **principal** subhafnian energies. It is different from lower bounds
for repeated-index PSD hafnians arising from Gaussian product inequalities, and
from Roos's log-subadditivity theorem. The three results should not be conflated.

## 1. Definitions

Let `A=A^T` be a real symmetric `N x N` matrix and put

\[
d(A)=\lambda_{\max}(A)-\lambda_{\min}(A).
\]

For `r>=1`, define the order-`2r` principal subhafnian energy

\[
H_{2r}(A)=\sum_{\substack{S\subseteq[N]\\|S|=2r}}
|\operatorname{haf}(A[S])|^2,
\qquad
\overline H_{2r}(A)=\frac{H_{2r}(A)}{\binom N{2r}}.
\]

Following the normalization in Roos's Theorem 2.3, write

\[
G_r(A)=\frac{\overline H_{2r}(A)}{((2r-1)!!)^2}.
\]

Thus `G_1=H_2/binom(N,2)` is the average squared off-diagonal entry.

## 2. Critical-dimension theorem

Let `r>=2` and let `B` have order `4r-2`. Then

\[
\boxed{H_{2r}(B)\le\frac{d(B)^2}{4}H_{2r-2}(B).}\tag{1}
\]

The constant `1/4` is sharp for every `r`.

### Proof from the weighted Kneser dilation

Set `m=2r-1`. Index a vector `v` by the `(2r-2)`-subsets `E` and put

\[
v_E=\operatorname{haf}(B[E]).
\]

The operator `T_m(B)` from `matching_dilation.md` is indexed by these same
subsets. For fixed `E`, the complement has size `2r`, and the hafnian Laplace
expansion gives

\[
(T_m(B)v)_E
=\sum_{\{i,j\}\subseteq E^c}B_{ij}\operatorname{haf}
  (B[E^c\setminus\{i,j\}])
=r\,\operatorname{haf}(B[E^c]).\tag{2}
\]

Every perfect matching of `E^c` is counted once for each of its `r` edges.
Hence

\[
\|v\|^2=H_{2r-2}(B),\qquad
\|T_m(B)v\|^2=r^2H_{2r}(B).
\]

The sharp dilation theorem gives

\[
\|T_m(B)\|\le\frac{m+1}{4}d(B)=\frac r2d(B).
\]

Substitution and cancellation of `r^2` proves (1).

### Sharp families

At order `N=4r-2`, all three of the following PSD contractions attain equality:

1. the balanced signed rank-one projection `ss^T/N`;
2. its rank-`N-1` complement `I-ss^T/N`;
3. `(I+M)/2`, where `M` is a signed perfect-matching matrix.

For the rank-one family,

\[
H_{2s}=\binom N{2s}\frac{((2s-1)!!)^2}{N^{2s}}.
\]

At the critical order `binom(N,2r)=binom(N,2r-2)` and
`(2r-1)^2/N^2=1/4`. The complement has the same even-order energies.
For the matching family,

\[
H_{2s}=\binom{N/2}{s}4^{-s};
\]

when `N/2=2r-1`, the consecutive binomial coefficients at `s=r-1,r`
are equal, again giving the ratio `1/4`.

These examples prove sharpness. No complete equality classification for (1)
is claimed here when `r>=3`.

## 3. Averaged hierarchy in arbitrary ambient dimension

Let `N>=4r-2`. Applying (1) to every principal `(4r-2)`-subset `U` and
using Cauchy interlacing, `d(A[U])<=d(A)`, yields

\[
\boxed{
\overline H_{2r}(A)
\le\frac{d(A)^2}{4}\,\overline H_{2r-2}(A).
}\tag{3}
\]

Indeed, each order-`2r` subset is contained in
`binom(N-2r,2r-2)` critical subsets, whereas each order-`2r-2` subset is
contained in `binom(N-2r+2,2r)` of them. Their ratio is exactly

\[
\frac{\binom{N-2r+2}{2r}}{\binom{N-2r}{2r-2}}
=\frac{\binom N{2r}}{\binom N{2r-2}},
\]

which is why binomial normalization removes all ambient-dimension factors.

Iterating (3) gives

\[
\boxed{
\overline H_{2r}(A)
\le\left(\frac{d(A)^2}{4}\right)^{r-1}\overline H_2(A),
\qquad N\ge4r-2.
}\tag{4}
\]

The dimension-free one-step constant `1/4` is optimal over the stated family
of admissible `(N,r)`, because equality occurs at every critical `N=4r-2`.
The theorem does **not** assert that `1/4` is the best fixed-`N` constant when
`N>4r-2`.

## 4. A spectral recurrence for Roos's normalized energy

Since `((2r-1)!!)^2=(2r-1)^2((2r-3)!!)^2`, (3) becomes

\[
\boxed{
G_r(A)\le\frac{d(A)^2}{4(2r-1)^2}G_{r-1}(A),
\qquad N\ge4r-2.
}\tag{5}
\]

Roos's Theorem 2.3 gives the distinct log-subadditive recurrence

\[
G_r(A)\le G_{r-1}(A)G_1(A).\tag{6}
\]

Combining (5) and (6) gives the certified hybrid bound

\[
\boxed{
G_r(A)\le G_{r-1}(A)
\min\left\{G_1(A),\frac{d(A)^2}{4(2r-1)^2}\right\}.
}\tag{7}
\]

Consequently

\[
\boxed{
G_r(A)\le G_1(A)\prod_{s=2}^r
\min\left\{G_1(A),\frac{d(A)^2}{4(2s-1)^2}\right\}.
}\tag{8}
\]

Because the spectral factor decreases strictly with `s`, the minimum in (8)
has **at most one switch**. Let `k` be the largest order at which the Roos
factor is still selected (ties assigned to Roos). Then

\[
\boxed{
G_r(A)\le G_1(A)^k
\left(\frac{d(A)^2}{4}\right)^{r-k}
\left(\frac{(2k-1)!!}{(2r-1)!!}\right)^2.
}\tag{9}
\]

The resulting upper-bound sequence is itself log-subadditive: if `U_s` denotes
the right side of (8), then `U_{a+b}<=U_aU_b`. Therefore applying a larger
Roos partition to already bounded lower orders cannot improve (8) using only
these same ingredients.

## 5. Spectral-density crossover

Define the scale-invariant edge-density parameter

\[
\theta(A)=\frac{4G_1(A)}{d(A)^2}
\]

when `d(A)>0`. The spectral recurrence (5) is strictly stronger than the
Roos step (6) exactly when

\[
\boxed{\theta(A)>\frac1{(2r-1)^2}.}\tag{10}
\]

There is also the sharp elementary universal bound

\[
\boxed{\theta(A)\le\frac1{N-1}.}\tag{11}
\]

To prove it, shift and scale without changing off-diagonal entries so that
`0<=B<=I`. If `t=tr(B)`, then

\[
2H_2(B)=\operatorname{tr}(B^2)-\sum_iB_{ii}^2
\le t-\frac{t^2}{N}\le\frac N4.
\]

Thus `H_2<=N/8`, so `G_1<=1/[4(N-1)]`. Equality in every step occurs
**exactly** when `N` is even and `B` is a rank-`N/2` orthogonal projection
with constant diagonal `1/2`.

This identifies the natural window in which the new step can improve Roos:

\[
4r-2\le N\quad\text{and}\quad N-1<(2r-1)^2.
\]

At the critical dimension `N=4r-2`, the balanced rank-one projection has
`theta=1/(2r-1)^2`, so it lies exactly on the crossover and attains both
recurrences. The matching projection has the maximal `theta=1/(4r-3)`, and
the spectral coefficient improves the Roos coefficient by

\[
\boxed{\frac{(2r-1)^2}{4r-3}.}\tag{12}
\]

This factor is `9/5` at `r=2`, `25/9` at `r=3`, and grows asymptotically like
`r`. Thus the new estimate is not merely a reparameterization of Roos's bound:
each theorem is sharp on a different extremal geometry.

## 6. Collision-free Gaussian boson sampling corollary

For a pure zero-displacement Gaussian boson sampler, a collision-free outcome
`S` has probability

\[
p_B(S)=C(B)\,|\operatorname{haf}(B[S])|^2,
\]

where the normalization `C(B)` is independent of the collision-free pattern.
Therefore our hierarchy immediately becomes a statement about **unconditional
average per-pattern probability** in low collision-free photon sectors. If `B`
is a real symmetric sampling matrix on `N` modes and `N>=4r-2`, then

\[
\boxed{
\frac{1}{\binom N{2r}}\sum_{|S|=2r}p_B(S)
\le\frac{d(B)^2}{4}
\frac{1}{\binom N{2r-2}}\sum_{|S|=2r-2}p_B(S).
}\tag{13}
\]

Equivalently, if `P_{2r}^{cf}` is the total collision-free probability mass of
the `2r`-photon sector,

\[
P_{2r}^{cf}\le\frac{d(B)^2}{4}
\frac{\binom N{2r}}{\binom N{2r-2}}P_{2r-2}^{cf}.\tag{14}
\]

For any real symmetric physical sampling matrix, `||B||<1` and
`d(B)<=2||B||`, so (13) contracts average per-pattern weight by at most
`||B||^2<1`. In the real-interferometer subclass
`B=O diag(tanh s_j) O^T`, one has `B>=0` and thus `d(B)<=||B||`; the contraction
factor improves to `||B||^2/4<1/4`.

This is deliberately an average-sector statement. The total sector mass need
not decrease, because the number of collision-free patterns grows with photon
number. No claim about threshold-detector Torontonians or collisionful output
patterns is made.

## 7. Prior-art boundary

Roos, *New inequalities for permanents and hafnians and some generalizations*,
Linear and Multilinear Algebra 73 (2025), Theorem 2.3, defines the same `G_r`
normalization and proves log-subadditivity. In particular,
`G_r<=G_{r-1}G_1`, with equality when all off-diagonal entries are identical.
The spectral recurrence (5), its edge-density crossover, and the hybrid bound
use eigenvalue information absent from that statement.

Ouimet and Greaves, *A proof of the strong Gaussian product inequality
conjecture* (July 2026 preprint), Corollary 4.7, prove a lower bound for the
hafnian of a **repeated-index PSD covariance matrix**, with equality
characterized by a zero row or diagonality. That result concerns a single
repeated-index hafnian and a lower bound. Equations (3)-(9) here concern upper
bounds for averaged squared **principal** subhafnians, and therefore address a
different quantity and direction.

Kummer's ordinary Kneser adjacency construction in *Spectral linear matrix
inequalities* and Feinsilver's Boolean-lattice `sl(2)` representation are
relevant background. Neither source was found to state (3) or (5). A final
specialist priority review remains required before claiming first discovery.

Primary sources checked during this pass:
- B. Roos, arXiv:1906.06176v2 / DOI 10.1080/03081087.2024.2436061.
- F. Ouimet and D. Greaves, *A proof of the strong Gaussian product inequality
  conjecture*, July 2026 preprint, especially Corollary 4.7.
- M. Kummer, arXiv:2008.13452 / Adv. Math. 384 (2021), 107749.
- P. Feinsilver, arXiv:1102.0368.
- C. S. Hamilton et al., Phys. Rev. Lett. 119, 170501 (2017).
- J. M. Arrazola and T. R. Bromley, Phys. Rev. Lett. 121, 030503 (2018).

## 8. Verification

`tests/test_subhafnian_hierarchy.py` checks:

- the exact middle-operator recursion `T_m v = r h` on integer matrices;
- the critical equality families through order eight;
- the global double-counting multiplicities before the interlacing step;
- random real symmetric instances of (3) through `r=4`;
- the hybrid bound against both pure chains;
- the crossover formulas, one-switch hybrid form, and sharp edge-density bound;
- real-interferometer GBS-like PSD matrices, where the common collision-free
  normalization cancels from the sector recurrence.

`src/no_free_collapse/subhafnian_hierarchy.py` computes all principal hafnians
for small matrices by a memoized `2^N` dynamic program. It is a diagnostic
implementation; the analytic proof above is the universal certificate.
