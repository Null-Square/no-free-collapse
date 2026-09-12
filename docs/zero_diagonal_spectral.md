# Sharp spectral subhafnian bounds and matching stability

**Status: analytic theorems, with exact polynomial and rational checks.**
This note does not settle the unrestricted PSD-contraction conjecture or the
Boolean-cube-normalized `1/216` conjecture. The zero-diagonal hypothesis applies
to the matrix whose operator norm is bounded; a diagonal cannot be removed for
free when transferring a norm bound.

## 1. Definitions and sharp bound

Let `B=B^T` be a complex symmetric `n x n` matrix, `n>=4`, with zero diagonal.
The transpose in this assumption is not the conjugate transpose. Put

\[
T=\sum_{i<j}|b_{ij}|^2,\qquad
H_4=\sum_{|S|=4}|\operatorname{haf}(B[S])|^2,
\]

\[
w_{ij}=|b_{ij}|^2,\quad r_i=\sum_{j\ne i}w_{ij},\quad
W=\sum_i\sum_{\substack{j<k\\j,k\ne i}}w_{ij}w_{ik},\quad G=BB^*.
\]

For any `ell>=||B||_op`,

\[
\boxed{H_4\le\frac{n-2}{4}\ell^2T.}\tag{1}
\]

The constant is optimal for each even `n>=4`. For `ell>0`, equality holds
exactly for `B=0` and for permutation-conjugates of direct sums of blocks

\[
\begin{pmatrix}0&\ell\omega_j\\\ell\omega_j&0\end{pmatrix},
\qquad |\omega_j|=1.
\]

In particular, odd orders have no nonzero equality case. This statement does
not claim that the constant is optimal in odd order. In the real case the
phases reduce to signs.

## 2. Exact four-term defect identity

The proof is the identity

\[
\boxed{\begin{aligned}
\frac{n-2}{4}\ell^2T-H_4
={}&\frac{n-4}{2n}T\left(\frac{n\ell^2}{2}-T\right)
+\frac12\sum_i\left(r_i-\frac{2T}{n}\right)^2\\
&+W+\frac14\operatorname{tr}\bigl(G(\ell^2I-G)\bigr).
\end{aligned}}\tag{2}
\]

Every term is nonnegative: `T=tr(G)/2<=n ell^2/2`, `G` has eigenvalues in
`[0,ell^2]`, and the other terms are sums of nonnegative real quantities.

Here is a direct derivation. Let `F=sum_{i<j}w_ij^2`, let `D` be the sum of
`w_e w_f` over unordered disjoint edges, and let `C` be the sum, over four-sets,
of the three pairwise products of matching monomials, conjugating the second
factor. Expanding the squared hafnians and the off-diagonal entries of `G`
gives

\[
H_4=D+2\operatorname{Re}C,\qquad
\sum_{i<j}|G_{ij}|^2=W+4\operatorname{Re}C.
\]

Every wedge is counted once in the second identity; every matching cross-term
is counted through two opposite pairs. Also

\[
T^2=F+2D+2W,\qquad \sum_i r_i^2=2F+2W.
\]

Eliminating `C,D,W` yields the trace identity

\[
\boxed{H_4=\frac12T^2+F-\sum_i r_i^2+\frac14\operatorname{tr}(G^2).}\tag{3}
\]

Substitute `F=(sum r_i^2)/2-W` and
`sum r_i^2=4T^2/n+sum(r_i-2T/n)^2` to obtain (2).

For equality in (1), `W=0` means each row has at most one nonzero entry.
The row-variance term forces every row to have the same squared norm. Unless
`B=0`, every vertex therefore has one partner, so `n` is even and the support
is a perfect matching. The spectral term then forces each edge modulus to be
`ell`. This argument also covers `n=4`, where the first term of (2) vanishes.

## 3. Quantitative stability in all even orders at least six

Normalize `||B||_op<=1` and write `delta=(n-2)T/4-H_4`.
Assume

\[
T\ge n/4,\qquad
\delta\le\min\{(n-4)/64,1/16\}.
\]

Then a unit-modulus perfect-matching matrix `M` can be recovered by selecting
the largest squared entry of every row and retaining its phase, and

\[
\boxed{\|B-M\|_F^2\le\left(\frac{16}{n-4}+\frac{16}{3}\right)\delta.}\tag{4}
\]

**Proof.** Let `E=n-2T=sum_i(1-r_i)`. The first term of (2), together with
`T>=n/4`, gives `E<=16 delta/(n-4)<=1/4`. Each `r_i<=1`, so each `r_i>=3/4`.
Let `W_i` be the row contribution to `W` and `m_i=max_j w_ij`. Then

\[
m_i\ge\frac{\sum_jw_{ij}^2}{r_i}
=r_i-\frac{2W_i}{r_i}
\ge\frac34-\frac83\delta\ge\frac7{12}>\frac12.
\]

Thus dominant edges are mutual: a row of squared norm at most one cannot
contain two entries of squared modulus greater than one half. They form a
perfect matching. The squared distance contributed by row `i` is

\[
r_i-m_i+(1-\sqrt{m_i})^2
\le1-r_i+2(r_i-m_i)
\le1-r_i+4W_i/r_i.
\]

Summing proves (4). If instead `T<n/4`, the first term of (2) gives
`||B||_F^2=2T<=16 delta/(n-4)`. Thus a small defect gives an explicit dichotomy:
proximity to zero or proximity to a unit-modulus matching. A nondegeneracy
condition cannot be omitted when demanding proximity specifically to a unit
matching.

## 4. A sharp six-variable hafnian theorem

For any zero-diagonal complex symmetric `6 x 6` matrix,

\[
\boxed{|\operatorname{haf}(B)|\le\|B\|_{op}^3.}\tag{5}
\]

Euler's homogeneous identity and Cauchy-Schwarz give

\[
3|\operatorname{haf}(B)|
=\left|\sum_{i<j}b_{ij}\operatorname{haf}(B_{\widehat i,\widehat j})\right|
\le\sqrt{T H_4}\le\ell T\le3\ell^3.
\]

For nonzero `B`, equality holds exactly for the matching matrices in Section 1,
with `ell=||B||_op`.

There is also full-hafnian stability. If `||B||_op<=1` and

\[
|\operatorname{haf}(B)|\ge1-\varepsilon,
\qquad 0\le\varepsilon\le1/48,
\]

then the same rowwise phase-preserving rounding gives

\[
\boxed{\|B-M\|_F^2\le34\varepsilon.}\tag{6}
\]

Indeed `T>=3(1-epsilon)`, so `E=6-2T<=6 epsilon`, and

\[
\delta=T-H_4\le T-9(1-\varepsilon)^2/T
\le3-3(1-\varepsilon)^2\le6\varepsilon.
\]

Consequently `r_i>=1-6 epsilon>=7/8` and
`m_i>=1-6 epsilon-12 epsilon/(1-6 epsilon)>=33/56>1/2`.
The distance estimate from Section 3 is at most

\[
6\varepsilon+\frac{24\varepsilon}{1-6\varepsilon}
\le\frac{234}{7}\varepsilon<34\varepsilon.
\]

## 5. What this adds to the existing project

For a real PSD contraction `0<=A<=I` with **every `A_ii=1/2`**, set `B=2A-I`.
Then `B` is zero-diagonal and `||B||_op<=1`. Since `T(B)=4q1(A)` and
`H_4(B)=16q2(A)`, the six-variable conclusion is

\[
q_2(A)\le q_1(A)/4.
\]

This covers non-idempotent contractions, not only the equal-diagonal rank-three
projections previously handled in the repository. The new result also provides
all-even-dimensional complex-symmetric bounds and explicit stability.

It is **not** the general PSD-contraction theorem. For instance, for the
rank-one projection `P=J_6/6`, `||2P-I||=1`, but
`||2 offdiag(P)||=5/3`. Deleting the diagonal preserves hafnians but invalidates
the required norm hypothesis.

## 6. Prior-art comparison and limits of the novelty claim

Roos, *New inequalities for permanents and hafnians and some generalizations*,
Linear and Multilinear Algebra 73(8), 1634-1667 (2025), Theorem 2.3, implies
for even `n`

\[
H_4\le\alpha_n T^2,\qquad
\alpha_n=\frac{3(n-2)(n-3)}{2n(n-1)}.
\]

The spectral bound (1) is smaller whenever
`T>n(n-1) ell^2/[6(n-3)]`. At a unit matching it is exact; the displayed Roos
specialization is larger by a factor `3(n-3)/(n-1)`. At `n=6`, the comparison
is `H_4<=3` versus `H_4<=27/5`. Conversely, the constant off-diagonal matrix
`(J-I)/(n-1)` attains Roos's bound, so neither displayed bound uniformly
supersedes the other.

General zeon norm inequalities are established prior art (Lindell and Staples,
2019). The contribution proposed here is the spectral-sensitive defect identity,
its equality classification, and the explicit constructive stability results.
The full text of the latter paper was not available in this audit. A specialist
comparison with that work and adjacent hafnian inequalities remains a
pre-submission novelty check; absence from a search is not proof of priority.

Primary sources checked:
- https://arxiv.org/html/1906.06176v2 (Roos, especially Theorem 2.3).
- https://doi.org/10.1080/03081087.2024.2436061 (published Roos article).
- https://doi.org/10.1007/s00006-018-0934-z (Lindell and Staples, abstract and metadata).

## 7. Reproducibility and status

Run `python -m pytest tests/test_zero_diagonal.py`. Exact sparse-polynomial tests
compare all coefficients of the trace and defect identities for independent real
and imaginary entries. Rational-instance certificates verify both spectral
hypotheses by exact semidefinite elimination, then check nonnegativity and zero
identity residual. Floating-point checks are separately labeled diagnostics.

Run `python experiments/e20_spectral_stability.py --samples 250 --output results/spectral_stability.json`
after installing the package (or set `PYTHONPATH=src`). The script records its seed,
environment, sharp certificates, matching-recovery errors, and numerical identity
residuals. Numerical sampling is not used as a proof of any universal inequality.

The complete analytic proofs above are the research result. Independent proof
review, a full novelty comparison, author metadata, and venue-specific formatting
remain before external submission. Neither repository CI nor a numerical search
can substitute for those publication checks.
