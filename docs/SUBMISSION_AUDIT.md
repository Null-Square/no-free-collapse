# Submission-candidate audit

**12 September 2026.** This is a self-audit with additional derivations and tests, not an independent human referee report.

## Outcome

The complementary-subset dilation, critical recursion and global double-counting argument survived this pass. No counterexample to the stated real spectral hierarchy was found. That numerical statement is not its proof; the manuscript contains the analytic argument. The manuscript has been narrowed to a focused, self-contained claim set.

## 1. The missing odd-dimension Roos justification

Roos's Theorem 2.3 explicitly assumes even ambient order. It was therefore insufficient to cite that statement alone for arbitrary odd N. His Corollary 3.1 is parity-free. Set `F_j(S)=haf(A[S])/(2j-1)!!`. For `r=a+b`, the matching-partition identity gives

\[
\sum_{T\subset S,\ |T|=2a}\operatorname{haf}(A[T])\operatorname{haf}(A[S\setminus T])
=\binom ra\operatorname{haf}(A[S]).
\]

The factorial identity converts this to

\[
F_r(S)=\binom{2r}{2a}^{-1}\sum_{T\subset S,\ |T|=2a}F_a(T)F_b(S\setminus T).
\]

Apply Roos Corollary 3.1 with subset size 2r and block sizes (2a,2b). Its conclusion is `G_r<=G_a G_b`, without a parity restriction on N. This closes the citation gap; it does not claim a new replacement for Roos's theorem. The main manuscript now supplies the derivation. Primary source: https://arxiv.org/html/1906.06176v2 , Corollary 3.1.

## 2. Sharp density in every fixed dimension

For non-scalar real symmetric A of order N>=2, put `B=(A-lambda_min I)/d(A)` and `theta=4G_1(A)/d(A)^2`. Then

\[
\boxed{\theta\le\frac{4\lfloor N^2/4\rfloor}{N^2(N-1)}}.
\]

The right side is 1/(N-1) for even N and (N+1)/N^2 for odd N. If the eigenvalues of B are lambda_i in [0,1], Cauchy on the diagonal gives

\[
2H_2(B)\le\sum_i\lambda_i^2-\frac{(\sum_i\lambda_i)^2}{N}.
\]

The last expression is strictly convex separately in each coordinate. Its maximum is attained only at a cube vertex with k ones, where its value is k(N-k)/N. The maximizing ranks are floor(N/2) and ceil(N/2). Equality in the diagonal Cauchy bound requires every diagonal entry k/N. Thus equality is exactly a constant-diagonal projection of either maximizing rank. Matching projections give existence for even N; real Fourier sine/cosine pairs, with the constant column for odd rank, give existence for odd N. This elementary sharpening is supporting structure, not the headline novelty claim.

## 3. Boundary and presentation corrections

All recurrence statements require integer r>=2 and N>=4r-2. Division by the spectral diameter is restricted to non-scalar matrices; scalar matrices have zero positive-order principal hafnians. If G_1=0 the hybrid is handled before dividing by an upper bound. The PSD real-interferometer GBS subclass explicitly requires nonnegative squeezing parameters, a pure zero-displacement state and collision-free number-resolving patterns.

The six-variable equality proof now spells out support propagation, the three-edge matching alternative and the sign consistency argument. Its finite graph regression checks all 2^15 unweighted supports, while the weighted theorem remains analytic.

## 4. Implementation hardening

The floating principal-hafnian routine now returns a read-only matrix-bound cache and enforces an order-20 cap. Incomplete dictionaries and caches for a different off-diagonal matrix are rejected. Floating-point routines are labeled as such. A separate Fraction-only order-14 certificate verifies both spectral interval hypotheses and the averaged inequality exactly. Noninteger orders, bools, unsupported sizes, nonsymmetry and invalid spectral intervals receive explicit errors.

## 5. Remaining review limits

The full Lindell--Staples norm paper was not accessible; only its abstract and publisher preview were examined. The comparison with Roos is proposition-level, while the survey of all possible weighted subset operators is not exhaustive. Human verification of the proofs and a specialist priority check remain advisable before submission. There is no external review certificate and no claim that journal acceptance follows from test results.

The separate cube-normalized 1/216 conjecture, the complete critical equality locus for r>=3, and optimal fixed-N hierarchy constants beyond critical dimension remain open. None is assumed by the primary paper.
