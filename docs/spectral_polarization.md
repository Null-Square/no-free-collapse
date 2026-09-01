# Nested spectral polarization reduction

This note gives an exact bridge from the completed projection-gradient theorem to the still-open PSD-contraction theorem.

For a real symmetric `6 x 6` matrix `A`, define

\[
q_1(A)=\sum_{i<j}A_{ij}^2,
\qquad
q_2(A)=\sum_{i<j}\operatorname{haf}(A_{\widehat i,\widehat j})^2,
\]

and

\[
\Delta(A)=\frac14q_1(A)-q_2(A).
\]

The projection theorem now proves `Delta(P)>=0` for every orthogonal projection `P`.  The remaining target is

\[
\Delta(A)\ge0
\qquad\text{for every }0\preceq A\preceq I.
\tag{1}
\]

A naive extreme-point argument does not work: for a fixed eigenbasis, varying one eigenvalue can give a strict interior maximum of `q2-q1/4`. Thus separate convexity / endpoint maximization is not available.

The correct exact reduction instead uses a symmetric four-slot kernel on the nested spectral projections of `A`.  A terminology point is important: because `Delta` contains both a quadratic and a quartic part, this kernel is obtained by homogenizing the quadratic part with the spectral barycentric identity `sum delta_k=1`.  It is therefore not literally four-linear in the projection matrices themselves.

## 1. Square-free product polarization

Let commuting square-free generators satisfy `z_i^2=0`, and put

\[
\Omega_A=\sum_{i<j}A_{ij}z_i z_j.
\]

For two symmetric matrices `A,B`, the degree-four coefficient vector of

\[
\Omega_A\Omega_B
\]

is bilinear in `(A,B)`.  When `A=B`,

\[
\|\Omega_A^2\|^2=4q_2(A).
\]

For four projections `P_1,...,P_4`, define the symmetric homogenized defect kernel

\[
\boxed{
\mathcal D(P_1,P_2,P_3,P_4)
=
\frac1{24}\sum_{a<b}
\langle\Omega_{P_a},\Omega_{P_b}\rangle
-
\frac1{12}\sum_{\text{3 pairings}}
\langle
\Omega_{P_a}\Omega_{P_b},
\Omega_{P_c}\Omega_{P_d}
\rangle.
}
\tag{2}
\]

The first sum runs over the six pairs among the four arguments. The second sum runs over the three pairings

\[
(12)(34),\qquad(13)(24),\qquad(14)(23).
\]

Setting all four arguments equal gives exactly

\[
\mathcal D(P,P,P,P)=\Delta(P).
\tag{3}
\]

Thus the completed projection theorem is the diagonal specialization of `D`.  Equation (2) is symmetric in its four slots, but the first term is only pairwise linear in the projection matrices; the degree-four structure appears after homogenization in the spectral weights below.  In particular, one must not expand `D` termwise under a decomposition of a projection into rank-one summands as though `D` were matrix-multilinear.

## 2. Spectral chain

Let `A` be a nonzero PSD contraction and put

\[
\lambda=\lambda_{\max}(A),
\qquad
B=A/\lambda.
\]

Since `lambda<=1`, proving `Delta(B)>=0` is enough: indeed

\[
q_2(A)=\lambda^4q_2(B)
\le\frac{\lambda^4}{4}q_1(B)
\le\frac{\lambda^2}{4}q_1(B)
=\frac14q_1(A).
\]

So normalize `lambda_max(B)=1`. Order the eigenvalues

\[
1=\mu_1\ge\mu_2\ge\cdots\ge\mu_6\ge0,
\]

and let `P_k` be the projection onto the first `k` eigenvectors. Define

\[
\delta_k=\mu_k-\mu_{k+1},
\qquad \mu_7=0.
\]

Then

\[
\boxed{
B=\sum_{k=1}^6\delta_kP_k,
\qquad
\delta_k\ge0,
\qquad
\sum_k\delta_k=1,
}
\tag{4}
\]

and

\[
P_1\le P_2\le\cdots\le P_6=I.
\]

This is the finite spectral version of the layer-cake formula

\[
B=\int_0^1 P_t\,dt,
\qquad
P_t=1_{B\ge t}.
\]

The spectral decomposition/layer-cake representation is standard functional calculus and is not a novelty claim.

## 3. Exact degree-four reconstruction

Because the first part of `Delta` is quadratic and the second part quartic, homogenize the quadratic term with

\[
\left(\sum_k\delta_k\right)^2=1.
\]

Expanding (4) gives the exact identity

\[
\boxed{
\Delta(B)
=
\sum_{i,j,k,l=1}^6
\delta_i\delta_j\delta_k\delta_l
\,\mathcal D(P_i,P_j,P_k,P_l).
}
\tag{5}
\]

Equivalently, in layer-cake form,

\[
\Delta(B)
=
\int_{[0,1]^4}
\mathcal D(P_{t_1},P_{t_2},P_{t_3},P_{t_4})
\,dt_1dt_2dt_3dt_4.
\tag{6}
\]

Since the four threshold projections are automatically nested after sorting the four thresholds, the contraction theorem follows from one precise mixed statement:

\[
\boxed{
P_1\le P_2\le P_3\le P_4
\quad\Longrightarrow\quad
\mathcal D(P_1,P_2,P_3,P_4)\ge0.
}
\tag{N}
\]

Equation (N) is **not yet proved**. It is the next sharp theorem target.

## 4. Finite coefficient form

For a fixed eigenbasis, (5) is a homogeneous degree-four polynomial in the six nonnegative spectral increments `delta_k`. There are

\[
\binom{6+4-1}{4}=126
\]

multinomial coefficients. Each normalized coefficient is exactly `D` evaluated on one nondecreasing rank quadruple

\[
1\le r_1\le r_2\le r_3\le r_4\le6.
\]

Thus the continuum problem has become a finite family of mixed nested-projection inequalities, uniform over the orthogonal eigenbasis.

Several simplifications are immediate:

- `r_4=6` is trivial because `P_6=I` has zero off-diagonal / zeon degree-two component, so every quartic mixed term vanishes and only a nonnegative pair-inner-product term remains;
- simultaneous complementation preserves `D` (off-diagonal signs flip in every projection, but all terms have even total degree);
- `r_1=r_2=r_3=r_4` is exactly the already-proved projection theorem.

The genuinely new content is therefore only the mixed nested ranks.

## 5. Numerical diagnostic

Deterministic random tests currently find all mixed nested coefficients nonnegative to floating-point precision. In a direct scan of 200 random orthogonal eigenbases, all 126 coefficients were nonnegative; the worst values were on the order of `5e-17`, consistent with exact zeros from trivial identity-containing patterns.

This is evidence only. The repository does not promote (N) or the PSD-contraction theorem to a theorem until an analytic or exact finite certificate is obtained.

A separate tempting route has already been rejected: the defect is not separately convex in the eigenvalues, and an interior eigenvalue can give a larger value of `q2-q1/4` than both endpoints. The nested homogenized kernel inequality is therefore the preferred next direction.

## 6. Consequence if (N) is proved

If (N) holds, then (5) immediately gives

\[
q_2(A)\le\frac14q_1(A)
\qquad(0\preceq A\preceq I).
\]

Scaling an arbitrary PSD matrix by `lambda_max(A)` yields

\[
\boxed{
q_2(A)
\le
\frac{\lambda_{\max}(A)^2}{4}
q_1(A).
}
\tag{7}
\]

Euler plus Cauchy would then give the prospective six-variable spectral hafnian-energy inequality, feeding directly back into the sharp PSD hafnian-capacity problem.
