# Complementary-subset dilation: the full PSD gradient theorem

**Status: proved analytically.** This note closes results-ledger E4 and E5.
It also proves the stronger six-variable PSD operator bound that older notes
incorrectly labeled false. The Boolean-cube-normalized `1/216` extremal problem
is a different constraint problem and remains open in this work.

## 1. A sharp theorem in every even dimension

Let `n=2m`, `m>=2`, and index a matrix `T_m(A)` by the `(m-1)`-subsets of `[n]`.
For a real symmetric matrix `A`, set its `(E,F)` entry to zero when `E,F`
intersect; otherwise set it to `A_ij`, where `{i,j}` is the complement of their union.
Then

\[
\boxed{\|T_m(A)\|_{\rm op}\le\frac{m+1}{4}
(\lambda_{\max}(A)-\lambda_{\min}(A)).}\tag{1}
\]

The coefficient is sharp for every `m>=2`. In particular,

\[
0\preceq A\preceq LI\quad\Longrightarrow\quad
-\frac{m+1}{4}LI\preceq T_m(A)\preceq\frac{m+1}{4}LI.\tag{2}
\]

### Explicit dilation and proof

For an `m`-subset `R` and vector `v` indexed by `(m-1)`-subsets define

\[
(c_R(v))_i=\begin{cases}v_{R\setminus\{i\}},&i\in R,\\0,&i\notin R.\end{cases}
\]

Choose one representative from each pair `{R,R^c}`. The number of pairs is
`N=binom(2m,m)/2`. Define rectangular maps `L_+` and `L_-` by stacking
`c_R(v)+c_Rc(v)` and `c_R(v)-c_Rc(v)`, respectively. Each row has one nonzero
coefficient, equal to `+1` or `-1`. Each input coordinate appears `m+1` times.
Consequently

\[
L_+^TL_+=L_-^TL_-=(m+1)I.\tag{3}
\]

Expansion of the quadratic forms gives the exact identity

\[
\boxed{T_m(A)=\frac14\left[
L_+^T(I_N\otimes A)L_+-L_-^T(I_N\otimes A)L_-
\right].}\tag{4}
\]

Indeed the difference for one complementary pair is `4 c_R(v)^T A c_Rc(v)`.
Summing these cross terms counts each unordered disjoint pair `E,F` twice,
exactly as in `v^T T_m(A)v`: the two complementary partitions are obtained by
assigning either of the two remaining vertices to `E`.

Put `c=(m+1)/4`. For `0<=A<=LI`, (3)-(4) give the two explicit PSD certificates

\[
\begin{aligned}
cLI-T_m(A)&=\tfrac14\bigl[
L_+^T(I_N\otimes(LI-A))L_++L_-^T(I_N\otimes A)L_-\bigr],\\
cLI+T_m(A)&=\tfrac14\bigl[
L_-^T(I_N\otimes(LI-A))L_-+L_+^T(I_N\otimes A)L_+\bigr].
\end{aligned}\tag{5}
\]

This proves (2). Since diagonal entries of `A` do not enter `T_m`, subtracting
`lambda_min(A) I` proves (1). For sharpness, take `A=J_n/n`. Each row of
`T_m(A)` has `binom(m+1,2)` nonzero entries, each `1/n`. The all-ones vector
therefore has eigenvalue `(m+1)/4`, while the spectral diameter of `A` is one.
No numerical optimization enters this proof.

## 2. The six-variable global hafnian-gradient theorem

For `n=6`, `T_3(A)` is exactly the repository's 15-edge matching operator. Put

\[
a_{ij}=A_{ij},\qquad h_{ij}=\operatorname{haf}(A_{\widehat i,\widehat j}),
\quad q_1=\sum_{i<j}a_{ij}^2,\quad q_2=\sum_{i<j}h_{ij}^2.
\]

The elementary identities `T_3(A)a=2h` and `3 haf(A)=a^T h` yield, with
`d= lambda_max(A)-lambda_min(A)`,

\[
\boxed{q_2(A)\le\frac{d^2}{4}q_1(A),\qquad
|\operatorname{haf}(A)|\le\frac d6q_1(A).}\tag{6}
\]

These hold for **every real symmetric** six-dimensional matrix. For PSD `A`,
`d<=lambda_max(A)`, proving the formerly open statements

\[
\boxed{q_2(A)\le\frac{\lambda_{\max}(A)^2}{4}q_1(A)}
\]

and, in particular, `q2(A)<=q1(A)/4` for every PSD contraction `0<=A<=I`.
The proof bypasses nested spectral coefficientwise positivity; it does not
assert that every coefficient of that stronger auxiliary conjecture is positive.

## 3. Complete equality classification for PSD matrices

Let `A` be nonzero PSD of order six and let `L=lambda_max(A)`.
Equality in `||T_3(A)||<=L` holds exactly when `A/L` has one of these forms:

\[
\boxed{\frac{ss^T}{6},\qquad I-\frac{ss^T}{6},\qquad\frac{I+M}{2},}\tag{7}
\]

where `s_i` are arbitrary signs and `M` is a signed perfect-matching matrix
(one signed unit edge in each row). The respective ranks are one, five, three.
For either `q2<=L^2 q1/4` or `|haf(A)|<=L q1/6`, the equality cases are (7)
and all diagonal PSD matrices, the latter having `q1=q2=haf=0`.

### Proof of necessity

Normalize `L=1` and choose a unit real vector `v` with
`v^T T_3(A)v=+1` or `-1`. Regard `v` as a symmetric zero-diagonal edge array.
Let `C` be the `6 x 20` matrix with columns `c_R(v)` for all triples, let `J`
swap each triple with its complement, and put `P_+=(I+J)/2`, `P_-=(I-J)/2`.
The matrix

\[
H=\frac12 CJC^T=H_+-H_-,\qquad H_\pm=\frac12 CP_\pm C^T
\]

has diagonal zero and off-diagonal entries `h_ij(v)`, the complementary
four-variable hafnians of the edge array `v`. Further,

\[
\operatorname{tr}H_+=\operatorname{tr}H_-=1,\qquad
v^T T_3(A)v=\operatorname{tr}(AH).
\]

To check the traces, `||C||_F^2=4||v||^2=4` and `tr(H)=0`.
Equality forces `A` to act as one on the range of `CP_+` and as zero on the
range of `CP_-`, or conversely. These ranges must be orthogonal. Thus
`C^T C` commutes with `J`.

The diagonal entries of this commutation relation give, for every triple `R`,

\[
\sum_{\{i,j\}\subset R}v_{ij}^2
=\sum_{\{i,j\}\subset R^c}v_{ij}^2.\tag{8}
\]

Using triples `{i,j,k}` and `{i,j,l}` gives

\[
v_{ik}v_{il}+v_{jk}v_{jl}=v_{pk}v_{pl}+v_{qk}v_{ql},
\]

where `p,q` are the remaining indices. Varying the split of the four indices
outside `{k,l}` implies

\[
\boxed{v_{ak}v_{al}\text{ is independent of }a\notin\{k,l\}.}\tag{9}
\]

If a vertex is adjacent to two others `k,l` in the support graph of `v`,
(9) forces both `k,l` adjacent to all four remaining vertices. Applying (9)
to pairs in those four vertices makes them a clique, and another application
forces the edge `k,l`. The support is therefore complete. Otherwise it is a
matching. In the matching case, (8), applied to a nonzero edge and each possible
third vertex, forces two further disjoint edges of the same squared weight.
Thus the support is a perfect matching with equal absolute edge weights.

In the complete case, fix `a,b` and put `r_k=v_ak/v_bk` for the other four
indices. Equation (9) gives `r_k r_l=1` for every distinct `k,l`; hence all
four ratios are the same sign. It follows that every incident edge has the
same modulus, and connectedness makes this modulus common to all edges.
The row-sign ratios and symmetry then give `v_ij=alpha s_i s_j` for `i!=j`,
with a common nonzero real `alpha` and signs `s_i`.

In the dense case,

\[
H=3\alpha^2\left(\prod_i s_i\right)D_s(J_6-I)D_s.
\]

Its positive and negative spectral subspaces have dimensions one and five,
in either order. In the matching case,

\[
H=c^2\left(\prod_{e\in M}\operatorname{sign}v_e\right)M,
\]

where `c` is the common edge modulus. This matrix has signatures three and
three. Both matrices are invertible. Therefore the orthogonal ranges of
`H_+,H_-` fill the entire space, and `A` must be precisely the positive or
negative spectral projection of `H`. These are exactly (7).

Conversely, direct substitution in all three families gives equality. For the
gradient inequality with `q1>0`, equality forces operator-norm equality via
`||T a||<=||T|| ||a||`, hence necessity follows. For the hafnian-energy bound,
equality similarly forces equality in the Rayleigh bound. When `q1=0`, `A`
is diagonal and equality is degenerate. This completes the classification.

## 4. Correction of the old operator-norm warning

The historical conference example has `K=S/sqrt(5)`, `S^2=5I`, and
`||T_3(K)||=3/sqrt(5)>1`. It refutes the bound `||T_3(K)||<=1` for symmetric
orthogonal involutions. It does **not** refute the PSD statement: for
`P=(I+K)/2`, one has

\[
\|T_3(P)\|=\tfrac12\|T_3(K)\|=\frac{3}{2\sqrt5}<1.
\]

The old results-ledger assertion that the PSD operator inequality itself was
false must therefore be withdrawn. Equation (5) is an explicit proof of it.
A regression test preserves both scaled quantities to prevent recurrence.

## 5. Scope, prior art, and verification

The result is a spectral theorem, not a proof of the separate maximum over
`B>=0` satisfying `max_x x^T B x<=1`. In particular, the cube-normalized
`1/216` conjecture and its proposed all-even-dimensional extension remain open.
The older rank-by-rank proofs remain valid but are no longer needed to prove
the common `1/4` bound. Sharp smaller constants within rank two remain distinct.

Kneser adjacency, square-free subset operators, and PSD compression are established
mathematical tools. Kummer, *Spectral linear matrix inequalities* (Adv. Math. 384,
107749, 2021), Lemma 3.5, uses an unweighted Kneser adjacency matrix to prove a
square-free derivative isomorphism; its main Theorem 2.2 concerns spectrahedral
representations. Those displayed statements are different from the weighted
spectral-diameter estimate and equality theorem here. Roos's averaged hafnian
inequalities and Lindell--Staples's zeon norm work are also relevant. This comparison
is not an exhaustive priority certification; specialist review remains appropriate.
Primary sources: https://arxiv.org/abs/2008.13452 and
https://arxiv.org/html/1906.06176v2 .

The implementation uses only NumPy and Python's standard library. Tests compare
(4) coefficient-by-coefficient over every symmetric matrix basis element for
orders 4,6,8,10, using integers. Rational certificates check both input PSD
hypotheses, both operator PSD slacks, and the corollaries exactly. Numerical
checks are separate diagnostics, never universal proofs. Explicit enumeration
is capped at order ten in code to avoid accidental exponential allocations;
there is no such cap in the analytic theorem.
