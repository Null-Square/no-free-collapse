# Sharp equal-diagonal rank-three projection theorem

This note closes the sharp hafnian-gradient inequality on the critical
equal-diagonal rank-three projection stratum.

The general projection-gradient notation is defined in
[`projection_gradient.md`](projection_gradient.md):

\[
q_1(A)=\sum_{i<j}A_{ij}^2,
\qquad
q_2(A)=\sum_{i<j}\operatorname{haf}(A_{\widehat i,\widehat j})^2.
\]

The target for an orthogonal projection `P` is

\[
q_2(P)\le\frac14q_1(P).
\]

## Theorem

Let `P` be a real `6 x 6` orthogonal projection of rank three with constant
diagonal

\[
P_{11}=\cdots=P_{66}=\frac12.
\]

Then

\[
\boxed{q_2(P)\le\frac14q_1(P).}
\]

In fact

\[
q_1(P)=\frac34
\]

and, writing `K=2P-I`,

\[
\boxed{
16q_2(P)
=q_2(K)
=\sum_{i<j}K_{ij}^4.
}
\tag{1}
\]

Therefore

\[
q_2(P)\le\frac3{16}.
\]

Equality holds exactly when the support of the off-diagonal part of `K` is a
perfect matching and every nonzero entry has absolute value one.  Equivalently,
up to coordinate permutation and independent sign switches, `P` is the three
block disjoint-pair projection

\[
\frac12\operatorname{diag}(J_2,J_2,J_2).
\]

Thus the disjoint-pair equality geometry is the unique equality type on the
entire equal-diagonal rank-three projection stratum.

## Reduction to a zero-diagonal involution

Set

\[
K=2P-I.
\]

Since `P` is a rank-three orthogonal projection,

\[
K^T=K,
\qquad
K^2=I,
\qquad
\operatorname{tr}K=0.
\]

The constant diagonal of `P` gives

\[
K_{ii}=0
\]

for every `i`.  Conversely, every symmetric orthogonal zero-diagonal `6 x 6`
matrix gives such a projection through `P=(I+K)/2`.

For `i != j`, `P_ij=K_ij/2`, so

\[
q_1(P)=\frac14q_1(K),
\qquad
q_2(P)=\frac1{16}q_2(K).
\]

It is therefore enough to prove

\[
q_2(K)=\sum_{i<j}K_{ij}^4\le q_1(K).
\]

## Edge-square bookkeeping

Put

\[
x_{ij}=K_{ij}^2.
\]

Because the diagonal of `K` vanishes and every row of the orthogonal matrix
`K` has Euclidean norm one,

\[
\sum_{j\ne i}x_{ij}=1
\tag{2}
\]

for every vertex `i`.  Define

\[
A=\sum_{i<j}x_{ij}^2,
\]

let `B` be the sum of `x_e x_f` over unordered pairs of adjacent edges, and
let `C` be the corresponding sum over unordered pairs of disjoint edges.

Squaring (2) and summing over the six vertices gives

\[
6=2A+2B,
\]

hence

\[
A+B=3.
\tag{3}
\]

Also

\[
\sum_{i<j}x_{ij}=3,
\]

so squaring this identity gives

\[
9=A+2B+2C.
\]

Using (3),

\[
C=\frac{3+A}{2}.
\tag{4}
\]

## Orthogonality controls the four-cycle cross terms

For distinct vertices `i,j`, row orthogonality and zero diagonal imply

\[
\sum_{k\ne i,j}K_{ik}K_{jk}=0.
\]

Squaring,

\[
0=
\sum_{k\ne i,j}x_{ik}x_{jk}
+2\sum_{\{k,l\}\subset[6]\setminus\{i,j\}}
K_{ik}K_{jk}K_{il}K_{jl}.
\tag{5}
\]

Sum (5) over all unordered pairs `{i,j}`.  The first term counts every pair of
adjacent edges exactly once, and therefore contributes `B`.

Let `D` denote the sum of the signed products around the three four-cycles on
each four-vertex subset.  Every such four-cycle product appears in (5) for
its two pairs of opposite vertices, and each occurrence carries the factor
`2`.  Consequently

\[
0=B+4D,
\qquad
D=-\frac B4.
\tag{6}
\]

## Exact gradient-energy identity

For an omitted edge `e`, the complementary `4 x 4` hafnian is the sum of the
three perfect-matching products on the remaining four vertices.  On squaring
and summing over all omitted edges:

- the squares of the three matching terms sum to `C`, because every unordered
  pair of disjoint edges determines a unique omitted complementary edge;
- the pairwise cross terms give `2D`.

Thus

\[
q_2(K)=C+2D.
\]

Using (3), (4), and (6),

\[
q_2(K)
=\frac{3+A}{2}-\frac B2
=\frac{3+A-(3-A)}2
=A.
\]

This proves the exact identity

\[
\boxed{q_2(K)=\sum_{i<j}K_{ij}^4.}
\]

## Sharp bound and equality

Equation (2) implies

\[
q_1(K)=\sum_{i<j}x_{ij}=3.
\]

Moreover `0<=x_ij<=1`, so

\[
q_2(K)=\sum_{i<j}x_{ij}^2
\le
\sum_{i<j}x_{ij}
=3.
\]

Scaling back to `P` yields

\[
q_2(P)\le\frac3{16}=\frac14q_1(P).
\]

Equality requires `x_ij^2=x_ij` for every edge, hence every `x_ij` is zero or
one.  The row constraints (2) then force exactly one unit edge incident to
each vertex.  The support is therefore a perfect matching.  Conversely every
signed perfect matching gives equality.

## Significance

The rank-three projection problem is still open away from constant diagonal,
but the stratum containing the sharp disjoint-pair equality point is now
solved globally, with a complete equality characterization.  Any rank-three
counterexample to the proposed gradient contraction must therefore have
nonconstant diagonal.
