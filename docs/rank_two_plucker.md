# Exact rank-two Pluecker reduction

Let `P` be a real `6 x 6` rank-two orthogonal projection.  Choose a `6 x 2` Parseval frame `U` with `P=UU^T`, and write its row vectors as `u_i in R^2`.

Define

\[
d_i=P_{ii}=\|u_i\|^2,
\qquad
m_{ij}=\det(u_i,u_j)^2.
\]

By Lagrange's identity,

\[
\boxed{m_{ij}=d_i d_j-P_{ij}^2.}
\tag{1}
\]

The `m_ij` are the squared Pluecker coordinates of the underlying two-plane.

## Edge-probability structure

Because `U` is Parseval,

\[
\sum_{j\ne i}m_{ij}=d_i,
\tag{2}
\]

and Cauchy-Binet gives

\[
\sum_{i<j}m_{ij}=1.
\tag{3}
\]

Thus `m` is a probability distribution on the 15 edges of `K_6`, with vertex marginals `d_i`.

The principal-minor inequalities for `P` and `I-P` also give

\[
\max(0,d_i+d_j-1)\le m_{ij}\le d_i d_j.
\tag{4}
\]

These are necessary but not sufficient for decomposability; the signed square roots of `m_ij` additionally satisfy the quadratic Pluecker relations.

## Tensor-Parseval identity

Let

\[
q_1(P)=\sum_{i<j}P_{ij}^2,
\qquad
q_2(P)=\sum_{i<j}\operatorname{haf}(P_{\widehat i,\widehat j})^2.
\]

Since `P` is a projection of trace two,

\[
q_1=1-\frac12D_2,
\qquad
D_k:=\sum_i d_i^k.
\tag{5}
\]

For the underlying two-dimensional Parseval frame, the ordered fourth Gaussian moment tensor has squared norm `24`.  Partitioning the ordered index quadruples according to collision type `4`, `3+1`, `2+2`, `2+1+1`, and `1+1+1+1`, and simplifying with projection idempotence, gives

\[
\boxed{
q_2
=1-\frac52D_2+3D_3-\frac98D_4
+\frac14\sum_{i<j}(d_i d_j+2P_{ij}^2)^2.
}
\tag{6}
\]

This identity has been checked independently against direct complementary-hafnian evaluation.

## Pure Pluecker defect identity

Using `P_ij^2=d_i d_j-m_ij`, define

\[
C=\sum_{i<j}d_i d_jm_{ij},
\qquad
M_2=\sum_{i<j}m_{ij}^2.
\]

Then (5)--(6) reduce the desired projection-gradient defect to

\[
\boxed{
8\left(\frac14q_1-q_2\right)
=-6+19D_2-24D_3+18D_4-9D_2^2+24C-8M_2.
}
\tag{7}
\]

So the global rank-two problem is no longer a hafnian problem: it is a quartic inequality for the squared Pluecker edge distribution of a decomposable unit bivector.

## Concavity and chord linearization

For fixed diagonal `d`, the Pluecker-dependent part is

\[
24C-8M_2
=\sum_{i<j}8m_{ij}(3d_i d_j-m_{ij}).
\]

Each edge term is concave and increasing on the admissible interval (4).  Writing

\[
a_{ij}=d_i d_j,
\qquad
\ell_{ij}=\max(0,d_i+d_j-1),
\]

its chord gives the rigorous linear lower bound

\[
8m(3a-m)
\ge
8\left[a\ell+(2a-\ell)m\right],
\qquad \ell\le m\le a.
\tag{8}
\]

Large random searches over rank-two projections found the resulting linearized defect nonnegative, nearly tight only near coordinate projections.  Equation (8) is therefore a concrete route to the remaining global theorem, but no global proof from the marginal constraints alone is claimed here.

## Two-direction extremizer

When the planar frame uses only two orthogonal directions, the Pluecker distribution is supported on a complete bipartite graph between the two coordinate supports.  The companion note [`rank_two_two_direction.md`](rank_two_two_direction.md) solves that class sharply: the unique nondegenerate maximizer has support partition `2+4` with equal weights and

\[
\frac{q_2}{q_1}=\frac{33}{160}<\frac14.
\]

Therefore any obstruction to the global rank-two `1/4` contraction must use a genuinely multi-directional planar frame and the nontrivial Pluecker relations.
