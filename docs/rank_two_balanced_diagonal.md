# Rank-two gradient theorem with balanced diagonals

This note proves the six-variable projection-gradient inequality on a genuinely multi-directional rank-two stratum.

Let `P` be a real `6 x 6` rank-two orthogonal projection and define

\[
q_1(P)=\sum_{i<j}P_{ij}^2,
\qquad
q_2(P)=\sum_{i<j}\operatorname{haf}(P_{\widehat i,\widehat j})^2.
\]

Assume

\[
\boxed{P_{ii}\le\frac12\quad\text{for every }i.}
\tag{1}
\]

Then

\[
\boxed{q_2(P)\le\frac14 q_1(P).}
\tag{2}
\]

No restriction is placed on the number of planar frame directions.

## 1. Exact Pluecker defect

Write

\[
d_i=P_{ii},
\qquad
D_k=\sum_i d_i^k,
\]

and let

\[
m_{ij}=d_i d_j-P_{ij}^2
\]

be the squared Pluecker coordinates.  From [`rank_two_plucker.md`](rank_two_plucker.md),

\[
8\left(\frac14q_1-q_2\right)
=-6+19D_2-24D_3+18D_4-9D_2^2+24C-8M_2,
\tag{3}
\]

where

\[
C=\sum_{i<j}d_i d_jm_{ij},
\qquad
M_2=\sum_{i<j}m_{ij}^2.
\]

Since `0<=m_ij<=d_i d_j`,

\[
M_2\le C,
\]

so

\[
24C-8M_2\ge16C.
\tag{4}
\]

## 2. A compression lower bound for `C`

Choose a `6 x 2` Parseval frame `U` with `P=UU^T` and put

\[
D=\operatorname{diag}(d_1,\ldots,d_6),
\qquad
A=U^TDU.
\]

Cauchy-Binet gives

\[
\det A
=\sum_{i<j}d_i d_j\det(u_i,u_j)^2
=C.
\tag{5}
\]

Also

\[
\operatorname{tr}A
=\operatorname{tr}(DP)
=D_2.
\]

Because `P<=I`,

\[
DPD\preceq D^2.
\]

Therefore

\[
\operatorname{tr}(A^2)
=\operatorname{tr}(DPDP)
\le\operatorname{tr}(D^2P)
=D_3.
\]

For a `2 x 2` matrix,

\[
\det A
=\frac{(\operatorname{tr}A)^2-\operatorname{tr}(A^2)}2,
\]

hence

\[
\boxed{
C\ge\frac{D_2^2-D_3}{2}.
}
\tag{6}
\]

Combining (3), (4), and (6), it is enough to prove

\[
F(d):=-6+19D_2-32D_3+18D_4-D_2^2\ge0.
\tag{7}
\]

## 3. Reduction to a scalar moment inequality

Under (1), define

\[
x_i=2d_i\in[0,1].
\]

Since `Tr(P)=2`,

\[
\sum_i x_i=4.
\]

Now put

\[
y_i=1-x_i\in[0,1].
\]

Then

\[
\sum_i y_i=2.
\]

Write

\[
S_k=\sum_i y_i^k.
\]

A direct substitution gives

\[
16F
=16-12S_2-8S_3+18S_4-S_2^2.
\tag{8}
\]

Cauchy gives

\[
S_2^2\le S_1S_3=2S_3,
\qquad
S_3^2\le S_2S_4.
\tag{9}
\]

Also

\[
S_2\ge\frac{S_1^2}{6}=\frac23.
\]

Thus

\[
S_3\ge\frac{S_2^2}{2}\ge\frac{S_2}{3}>rac{2S_2}{9}.
\]

For fixed `S2`, the function

\[
t\mapsto \frac{18t^2}{S_2}-8t
\]

is therefore increasing throughout the admissible range of `S3`.  Using (9),

\[
18S_4-8S_3
\ge
\frac{18S_3^2}{S_2}-8S_3
\ge
\frac92S_2^3-4S_2^2.
\]

Hence

\[
16F
\ge
H(S_2),
\qquad
H(s)=16-12s-5s^2+\frac92s^3.
\tag{10}
\]

Because `0<=y_i<=1` and `sum y_i=2`,

\[
\frac23\le s=S_2\le2.
\]

For `2/3<=s<=4/3`, `H` is decreasing, so

\[
H(s)\ge H(4/3)=\frac{16}{9}>0.
\]

For `4/3<=s<=2`, write `s=4/3+t`, `0<=t<=2/3`.  Then

\[
H(s)
=\frac92t^3+13t^2-\frac43t+\frac{16}{9}
\ge
13t^2-\frac43t+\frac{16}{9}.
\]

The last quadratic has minimum

\[
\frac{68}{39}>0.
\]

Thus `F>=0`, proving (2).

## Significance

The rank-two projection program now has two large exact pieces:

1. every two-orthogonal-direction frame is solved sharply, with maximum `q2/q1=33/160`;
2. every rank-two projection with `max_i P_ii<=1/2` satisfies the required global `1/4` contraction, regardless of its planar direction structure.

Therefore any remaining rank-two obstruction must simultaneously be genuinely multi-directional **and** contain at least one leverage score `P_ii>1/2`.  This pushes the unresolved geometry into a much thinner high-leverage boundary regime, where a coordinate-splitting or stability argument is the natural next target.
