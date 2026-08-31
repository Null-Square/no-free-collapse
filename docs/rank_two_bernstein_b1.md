# First interior Bernstein coefficient in the rank-two heavy interpolation

This note proves one of the three remaining interior Bernstein inequalities in the final high-leverage rank-two projection problem.

Let `p,q in R^5` be orthonormal and define the heavy-coordinate path

\[
P(\varepsilon)=
\begin{pmatrix}
1-\varepsilon & \sqrt{\varepsilon(1-\varepsilon)}p^T\\
\sqrt{\varepsilon(1-\varepsilon)}p &
\varepsilon pp^T+qq^T
\end{pmatrix},
\qquad 0\le\varepsilon\le\frac12.
\]

Put `x=2 eps` and

\[
\Delta(x)=\frac14 q_1(P(x/2))-q_2(P(x/2)).
\]

Write the degree-four Bernstein expansion

\[
\Delta(x)=\sum_{k=0}^4 b_k\binom4k x^k(1-x)^{4-k}.
\]

The endpoint coefficients `b0,b4` are nonnegative by the already-proved endpoint theorems.  This note proves

\[
\boxed{b_1\ge0.}
\tag{1}
\]

Numerically the sharp value appears to be `1/32`; the proof below deliberately uses a looser positivity certificate and does **not** claim the sharp constant globally.

## 1. Exact five-frame moments

Set

\[
a_i=p_i^2,\qquad b_i=q_i^2,
\]

and abbreviate

\[
B_k=\sum_i b_i^k,\qquad
AB^j=\sum_i a_i b_i^j,
\qquad
T_b=\sum_i p_i q_i^3.
\]

A direct expansion of the exact rank-two defect gives

\[
32 b_1
=-\Bigl[
6(AB)B_2-13AB+36AB^2-36AB^3
+36B_2^2-70B_2+96B_3-72B_4
+12T_b^2+9
\Bigr].
\tag{2}
\]

Because `p` is unit and orthogonal to `q`,

\[
T_b
=\langle p,q^{\circ3}-B_2q\rangle,
\]

so

\[
T_b^2\le \|q^{\circ3}-B_2q\|^2
=B_3-B_2^2.
\tag{3}
\]

Substituting (3) into the bracket in (2) gives the upper bound

\[
\sum_i a_i h(b_i)
+24B_2^2-70B_2+108B_3-72B_4+9,
\tag{4}
\]

where

\[
h(t)=(6B_2-13)t+36t^2-36t^3.
\]

Since the `a_i` form a probability vector, (4) is at most its largest coordinate choice.  Thus it is enough to prove the following purely scalar statement.

## 2. Quartic simplex lemma

Let `b_1,...,b_5>=0` with `sum b_i=1`.  For every distinguished coordinate `x=b_i`, let the other four coordinates be `y_1,...,y_4` and write

\[
R_k=\sum_{j=1}^4y_j^k,\qquad \sum_jy_j=1-x.
\]

Then

\[
H_x(y)
=24R_2^2+(48x^2+6x-70)R_2+108R_3-72R_4
-48x^4+78x^3-34x^2-13x+9
\le0.
\tag{5}
\]

For fixed `x`, `H_x` is a symmetric quartic polynomial in four nonnegative variables on a simplex.  We use the established half-degree principle for symmetric optimization (Timofte; Riener): a quartic symmetric polynomial on such a fixed-first-power-sum slice need only be checked on points with at most two distinct components.  This is prior art and is not a novelty claim of the project.

Therefore only two multiplicity patterns remain.

### Pattern `2+2`

Write

\[
(y_1,y_2,y_3,y_4)
=\frac{1-x}{2}(t,t,1-t,1-t),
\qquad 0\le x,t\le1.
\]

After substitution, (5) becomes the bivariate quartic `P22(x,t)<=0`.  An exact rational Bernstein certificate needs only the two cells

| `x` interval | `t` interval | largest degree-`(4,4)` Bernstein coefficient |
|---|---|---:|
| `[0,1/2]` | `[0,1]` | `-1/4` |
| `[1/2,1]` | `[0,1]` | `-17/12` |

Hence `P22<=0` everywhere.

### Pattern `1+3`

Write

\[
(y_1,y_2,y_3,y_4)
=(1-x)\left(1-t,\frac t3,\frac t3,\frac t3\right),
\qquad 0\le x,t\le1.
\]

The exact rational Bernstein certificate uses three cells:

| `x` interval | `t` interval | largest Bernstein coefficient |
|---|---|---:|
| `[1/2,1]` | `[0,1]` | `-7/6` |
| `[0,1/2]` | `[1/2,1]` | `-7/8` |
| `[0,1/2]` | `[0,1/2]` | `0` |

Thus `P13<=0` everywhere as well.  The test suite reconstructs these coefficients exactly with Python `Fraction`, so no floating-point optimization enters the certificate.

Equations (2)--(5) therefore imply `b1>=0`, proving (1).

## 3. Exact high-leverage stability split

The same rank-two algebra gives a second useful identity.  Let

\[
d_i=P_{ii},\qquad D_k=\sum_i d_i^k,
\]

and define

\[
F(d)=-6+19D_2-32D_3+18D_4-D_2^2.
\]

Let `D=diag(d)` and let

\[
m_{ij}=d_i d_j-P_{ij}^2
\]

be the squared Pluecker coordinates.  Then

\[
\boxed{
8\Delta(P)
=F(d)
+8\|(I-P)DP\|_F^2
+8\sum_{i<j}m_{ij}P_{ij}^2.
}
\tag{6}
\]

Both correction terms are manifestly nonnegative.  To derive (6), write

\[
C=\sum_{i<j}d_i d_jm_{ij},\qquad
J=\sum_{i<j}m_{ij}P_{ij}^2.
\]

Since `d_i d_j=m_ij+P_ij^2`,

\[
24C-8\sum m_{ij}^2=16C+8J.
\]

If `U` is a Parseval frame with `P=UU^T`, then

\[
C=\det(U^TDU)
\]

and

\[
C-\frac{D_2^2-D_3}{2}
=\frac12\|(I-P)DP\|_F^2.
\]

Substitution into the exact Pluecker defect yields (6).

The balanced-diagonal theorem is precisely the region where the scalar term `F(d)` is already nonnegative.  In the remaining high-leverage region, (6) identifies the only two geometric resources that must compensate a potentially negative `F`: failure of the leverage diagonal to commute with the projection, and genuine angle mixing between Pluecker and Gram weights.

## 4. Remaining target

The heavy-path Bernstein program is now

- `b0>=0`: proved at the coordinate-plus-rank-one endpoint;
- `b1>=0`: proved here;
- `b2,b3`: open;
- `b4>=0`: proved by the balanced-diagonal theorem.

So only the two middle coefficients remain.  Their numerical minima are still

\[
b_2\approx\frac{25}{768},\qquad
b_3\approx\frac7{256},
\]

attained by the same two-direction degeneration, but those constants remain conjectural until an analytic or exact computer certificate is obtained.

### Prior-art note

The half-degree principle used in the finite quartic reduction is due to Timofte and subsequent refinements/elementary proofs by Riener; see V. Timofte, *On the positivity of symmetric polynomial functions, Part I* (2003), and C. Riener, *On the degree and half-degree principle for symmetric polynomials* (2012).  The project contribution here is the reduction of the rank-two hafnian-gradient coefficient to this quartic certificate, not the general symmetric-polynomial principle itself.
