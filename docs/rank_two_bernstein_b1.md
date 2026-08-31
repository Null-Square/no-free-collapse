# First interior Bernstein coefficient in the rank-two heavy interpolation

This note proves a global positivity statement for the first interior Bernstein coefficient in the canonical high-leverage rank-two interpolation.

Let `p,q in R^5` be orthonormal and define

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

The `x=0` coefficient `b0` is nonnegative because that endpoint is a coordinate projection plus a five-coordinate rank-one projection, hence belongs to the solved two-direction class. This note proves

\[
\boxed{b_1\ge0.}
\tag{1}
\]

The fixed-path coefficient `b4` is **not** assumed known: at `eps=1/2`, the other diagonal entries are `p_j^2/2+q_j^2`, which need not be at most one half. Numerically the sharp value of `b1` appears to be `1/32`; the proof below deliberately uses a looser positivity certificate and does not claim that sharp constant globally.

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

Since the `a_i` form a probability vector, (4) is at most its largest coordinate choice. Thus it is enough to prove the following purely scalar statement.

## 2. Quartic simplex lemma

Let `b_1,...,b_5>=0` with `sum b_i=1`. For every distinguished coordinate `x=b_i`, let the other four coordinates be `y_1,...,y_4` and write

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

For fixed `x`, `H_x` is a symmetric quartic polynomial in four nonnegative variables on a simplex. We use the established nonnegative half-degree principle for symmetric optimization (Timofte; Riener): after fixing the first power sum, it is enough to inspect every possible support size with at most two distinct positive coordinate values. Up to permutation the required patterns are

\[
1,\qquad 1+1,\qquad 1+2,\qquad 1+3,\qquad 2+2.
\]

The test suite reconstructs every certificate below exactly with Python `Fraction`; no floating-point optimization enters the proof.

### Support `1`

Put

\[
(y_1,y_2,y_3,y_4)=(1-x,0,0,0).
\]

Then

\[
H_x=-48x^4+72x^3-32x^2+x-1.
\]

The two intervals have largest degree-four Bernstein coefficients

| `x` interval | largest coefficient |
|---|---:|
| `[0,1/2]` | `-7/8` |
| `[1/2,1]` | `-5/2` |

### Pattern `1+1`

Write

\[
(y_1,y_2,y_3,y_4)=(1-x)(t,1-t,0,0).
\]

Three cells suffice:

| `x` interval | `t` interval | largest Bernstein coefficient |
|---|---|---:|
| `[0,1/2]` | `[0,1/2]` | `0` |
| `[0,1/2]` | `[1/2,1]` | `0` |
| `[1/2,1]` | `[0,1]` | `-77/36` |

### Pattern `1+2`

Write

\[
(y_1,y_2,y_3,y_4)=(1-x)\left(1-t,\frac t2,\frac t2,0\right).
\]

Again three cells suffice:

| `x` interval | `t` interval | largest Bernstein coefficient |
|---|---|---:|
| `[1/2,1]` | `[0,1]` | `-59/32` |
| `[0,1/2]` | `[0,1/2]` | `0` |
| `[0,1/2]` | `[1/2,1]` | `-15/8` |

### Pattern `2+2`

Write

\[
(y_1,y_2,y_3,y_4)=\frac{1-x}{2}(t,t,1-t,1-t).
\]

The two cells have largest degree-`(4,4)` Bernstein coefficients

| `x` interval | `t` interval | largest coefficient |
|---|---|---:|
| `[0,1/2]` | `[0,1]` | `-1/4` |
| `[1/2,1]` | `[0,1]` | `-17/12` |

### Pattern `1+3`

Write

\[
(y_1,y_2,y_3,y_4)=(1-x)\left(1-t,\frac t3,\frac t3,\frac t3\right).
\]

The three cells have largest Bernstein coefficients

| `x` interval | `t` interval | largest coefficient |
|---|---|---:|
| `[1/2,1]` | `[0,1]` | `-7/6` |
| `[0,1/2]` | `[1/2,1]` | `-7/8` |
| `[0,1/2]` | `[0,1/2]` | `0` |

All support patterns are therefore nonpositive. Equations (2)--(5) imply `b1>=0`, proving (1).

## 3. Exact high-leverage stability split

Let

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

be the squared Pluecker coordinates. Then

\[
\boxed{
8\Delta(P)
=F(d)
+8\|(I-P)DP\|_F^2
+8\sum_{i<j}m_{ij}P_{ij}^2.
}
\tag{6}
\]

Both correction terms are manifestly nonnegative. The balanced-diagonal theorem is precisely the region where `F(d)>=0`; in the remaining high-leverage region, (6) identifies the two geometric terms that must compensate a possibly negative scalar defect.

## 4. Remaining target

The coefficientwise heavy-path program currently has

- `b0>=0`: proved at `eps=0` by the two-direction theorem;
- `b1>=0`: proved here;
- `b2,b3,b4`: not globally proved for a fixed arbitrary pair `p,q`.

The earlier statement that `b4` followed from the balanced-diagonal theorem was incorrect: `P(1/2)` need not have all diagonal entries at most `1/2`. For the global rank-two theorem, a stronger direct route is now preferable: apply the exact Pluecker defect together with the sharp interval bounds on each `m_ij`, yielding a finite linear relaxation that numerically remains nonnegative even after decomposability is dropped. An analytic dual certificate for that relaxation would close rank two without proving `b2,b3,b4` separately.

### Prior-art note

The half-degree principle used in the finite quartic reduction is due to Timofte and subsequent refinements/elementary proofs by Riener; see V. Timofte, *On the positivity of symmetric polynomial functions, Part I* (2003), and C. Riener, *On the degree and half-degree principle for symmetric polynomials* (2012). The project contribution here is the reduction of the rank-two hafnian-gradient coefficient to this quartic certificate, not the general symmetric-polynomial principle itself.
