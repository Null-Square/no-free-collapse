# Global rank-three projection gradient theorem

This note closes the rank-three projection case of the six-variable hafnian-gradient problem.

For a real `6 x 6` rank-three orthogonal projection `P`, define

\[
q_1(P)=\sum_{i<j}P_{ij}^2,
\qquad
q_2(P)=\sum_{i<j}\operatorname{haf}(P_{\widehat i,\widehat j})^2.
\]

## Theorem

Every rank-three orthogonal projection satisfies

\[
\boxed{q_2(P)\le\frac14 q_1(P).}
\tag{1}
\]

The constant `1/4` is sharp.  It is attained by the equal disjoint-pair projection.  Coordinate rank-three projections give the degenerate equality `q1=q2=0`.

The proof uses the exact defect identity from [`rank_three_defect_identity.md`](rank_three_defect_identity.md), the sharp two-coordinate PSD capacities of `P` and `I-P`, and an explicit finite LP-dual certificate.  No numerical optimization is used in the proof.

## 1. Exact defect identity

Write

\[
y_i=P_{ii}-\frac12,
\qquad
x_{ij}=P_{ij}^2,
\]

and define

\[
S_2=\sum_i y_i^2,
\qquad
S_4=\sum_i y_i^4,
\]

\[
W=\sum_i\sum_{\substack{j<k\\j,k\ne i}}x_{ij}x_{ik},
\qquad
L=\sum_{i<j}x_{ij}(y_i-y_j)^2.
\]

The previous exact identity is

\[
8\left(\frac14q_1-q_2\right)
=8W+4L+10S_4-S_2-S_2^2.
\tag{2}
\]

Thus it is enough to prove that the right-hand side of (2) is nonnegative.

Projection idempotence gives

\[
r_i:=\sum_{j\ne i}x_{ij}
=\frac14-y_i^2.
\tag{3}
\]

## 2. Sharp two-coordinate capacities

Both `P` and `I-P` are PSD.  Their `2 x 2` principal minors therefore give

\[
x_{ij}\le P_{ii}P_{jj}
\]

and

\[
x_{ij}\le(1-P_{ii})(1-P_{jj}).
\]

Hence

\[
\boxed{
x_{ij}\le c_{ij}:=
\frac14+y_i y_j-\frac12|y_i+y_j|.
}
\tag{4}
\]

Because `x_ij>=0`, (4) implies `x_ij^2<=c_ij x_ij`.

Using

\[
W=\frac12\sum_i r_i^2-\sum_{i<j}x_{ij}^2,
\]

and the identities

\[
\sum_{i<j}x_{ij}=\frac34-\frac12S_2,
\]

\[
\sum_{i<j}x_{ij}y_i y_j
=\frac12\left(\frac14S_2-S_4-L\right),
\]

one obtains

\[
W\ge
-\frac14S_2+S_4+\frac12L
+\frac12\sum_{i<j}x_{ij}|y_i+y_j|.
\tag{5}
\]

Substituting (5) into (2) shows that it suffices to prove

\[
\sum_{i<j}x_{ij}
\left(4|y_i+y_j|+8(y_i-y_j)^2\right)
+18S_4
\ge
3S_2+S_2^2.
\tag{6}
\]

For fixed diagonal `y`, this is now linear in the fifteen nonnegative edge variables `x_ij`.

## 3. Explicit LP-dual certificate

Set

\[
\alpha_i=16y_i^2.
\]

For an edge define

\[
d_{ij}=4|y_i+y_j|+8(y_i-y_j)^2
\]

and

\[
\beta_{ij}=\max\{0,\alpha_i+\alpha_j-d_{ij}\}.
\]

Then, tautologically,

\[
d_{ij}\ge\alpha_i+\alpha_j-\beta_{ij}.
\]

Multiplying by `x_ij`, summing, using the row sums (3), and then the capacity (4), gives

\[
\sum_{i<j}d_{ij}x_{ij}
\ge
\sum_i\alpha_i r_i
-\sum_{i<j}\beta_{ij}c_{ij}.
\tag{7}
\]

The first term is exact:

\[
\sum_i\alpha_i r_i
=4S_2-16S_4.
\tag{8}
\]

Moreover

\[
\alpha_i+\alpha_j-d_{ij}
=4u_{ij}(2u_{ij}-1),
\qquad
u_{ij}:=|y_i+y_j|.
\tag{9}
\]

Thus `beta_ij` is nonzero only when `u_ij>1/2`.  Since `|y_i|<=1/2`, such a pair must have the same sign.  Put

\[
a=|y_i|,
\qquad b=|y_j|,
\qquad u=a+b>\frac12.
\]

Then

\[
c_{ij}=\left(\frac12-a\right)\left(\frac12-b\right)
\]

and the dual correction on that edge is

\[
g(a,b)
:=4(a+b)(2a+2b-1)
\left(\frac12-a\right)
\left(\frac12-b\right).
\tag{10}
\]

Combining (6)--(10), it remains only to show

\[
\sum_{\text{same-sign }i<j\atop |y_i|+|y_j|>1/2}g(|y_i|,|y_j|)
\le
S_2+2S_4-S_2^2.
\tag{11}
\]

## 4. Two-variable correction lemma

Define

\[
f(t)=t^2(1-4t^2),
\qquad 0\le t\le\frac12.
\]

### Lemma

For `0<=a,b<=1/2` with `a+b>=1/2`,

\[
\boxed{
g(a,b)\le\frac14\left(f(a)+f(b)\right).}
\tag{12}
\]

### Proof

Set

\[
u=a+b,
\qquad v=a-b,
\qquad w=v^2.
\]

The domain is

\[
\frac12\le u\le1,
\qquad
0\le w\le(1-u)^2.
\]

A direct expansion gives

\[
\frac14(f(a)+f(b))-g(a,b)
=-\frac18N_u(w),
\]

where

\[
N_u(w)
=w^2+(-10u^2+8u-1)w
+u(u-1)(17u^2-23u+8).
\tag{13}
\]

For fixed `u`, `N_u(w)` is convex in `w`, so its maximum on the interval occurs at an endpoint.  At the endpoints,

\[
N_u(0)
=u(u-1)(17u^2-23u+8)\le0,
\]

because `17u^2-23u+8>0` (its discriminant is `-15`), and

\[
N_u((1-u)^2)
=2u(u-1)(2u-1)^2\le0.
\]

Therefore `N_u(w)<=0` throughout the interval, proving (12).

## 5. Reduction to three sign-count cases

Let the positive imbalances have magnitudes `a_1,...,a_m` and the negative imbalances have magnitudes `b_1,...,b_n`.  Zeros do not contribute to (11).  Since `sum_i y_i=0`, unless every `y_i=0` both signs occur.

If zeros are present, pad one sign class with zeros.  Since `f(0)=0`, this can only weaken the bound.  After possibly flipping every sign, it is enough to consider

\[
m+n=6,
\qquad
m\in\{1,2,3\},
\qquad n=6-m.
\]

Write

\[
A_2=\sum_{i=1}^m a_i^2,
\quad
A_4=\sum_{i=1}^m a_i^4,
\]

and similarly `B2,B4` for the negative class.

By (12), each variable appears in at most `m-1` correction pairs inside the positive class and at most `n-1` inside the negative class.  Therefore

\[
\text{correction}
\le
\frac{m-1}{4}(A_2-4A_4)
+
\frac{n-1}{4}(B_2-4B_4).
\tag{14}
\]

Subtracting (14) from the right-hand side of (11) gives

\[
D
=
\frac{n-1}{4}A_2
+
\frac{m-1}{4}B_2
+(m+1)A_4
+(n+1)B_4
-(A_2+B_2)^2.
\tag{15}
\]

Cauchy gives

\[
A_4\ge\frac{A_2^2}{m},
\qquad
B_4\ge\frac{B_2^2}{n}.
\]

Also

\[
0\le A_2\le\frac m4,
\qquad
0\le B_2\le\frac n4.
\]

Hence

\[
D\ge E_{m,n}(A_2,B_2),
\]

where

\[
E_{m,n}
=
\frac{n-1}{4}A_2
+
\frac{m-1}{4}B_2
+\frac{A_2^2}{m}
+\frac{B_2^2}{n}
-2A_2B_2.
\tag{16}
\]

Only three cases remain.

### Case `m=1,n=5`

\[
E_{1,5}
=A_2-4A_2^2
+\frac15(B_2-5A_2)^2
\ge0,
\]

because `A2<=1/4`.

### Case `m=2,n=4`

Here

\[
E_{2,4}
=\frac{2A_2^2-8A_2B_2+3A_2+B_2^2+B_2}{4}.
\]

For fixed `A2`, this is a convex quadratic in `B2` on `[0,1]`.

- If `A2<=1/8`, the minimum occurs at `B2=0`, giving
  \[
  E_{2,4}=\frac14A_2(2A_2+3)\ge0.
  \]
- If `1/8<=A2<=3/8`, the minimum occurs at `B2=4A2-1/2`, giving
  \[
  E_{2,4}
  =-\frac{56A_2^2-28A_2+1}{16}\ge0.
  \]
  The quadratic in the numerator is convex and is negative at both endpoints of this interval.
- If `3/8<=A2<=1/2`, the minimum occurs at `B2=1`, giving
  \[
  E_{2,4}=\frac{(A_2-2)(2A_2-1)}4\ge0.
  \]

### Case `m=n=3`

Put `R=A2+B2`.  Then

\[
E_{3,3}
=\frac12R+\frac13R^2-\frac83A_2B_2.
\]

Since `A2 B2<=R^2/4`,

\[
E_{3,3}
\ge
\frac12R-\frac13R^2
=\frac{R(3-2R)}6\ge0,
\]

because `R<=3/2`.

Thus `D>=0` in every possible sign-count case.  This proves (11), then (6), then (2), and therefore the rank-three contraction (1).

## Significance

The projection-gradient program is now completely solved in ranks one, three, and five; ranks zero and six are trivial, and rank four is the complement of rank two.  The only unresolved nontrivial projection rank is therefore **rank two** (equivalently rank four).

The proof also clarifies why the disjoint-pair geometry is sharp: equal diagonals remove the scalar imbalance term, while the PSD capacity argument prevents imbalanced rank-three projections from improving the `1/4` constant.
