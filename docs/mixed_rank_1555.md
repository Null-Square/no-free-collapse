# Mixed nested coefficient `(1,5,5,5)`

This note proves the first genuinely mixed coefficient in the nested spectral-kernel program.

Let

\[
P=uu^T,
\qquad
Q=I-vv^T,
\qquad
u,v\in\mathbb R^6,
\qquad
\|u\|=\|v\|=1,
\qquad
u\perp v.
\]

Then `P` has rank one, `Q` has rank five, and `P<=Q`.  For the homogenized spectral kernel `D` from [`spectral_polarization.md`](spectral_polarization.md), we prove

\[
\boxed{\mathcal D(P,Q,Q,Q)\ge0.}
\tag{1}
\]

By simultaneous complementation and symmetry of the four slots, the rank pattern `(1,1,1,5)` follows as well.

The proof is self-contained apart from the exact definition of `D`.

## 1. Exact quadratic-form reduction

Set

\[
t_i=v_i^2,
\qquad
\mu=\sum_i t_i^2.
\]

Thus `t_i>=0`, `sum t_i=1`.  Define

\[
\phi_\mu(t)
=1-\mu+(6\mu-5)t+24t^2-36t^3.
\tag{2}
\]

A direct expansion of the kernel gives

\[
\boxed{
16\mathcal D(P,Q,Q,Q)
=
\sum_i\phi_\mu(t_i)u_i^2
+12\left(\sum_i u_i v_i^3\right)^2.
}
\tag{3}
\]

The implementation `rank_one_rank_five_mixed_formula` is tested against the original kernel definition.

For reference, (3) can also be obtained from the square-free product.  Since `off(Q)=-off(vv^T)`, the pair-inner-product part contributes a rank-one quadratic term plus the constant rank-five edge energy, and the three quartic pairings are identical.  Collecting terms gives (2)--(3).

## 2. Light coordinates have positive diagonal coefficient

The first lemma is the main structural simplification.

### Lemma 1

Assume one coordinate has mass

\[
s=\max_i t_i>\frac12.
\]

For every other coordinate `x=t_j`,

\[
\boxed{\phi_\mu(x)\ge2x.}
\tag{4}
\]

If no coordinate exceeds `1/2`, then every `phi_mu(t_i)` is already nonnegative, so (3) is immediate.  Thus (4) reduces the only difficult case to one heavy atom.

### Proof for `x<=1/6`

Write `r=1-s`, so `x<=r<1/2`.  The identity

\[
\phi_\mu(x)-2x
=(1-\mu)(1-6x)+x(-1+24x-36x^2)
\tag{5}
\]

is useful.  Since

\[
1-\mu=2e_2(t)\ge2s(1-s)=2sr\ge x,
\]

and `1-6x>=0`, (5) gives

\[
\phi_\mu(x)-2x
\ge
x(1-6x)+x(-1+24x-36x^2)
=18x^2(1-2x)\ge0.
\]

### Proof for `1/6<=x<=r`

The four remaining coordinates have total mass `r-x`.  Cauchy gives

\[
\mu
\ge
s^2+x^2+\frac{(r-x)^2}{4}.
\]

Since `1-6x<=0`, (5) is bounded below by

\[
G(r,x)
=
\left(1-s^2-x^2-\frac{(r-x)^2}{4}\right)(1-6x)
+x(-1+24x-36x^2),
\tag{6}
\]

where `s=1-r`.

Parameterize the triangle

\[
x=\frac16+\frac p3,
\qquad
r=x+\left(\frac12-x\right)q,
\qquad
0\le p,q\le1.
\]

In the tensor Bernstein basis of degrees `(3,2)` on the unit square, the exact coefficient matrix of `G` is

\[
\begin{pmatrix}
1/3 & 1/3 & 1/3\\
16/27 & 11/27 & 17/54\\
7/9 & 2/3 & 5/9\\
0 & 0 & 0
\end{pmatrix}.
\tag{7}
\]

Every Bernstein basis function is nonnegative on `[0,1]^2`, so `G>=0`.  This proves Lemma 1.

## 3. Eliminate the unique heavy coordinate

Assume now

\[
s=t_1>\frac12,
\qquad
x_i=t_i\quad(i=2,\dots,6),
\qquad
\nu=\sum_{i=2}^6x_i^2.
\]

If `phi_mu(s)>=0`, equation (3) is already nonnegative by Lemma 1.  It remains to treat

\[
\phi_\mu(s)<0.
\]

Write `y=(u_2,...,u_6)`.  Let

\[
a_i=\sqrt{x_i},
\qquad
c_i=\sqrt{x_i}(x_i-s),
\qquad
D=\operatorname{diag}(\phi_\mu(x_i)).
\]

Orthogonality `u dot v=0` gives

\[
u_1=-\frac{a^Ty}{\sqrt s},
\]

while

\[
\sum_i u_i v_i^3=c^Ty.
\]

Therefore

\[
16\mathcal D
=
y^T(D+12cc^T)y
+\frac{\phi_\mu(s)}s(a^Ty)^2.
\tag{8}
\]

By Lemma 1, `D` is positive on every active light coordinate.

Set

\[
G=D+12cc^T.
\]

For `S=a^TG^{-1}a`, Cauchy in the `G` metric gives

\[
(a^Ty)^2\le S\,y^TGy.
\]

Because `phi_mu(s)<0`, (8) yields

\[
16\mathcal D
\ge
\left(1+\frac{\phi_\mu(s)}sS\right)y^TGy.
\tag{9}
\]

It remains to upper-bound `S`.

## 4. A one-parameter variational bound

For every real `lambda`, the elementary inequalities

\[
2(a-\lambda c)^Ty-y^TDy
\le
\sum_i\frac{(a_i-\lambda c_i)^2}{\phi_\mu(x_i)}
\]

and

\[
2\lambda(c^Ty)-12(c^Ty)^2\le\frac{\lambda^2}{12}
\]

imply

\[
S
\le
\sum_i\frac{(a_i-\lambda c_i)^2}{\phi_\mu(x_i)}
+\frac{\lambda^2}{12}.
\tag{10}
\]

Choose the specific value

\[
\lambda=-\frac1s.
\]

Then

\[
a_i-\lambda c_i
=\frac{x_i^{3/2}}s,
\]

so Lemma 1 gives

\[
S
\le
\frac1{s^2}
\left(
\sum_i\frac{x_i^3}{\phi_\mu(x_i)}+\frac1{12}
\right)
\le
\frac1{s^2}
\left(
\frac\nu2+\frac1{12}
\right).
\tag{11}
\]

Thus (9) is nonnegative if

\[
K(s,\nu)
:=
s^3+\phi_\mu(s)\left(\frac1{12}+\frac\nu2\right)
\ge0.
\tag{12}
\]

In fact we prove the stronger bound `K>=1/12`.

## 5. Exact scalar certificate

Put

\[
r=1-s\in[0,1/2],
\qquad
z=\frac\nu{r^2}\in[1/5,1],
\]

where the interval for `z` follows from Cauchy and concentration on the five light coordinates.  For `r=0`, direct substitution gives `K=1/12`.

For `r>0`, one finds

\[
K-\frac1{12}
=-\frac r{12}B(r,z),
\tag{13}
\]

where

\[
\begin{aligned}
B(r,z)={}&
36r^4z^2-180r^4z-30r^3z^2+402r^3z\\
&-288r^2z-18r^2+61rz+31r-13.
\end{aligned}
\tag{14}
\]

Map the rectangle to the unit square by

\[
r=\frac X2,
\qquad
z=\frac{1+4Y}{5},
\qquad
0\le X,Y\le1.
\]

The tensor Bernstein coefficient matrix of `B` of degrees `(4,2)` is

\[
\begin{pmatrix}
-13 & -13 & -13\\
-38/5 & -91/20 & -3/2\\
-107/20 & -81/20 & -11/4\\
-151/40 & -83/20 & -41/8\\
-64/25 & -19/5 & -6
\end{pmatrix}.
\tag{15}
\]

Every coefficient is strictly negative, hence `B<=0` on the full rectangle.  Therefore

\[
\boxed{K(s,\nu)\ge\frac1{12}>0.}
\tag{16}
\]

Combining (9), (11), and (16) proves (1).

## 6. Complementary rank pattern

Simultaneous complementation changes every off-diagonal projection entry by a sign.  The pair-inner-product and quartic product terms in `D` are unchanged.  Hence

\[
\mathcal D(P,Q,Q,Q)
=
\mathcal D(I-Q,I-Q,I-Q,I-P).
\]

The theorem therefore proves both nested patterns

\[
\boxed{(1,5,5,5)\quad\text{and}\quad(1,1,1,5).}
\]

These are the first non-diagonal mixed spectral coefficients proved analytically after the complete projection theorem.
