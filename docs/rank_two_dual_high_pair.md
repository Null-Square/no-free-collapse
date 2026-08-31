# Rank-two quadratic dual: global high-pair theorem

This note proves the rank-two gradient contraction on the entire diagonal regime in which the two largest leverage scores have sum at least `3/2`.  The proof is stronger than needed in one respect: after the exact rank-two Pluecker reduction, it discards decomposability and proves nonnegativity of an explicit fractional-matching dual certificate depending only on the projection diagonal.

Let `P` be a real `6 x 6` rank-two orthogonal projection and write

\[
\Delta(P)=\frac14q_1(P)-q_2(P),
\qquad d_i=P_{ii}.
\]

As always,

\[
0\le d_i\le1,
\qquad \sum_i d_i=2.
\]

The theorem proved here is:

> **High-pair theorem.** If the two largest diagonal entries satisfy
> \[
> d_{(1)}+d_{(2)}\ge\frac32,
> \]
> then
> \[
> \boxed{\Delta(P)\ge0.}
> \]
> In fact the explicit diagonal-only certificate below is strictly positive unless the diagonal is a permutation of `(1,1,0,0,0,0)`.

## 1. Chord relaxation and quadratic dual

For a rank-two projection, let

\[
m_{ij}=d_i d_j-P_{ij}^2
\]

be the squared Pluecker coordinates.  They satisfy

\[
\ell_{ij}:=(d_i+d_j-1)_+
\le m_{ij}\le
u_{ij}:=d_i d_j,
\]

and the exact marginals

\[
\sum_{j\ne i}m_{ij}=d_i.
\]

From the exact rank-two defect identity,

\[
8\Delta
=-6+19D_2-24D_3+18D_4-9D_2^2
+\sum_{i<j}8m_{ij}(3\nu_{ij}-m_{ij}),
\]

where `D_k=sum_i d_i^k`.

For fixed `d`, the edge term is concave in `m`.  Its chord on
`[ell,nu]` is

\[
8m(3\nu-m)
\ge
8\bigl[\nu\ell+(2\nu-\ell)m\bigr].
\tag{1}
\]

Now choose the node potential

\[
\boxed{\alpha(t)=-t^2+\frac32t-\frac14.}
\tag{2}
\]

For an edge with `s=x+y`, the coefficient residual after summing node
potentials is

\[
(2xy-\ell)-\alpha(x)-\alpha(y)
=
\begin{cases}
\frac12(s-1)(2s-1), & s\le1,\\[1mm]
\frac12(s-1)(2s-3), & s\ge1.
\end{cases}
\tag{3}
\]

Use the lower endpoint `m=ell` when the residual is nonnegative and the
upper endpoint `m=nu` when it is negative.  After summing all edges and using
`sum_j m_ij=d_i`, the resulting diagonal-only lower bound can be written very
compactly.

Set

\[
r_i=d_i(1-d_i),
\qquad
R=\sum_i r_i,
\qquad
E=\sum_i r_i^2.
\]

Then

\[
\boxed{8\Delta(P)\ge H(d),}
\tag{4}
\]

where

\[
\begin{aligned}
H(d)
={}&R(1-R)+2E\\
&-\sum_{d_i+d_j<1/2}
4d_i d_j(1-d_i-d_j)(1-2d_i-2d_j)\\
&-\sum_{d_i+d_j>3/2}
4(1-d_i)(1-d_j)(d_i+d_j-1)(2d_i+2d_j-3).
\end{aligned}
\tag{5}
\]

Equation (5) is the quadratic-dual certificate implemented in
`rank_two_quadratic_dual_certificate`.

## 2. Geometry of the high-pair regime

Sort

\[
d_1\ge d_2\ge\cdots\ge d_6
\]

and assume

\[
d_1+d_2\ge\frac32.
\]

Put

\[
a=1-d_1,
\qquad b=1-d_2,
\qquad z=(d_3,d_4,d_5,d_6).
\]

Since `sum d_i=2`,

\[
T:=a+b=\sum_{j=3}^6d_j\le\frac12.
\tag{6}
\]

There is only one possible `>3/2` edge, namely `(1,2)`.  Indeed if, say,
`d_1+d_3>3/2`, then `d_2>=d_3` would force the first three coordinates to
have sum strictly larger than two.  Likewise `d_2>=1/2`, so no edge joining
`d_1` or `d_2` to a small coordinate can have sum below `1/2`.

Conversely every pair among `d_3,...,d_6` has sum at most `T<=1/2`.
Therefore the two correction sums in (5) consist exactly of

- the complement pair `(a,b)`; and
- all six pairs inside `z`.

For `T=0` the diagonal is `(1,1,0,0,0,0)` and `H=0`.  Assume `T>0` and define
probability vectors

\[
u=(a,b)/T\in\Delta_2,
\qquad
v=z/T\in\Delta_4.
\]

Write

\[
A_k=\sum_{i=1}^2u_i^k,
\qquad
Z_k=\sum_{j=1}^4v_j^k.
\]

A direct substitution into (5) gives

\[
\begin{aligned}
H={}&2T+[3(A_2+Z_2)-8]T^2\\
&+16[(A_2-A_3)+(Z_2-Z_3)]T^3\\
&+Q\,T^4,
\end{aligned}
\tag{7}
\]

with

\[
Q=-9A_2^2-2A_2Z_2-8A_3+18A_4-9Z_2^2-8Z_3+18Z_4.
\tag{8}
\]

## 3. Bernstein positivity on `0<=T<=1/2`

Set

\[
x=2T\in[0,1]
\]

and write (7) in the degree-four Bernstein basis

\[
H=\sum_{k=0}^4\beta_k\binom4k x^k(1-x)^{4-k}.
\tag{9}
\]

The coefficients simplify to

\[
\beta_0=0,
\qquad
\beta_1=\frac14,
\tag{10}
\]

\[
\beta_2
=\frac16+\frac{A_2+Z_2}{8},
\tag{11}
\]

\[
\beta_3
=-\frac14+\frac78(A_2+Z_2)-\frac12(A_3+Z_3),
\tag{12}
\]

and

\[
\beta_4
=-1+\frac{11}{4}S_2-\frac52S_3+\frac98S_4
-\frac9{16}S_2^2+A_2Z_2,
\tag{13}
\]

where `S_k=A_k+Z_k`.

Because `u` is a probability vector on two atoms and `v` on four atoms,

\[
A_2\ge\frac12,
\qquad
Z_2\ge\frac14.
\]

Hence

\[
\boxed{\beta_2\ge\frac{25}{96}.}
\tag{14}
\]

Also `A_3<=A_2` and `Z_3<=Z_2`, so

\[
\beta_3
\ge -\frac14+\frac38(A_2+Z_2)
\ge\boxed{\frac1{32}}.
\tag{15}
\]

It remains only to control `beta4`.

## 4. The last endpoint coefficient

For a two-atom probability vector,

\[
A_3=\frac{3A_2-1}{2},
\qquad
A_4=\frac{A_2^2+2A_2-1}{2}.
\tag{16}
\]

Substitute (16) into (13).  The derivative with respect to `A2` is

\[
\frac{\partial\beta_4}{\partial A_2}
=\frac{1-Z_2}{8}\ge0.
\]

Thus `beta4` is minimized at `A2=1/2`, i.e. at the uniform two-atom vector.
After this substitution,

\[
\beta_4
\ge
\frac{-9Z_2^2+43Z_2-40Z_3+18Z_4-4}{16}.
\tag{17}
\]

Let `e_k` be the elementary symmetric polynomials of the four probabilities
`v_1,...,v_4`.  Newton's identities turn the numerator in (17) into

\[
2(4-e_2-24e_3-36e_4).
\]

Maclaurin gives

\[
e_2\le\binom42\left(\frac14\right)^2=\frac38,
\]

\[
e_3\le\binom43\left(\frac14\right)^3=\frac1{16},
\]

\[
e_4\le\left(\frac14\right)^4=\frac1{256}.
\]

Therefore

\[
\boxed{
\beta_4\ge
\frac{4-3/8-24/16-36/256}{8}
=\frac{127}{512}>0.
}
\tag{18}
\]

Equality in this endpoint bound occurs when both normalized groups are
uniform.

## 5. Conclusion

Every Bernstein basis function in (9) is nonnegative on `[0,1]`.  Equations
(10), (14), (15), and (18) therefore prove

\[
H(d)\ge0
\]

throughout the high-pair region.  If `T>0`, at least the `beta1` term is
strictly positive, so `H>0`.  The only equality diagonal is therefore

\[
(1,1,0,0,0,0)
\]

up to permutation.

Combining this with (4) proves the rank-two projection-gradient contraction
on the entire regime

\[
\boxed{d_{(1)}+d_{(2)}\ge\frac32.}
\]

The remaining diagonal-only target is the complementary region

\[
d_{(1)}+d_{(2)}\le\frac32.
\]

The same quadratic dual remains numerically positive there, with apparent
minimum `4/27` at the uniform diagonal `(1/3,...,1/3)`, but that low-pair
statement is not claimed in this note.
