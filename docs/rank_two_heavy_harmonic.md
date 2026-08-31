# Rank-two harmonic and heavy-coordinate reductions

This note records two exact reformulations of the remaining rank-two projection-gradient problem.  They are reductions, not a proof of the last high-leverage case.

Let `P` be a real `6 x 6` rank-two orthogonal projection and write

\[
q_1(P)=\sum_{i<j}P_{ij}^2,
\qquad
q_2(P)=\sum_{i<j}\operatorname{haf}(P_{\widehat i,\widehat j})^2,
\qquad
\Delta(P)=\frac14q_1(P)-q_2(P).
\]

The balanced-diagonal theorem already proves `Delta>=0` when every `P_ii<=1/2`.  Thus the unresolved region has at least one leverage score strictly above one half.

## 1. Planar harmonic form

Choose any Parseval frame `P=UU^T` with row vectors

\[
u_i=\sqrt{d_i}(\cos\theta_i,\sin\theta_i),
\qquad d_i=P_{ii}.
\]

Set

\[
z_i=e^{2i\theta_i}.
\]

The two Parseval equations for the latent coordinates are exactly

\[
\boxed{\sum_i d_i z_i=0.}
\tag{1}
\]

Define

\[
D_k=\sum_i d_i^k,
\qquad
A=\sum_i d_i^2 z_i,
\qquad
B=\sum_i d_i^2 z_i^2,
\qquad
S=\frac12\sum_i d_i z_i^2.
\]

A latent-plane rotation multiplies `A,S` by a unit phase and `B` by its square, so `|A|,|B|,|S|` are intrinsic.

If

\[
m_{ij}=d_i d_j-P_{ij}^2
=d_i d_j\sin^2(\theta_i-\theta_j),
\]

then the two Pluecker energies in the exact defect have the Fourier forms

\[
C:=\sum_{i<j}d_i d_jm_{ij}
=\frac{D_2^2-|A|^2}{4},
\tag{2}
\]

and

\[
J:=\sum_{i<j}m_{ij}P_{ij}^2
=\frac{D_2^2-|B|^2}{16}.
\tag{3}
\]

Since

\[
24C-8\sum m_{ij}^2=16C+8J,
\]

substitution into the exact Pluecker defect gives

\[
\boxed{
8\Delta
=-6+19D_2-24D_3+18D_4
-\frac92D_2^2
-4|A|^2-\frac12|B|^2.
}
\tag{4}
\]

Equation (4) is exact and has been checked in tests against direct complementary-hafnian evaluation.

### Bessel coupling

Put probability mass `d_i/2` on `z_i` and let the random variable `D` take value `d_i` at the same atom.  Equation (1) is `E Z=0`.  The functions

\[
1,\qquad \overline Z,\qquad
\frac{\overline Z^2-\overline S}{\sqrt{1-|S|^2}}
\]

are orthonormal whenever `|S|<1`.  Bessel therefore yields the exact useful constraint

\[
\boxed{
|A|^2+
\frac{|B-D_2S|^2}{1-|S|^2}
\le 2D_3-D_2^2.
}
\tag{5}
\]

The endpoint `|S|=1` is geometrically meaningful: equality in the triangle inequality forces all `z_i^2` with positive weight to have the same phase, so the frame uses two orthogonal latent directions.  That entire class is already solved sharply, with maximum `q2/q1=33/160`.

Thus `1-|S|^2` is a quantitative angular-mixing parameter connecting the known two-direction theorem to the genuinely multidirectional problem.

## 2. Canonical heavy-coordinate interpolation

Assume coordinate zero has leverage

\[
t=P_{00}>\frac12,
\qquad
\varepsilon=1-t\in(0,1/2).
\]

Rotate the latent plane so the corresponding row lies on the first latent axis.  Projection idempotence then gives the exact representation

\[
\boxed{
P(\varepsilon)=
\begin{pmatrix}
1-\varepsilon & \sqrt{\varepsilon(1-\varepsilon)}\,p^T\\
\sqrt{\varepsilon(1-\varepsilon)}\,p &
\varepsilon pp^T+qq^T
\end{pmatrix},
}
\tag{6}
\]

where `p,q in R^5` are orthonormal.

The endpoints are both solved:

- `eps=0` is a coordinate projection plus a five-coordinate rank-one projection;
- `eps=1/2` lies on the balanced-diagonal boundary already covered by the rank-two balanced theorem.

The squared Pluecker coordinates simplify especially strongly.  If

\[
n_{jk}=\det((p_j,q_j),(p_k,q_k))^2,
\]

then

\[
\boxed{
m_{0j}=(1-\varepsilon)q_j^2,
\qquad
m_{jk}=\varepsilon n_{jk}.
}
\tag{7}
\]

So the high-leverage interpolation is driven by a five-coordinate rank-two Pluecker distribution rather than a new six-variable object.

## 3. Bernstein diagnostic

For fixed orthonormal `p,q`, `Delta(P(eps))` is a degree-four polynomial in `eps`.  Put

\[
x=2\varepsilon\in[0,1]
\]

and write it in the degree-four Bernstein basis

\[
\Delta(P(x/2))
=\sum_{k=0}^4 b_k\binom4k x^k(1-x)^{4-k}.
\]

The endpoint coefficients `b0` and `b4` are nonnegative by the two already-proved endpoint theorems.

Extensive random scans and constrained Stiefel optimization have so far found

\[
b_1\ge\frac1{32},
\qquad
b_2\ge\frac{25}{768},
\qquad
b_3\ge\frac7{256},
\tag{8}
\]

with the apparent minima attained by two-direction degenerations.  **Equation (8) is conjectural and is not used as a theorem.**  If these three coefficient bounds are proved, the remaining high-leverage rank-two region closes immediately because every Bernstein basis function is nonnegative on `[0,1]`.

The next proof target is therefore finite and concrete: establish positivity of the three interior Bernstein coefficients, preferably after rewriting them in the five-coordinate Pluecker variables from (7), or replace (8) with an equally strong harmonic stability estimate using (5).
