# Exact diagonal-completion reduction for the six-variable PSD problem

The six-variable hafnian problem can be reduced exactly to the off-diagonal
quadratic coefficients.  This removes the diagonal entries of the PSD matrix as
optimization variables and exposes the precise geometric obstruction that
separates the PSD problem from the bounded degree-two relaxation.

Let `C` be a real symmetric `6 x 6` matrix with zero diagonal and define

\[
r_C(x)=\sum_{i<j}C_{ij}x_i x_j,
\qquad
s(C)=\max_{x\in\{\pm1\}^6} r_C(x).
\]

Because `r_C` has zero Boolean mean, `s(C)>=0`.

## Minimum-trace PSD completion

Define

\[
\boxed{
\tau(C)=\min\left\{\sum_i d_i:
\operatorname{diag}(d)+\frac12 C\succeq0\right\}.
}
\]

If

\[
B=\operatorname{diag}(d)+\frac12 C,
\]

then `2 offdiag(B)=C`, and for every Boolean vector

\[
x^TBx=\sum_i d_i+r_C(x).
\]

Therefore, writing `t=Tr(B)`,

\[
\max_x x^TBx=t+s(C).
\]

For fixed `C`, the hafnian objective and `s(C)` do not depend on the diagonal.
Moreover every PSD realization has

\[
t\ge\tau(C).
\]

Since

\[
t\mapsto t\bigl(t+s(C)\bigr)
\]

is increasing for `t>=0`, the hardest PSD realization of a fixed off-diagonal
matrix is always a minimum-trace completion.

## Exact equivalent conjecture

The proposed six-variable inequality

\[
54|\operatorname{haf}(C)|\le t\,s(C)\,[t+s(C)]
\]

for every PSD realization is therefore **equivalent** to the purely
off-diagonal statement

\[
\boxed{
54|\operatorname{haf}(C)|
\le
\tau(C)\,s(C)\,[\tau(C)+s(C)]
}
\tag{1}
\]

for every real symmetric zero-diagonal `C`.

The forward implication is obtained by applying the PSD statement to a
minimum-trace completion.  The reverse implication follows from `t>=tau(C)`
and monotonicity in `t`.

Thus the global six-variable theorem is a 15-variable inequality with one
finite cube maximum and one semidefinite support function.  No diagonal search
is needed.

If (1) holds, then

\[
\tau s\le\frac{(\tau+s)^2}{4}
\]

gives

\[
|\operatorname{haf}(C)|
\le\frac{[\tau(C)+s(C)]^3}{216}.
\]

After normalizing the Boolean cube maximum to one, this is exactly the desired
`1/216` bound.

## Elliptope dual

The completion problem has a particularly useful SDP dual.  Put

\[
A=\frac12 C.
\]

The primal is

\[
\min_d\;\mathbf 1^Td
\quad\text{subject to}\quad
\operatorname{diag}(d)+A\succeq0.
\]

Introducing a PSD dual matrix `Y` gives

\[
\boxed{
\tau(C)=
\max\left\{
-\sum_{i<j} C_{ij}Y_{ij}:
Y\succeq0,\;Y_{ii}=1
\right\}.
}
\tag{2}
\]

Strong duality holds because the primal is strictly feasible after adding a
sufficiently large positive diagonal.

Equation (2) is an elliptope optimization.  In contrast,

\[
s(C)=\max_{x\in\{\pm1\}^6}\sum_{i<j}C_{ij}x_i x_j
\]

optimizes the same quadratic edge functional over rank-one correlation
matrices `xx^T`.  The remaining theorem is therefore naturally a comparison
between a cut-type support function and an elliptope support function, coupled
through the cubic perfect-matching polynomial `haf(C)`.

This is a sharper description of the phrase "the gap is caused by PSD
geometry."  The PSD cost is exactly `tau(C)`.

## The two known equality geometries

The reduction preserves both important equality structures.

### Disjoint pairs

For the normalized three-pair matrix, the three nonzero entries of `C` are
`1/6`.  The diagonal `d_i=1/12` gives

\[
\tau(C)=\frac12,
\qquad
s(C)=\frac12,
\qquad
\operatorname{haf}(C)=\frac1{216}.
\]

A matching block dual certificate with blocks

\[
\begin{pmatrix}1&-1\\-1&1\end{pmatrix}
\]

has dual value `1/2`, proving that the completion trace is minimal.

### Rank one

For the equal all-parallel point, every off-diagonal coefficient is `1/18`.
The rank-one completion is `B=J_6/36`, so

\[
\tau(C)=\frac16,
\qquad
s(C)=\frac56,
\qquad
|\operatorname{haf}(C)|=\frac5{1944}.
\]

The dual certificate

\[
Y=\frac65\left(I-\frac16J\right)
\]

has unit diagonal, is PSD, and has objective value `1/6`.

Hence both equality families of the candidate inequality are also exact
minimum-trace completion points.

## Exact PSD cost of the non-PSD `1/200` witness

The bounded degree-two witness from
`experiments/e11_relaxed_six_variable_witness.py` has

\[
C=\frac1{10}S,
\]

where `S` is a symmetric zero-diagonal sign matrix satisfying

\[
S^2=5I.
\]

The relaxation uses trace `1/2`, equivalently diagonal `1/12`, but that matrix
is not PSD.

Let

\[
d=\frac{\sqrt5}{20}.
\]

Since

\[
\left(\frac C2\right)^2=d^2I,
\]

the equal-diagonal matrix

\[
B_\star=dI+\frac C2
\]

is PSD.  Its trace is

\[
6d=\boxed{\frac{3\sqrt5}{10}}\approx0.6708203932.
\]

This trace is exactly minimal.  Indeed

\[
Y=I-\frac{C/2}{d}
\]

has unit diagonal and eigenvalues `0` and `2`, hence is dual feasible.  Its
dual objective equals

\[
-\left\langle\frac C2,Y\right\rangle
=6d
=\frac{3\sqrt5}{10}.
\]

Thus primal and dual values coincide.

For this witness,

\[
s(C)=\frac12,
\qquad
|\operatorname{haf}(C)|=\frac1{200}.
\]

Its ratio in the completion-reduced candidate inequality is only

\[
\frac{54/200}{
(3\sqrt5/10)(1/2)(3\sqrt5/10+1/2)}
\approx0.68754<1.
\]

So the known `1/200` relaxation witness is not merely "non-PSD": its exact PSD
diagonal-completion cost can be certified in closed form, and that cost moves
it well inside the conjectured feasible bound.

## Consequence for the next proof step

The six-variable problem can now be attacked without carrying six free
diagonal variables.  A proof of (1) may use either side of (2):

1. construct a dual elliptope certificate `Y` from the edge pattern `C` and
   lower-bound `tau(C)` in terms of the matching polynomial;
2. characterize extremal minimum-trace completions using complementary
   slackness `(diag(d)+C/2)Y=0`;
3. combine the cut maximizer defining `s(C)` with a dual correlation matrix
   defining `tau(C)`.

This is also a natural formulation for a small exact SDP/SOS certificate if a
human-readable inequality does not close the final gap.
