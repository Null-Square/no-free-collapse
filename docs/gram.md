# Gram characterization and exact collapse optimization

The amplitude model can be reduced to a finite matrix problem that cleanly separates **preparation** from **collapse**.

Let `z_r(x)` be the vector of all Walsh monomials `x_S` with `|S|<=r`. Any order-`r` preparation can be written

\[
\tilde\psi(x)=Uz_r(x).
\]

For an effect `0<=M<=I`, define

\[
Q=U^*U,\qquad A=U^*MU.
\]

Then

\[
Q\succeq0,\qquad 0\preceq A\preceq Q,
\]

and

\[
p(x)=\frac{z_r(x)^*Az_r(x)}{z_r(x)^*Qz_r(x)}.
\]

## Converse and minimum latent dimension

Every pair `Q>=0`, `0<=A<=Q` with positive denominator on the input set has a Born realization. On the support of

\[
Q=V\Lambda V^*,
\]

choose

\[
U=\Lambda^{1/2}V^*,\qquad M=\Lambda^{-1/2}V^*AV\Lambda^{-1/2}.
\]

The inequality `0<=A<=Q` gives `0<=M<=I`. This realizes the same `(A,Q)`.

Any factorization `Q=U^*U` needs latent dimension at least `rank(Q)`, while the construction above uses exactly that many dimensions. Therefore the minimum latent dimension for a fixed Gram pair is

\[
\boxed{\operatorname{rank}(Q)}.
\]

## Exact optimal collapse for fixed preparation

Fix `Q` and a linear output objective

\[
L(p)=\sum_x w_xp(x).
\]

For example, `w_x=2^{-n}x_S` extracts one Walsh interaction coefficient. Define

\[
B_Q=\sum_x\frac{w_x}{z_x^*Qz_x}z_xz_x^*.
\]

Every feasible `A` can be written on the support of `Q` as

\[
A=Q^{1/2}CQ^{1/2},\qquad 0\preceq C\preceq I.
\]

Hence

\[
L(p)=\operatorname{Tr}(CH_Q),\qquad H_Q=Q^{1/2}B_QQ^{1/2}.
\]

The exact optimum is

\[
\max_{0\preceq A\preceq Q}L(p)=\operatorname{Tr}((H_Q)_+),
\]

attained by choosing `C` as the projector onto the positive eigenspace of `H_Q`. Replacing `w` by `-w` gives the minimum, so the maximum absolute objective is also solved exactly.

Therefore the final Born measurement never needs gradient training in our small exact experiments: for each preparation `Q`, the strongest possible collapse is obtained by eigendecomposition.

## Research consequence

The difficult optimization is isolated to the preparation Gram matrix `Q`, subject to constraints such as:

- `Q>=0`;
- positive denominator on every Boolean input;
- bounded norm condition number;
- bounded rank / latent dimension;
- fixed preparation order through the feature map `z_r`.

This gives a rigorous CPU methodology for searching for extremal constructions and testing the tightness of the conditioning bound before attempting a general proof.
