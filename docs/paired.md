# Matched-pair symmetry breaking

The permutation-symmetric mean-field construction is not globally optimal. A cleaner symmetry-breaking family can be solved exactly and is optimal in a broad matched pair-block class.

Let `n=2L>=4` and partition the variables into pairs. Restrict the variable Gram matrix to independent swap-symmetric blocks

\[
Q_j=\begin{pmatrix}u_j&v_j\\v_j&u_j\end{pmatrix}.
\]

After flipping one variable inside a pair if necessary, orient the block so its opposite-sign eigenvalue is at least its same-sign eigenvalue. Write

\[
a_j=\lambda_{j,+}\ge0,\qquad d_j=\lambda_{j,-}-\lambda_{j,+}\ge0.
\]

For pair-parity bits `y_j=x_{2j-1}x_{2j}`, the squared norm is

\[
q(y)=2\sum_j(a_j+d_j\mathbf 1\{y_j=-1\}).
\]

Scale so `sum_j a_j=1`. Then a condition-number budget `kappa` is exactly

\[
\sum_j d_j\le\kappa-1.
\]

## Exact capacity inside the class

Averaging over the within-pair sign choices diagonalizes the signed Helstrom operator. Its positive contribution from pair `j` is

\[
2^{-L}a_j\int_0^\infty e^{-t}
\prod_{\ell\ne j}(1-e^{-d_\ell t})\,dt.
\]

Therefore the exact optimal full-parity Born coefficient for fixed pair parameters is

\[
C(a,d)=2^{-L}\sum_j a_jF_j(d),
\qquad
F_j(d)=\int_0^\infty e^{-t}\prod_{\ell\ne j}(1-e^{-d_\ell t})dt.
\]

This immediately solves the class optimization.

1. For fixed `d`, the objective is linear in `a`, so all baseline norm is placed on one pair: the **readout pair**.
2. `F_j` does not depend on the variation `d_j` of that same readout pair, so assigning condition budget there is useless. Set `d_j=0`.
3. The integrand is increasing in every remaining `d`, so use the full budget `sum d=kappa-1`.
4. For every fixed `t>0`, `log(1-e^{-dt})` is concave in `d>0`. Jensen therefore makes the product maximal when the remaining `L-1` variations are equal.

Thus the optimum matched-pair preparation has one isotropic readout block and `L-1` equal rank-one normalizer blocks. In the feature order `[1,x1,...,xn]`, choose

\[
Q_{\rm read}=I_2,
\qquad
Q_{\rm norm}=\frac d2
\begin{pmatrix}1&-1\\-1&1\end{pmatrix},
\qquad
d=\frac{\kappa-1}{L-1}.
\]

It has cube condition number exactly `kappa` and rank

\[
L+1=n/2+1.
\]

The exact capacity is

\[
\boxed{
C_{\rm pair}(n,\kappa)
=2^{-L}\frac{\Gamma(1+\alpha)\Gamma(L)}{\Gamma(L+\alpha)},
\qquad
\alpha=\frac{L-1}{\kappa-1}.
}
\]

The formula follows from substituting `u=e^{-dt}` in

\[
2^{-L}\int_0^\infty e^{-t}(1-e^{-dt})^{L-1}dt.
\]

For `n=4`,

\[
C_{\rm pair}(4,\kappa)=\frac{\kappa-1}{4\kappa}.
\]

At `kappa=2`, this is exactly `1/8`, while the permutation-symmetric optimum is `3/40=0.075`. Hence symmetry breaking gives a strict analytic counterexample to unrestricted mean-field optimality.

## Asymptotic rate

For fixed `kappa>1`, Stirling gives

\[
C_{\rm pair}(n,\kappa)=\exp[-nI_{\rm pair}(\kappa)+O(\log n)],
\]

with

\[
\boxed{
I_{\rm pair}(\kappa)=\frac12\left[
\log2+(1+c)\log(1+c)-c\log c
\right],\qquad c=\frac1{\kappa-1}.
}
\]

The pair family has a better asymptotic exponent than mean-field at low conditioning; mean-field becomes better at sufficiently large conditioning. Numerically their rate crossover is near `kappa=5.31806`.

## Remaining problem

This theorem proves optimality only inside the matched pair-block class. Together with the symmetric-class theorem it shows that the unrestricted capacity has nontrivial competing extremizers. The next target is a structural lower/upper theory for arbitrary quadratic normalizers rather than assuming a single symmetry class.
