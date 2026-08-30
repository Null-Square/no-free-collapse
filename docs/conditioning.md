# Conditioning-controlled normalization

Let an unnormalized Boolean-order-`r` amplitude preparation be

\[
\tilde\psi(x)=\sum_{|S|\le r}x_Su_S,
\]

with a fixed effect `0 <= M <= I`. Define

\[
a(x)=\langle\tilde\psi,M\tilde\psi\rangle,\qquad q(x)=\|\tilde\psi\|^2,\qquad p(x)=a(x)/q(x).
\]

Both `a` and `q` have degree at most `2r`. Thus input-dependent normalization produces a rational function. This rational/postselection connection is established prior art (Mahadev and de Wolf, 2015); the question here is the conditioned version.

## Conditioning parameter

Let

\[
\kappa=q_{\max}/q_{\min},\qquad q_c=(q_{\max}+q_{\min})/2,
\]

and

\[
h(x)=q(x)/q_c-1,\qquad |h(x)|\le\delta=\frac{\kappa-1}{\kappa+1}<1.
\]

## Theorem: geometric de-normalization

For every integer `m>=0`, there is a polynomial `P_m` of degree at most `2rm` such that

\[
\|p-P_m\|_\infty\le\delta^m.
\]

For `m>=1`, take

\[
P_m=\frac{a}{q_c}\sum_{t=0}^{m-1}(-h)^t.
\]

The remainder is exact:

\[
p-P_m=\frac{a}{q_c}\frac{(-h)^m}{1+h}=p(-h)^m.
\]

Since `0<=p<=1`, the claimed bound follows immediately.

This cleanly separates exact degree from interaction strength: normalization can create full degree, but if the norm is well-conditioned, the new high-order component must be small.

## Spectral corollaries

For the Walsh expansion `p(x)=sum_S p_hat(S)x_S`, orthogonal projection gives

\[
\sum_{|S|>2rm}|\widehat p(S)|^2\le\delta^{2m}.
\]

For an individual order-`k` interaction, choose

\[
m=\left\lfloor\frac{k-1}{2r}\right\rfloor.
\]

Then

\[
|\widehat p(S)|\le\delta^m,\qquad |S|=k.
\]

Inverting the bound, an observed coefficient of magnitude `alpha` requires

\[
\kappa\ge\frac{1+\alpha^{1/m}}{1-\alpha^{1/m}}.
\]

## Prior-art boundary

The representation `p=a/q` as a low-degree rational function is not a novelty claim. Rational approximation is tightly connected to quantum query algorithms with postselection; see:

- U. Mahadev and R. de Wolf, *Rational approximations and quantum algorithms with postselection*, Quantum Information and Computation 15 (2015), 295-307.

Our current target is the extremal conditioned problem: determine the tight interaction spectrum achievable by physical Born-rational models as a function of preparation order, norm condition number, and latent Gram rank.

Define

\[
C_{n,r,k}(\kappa,d)=\sup |\widehat p(S)|,\qquad |S|=k,
\]

over order-`r` models with condition number at most `kappa` and Gram rank at most `d`. The theorem above gives a dimension-free upper bound, but current experiments show it is loose. Tightening this capacity is the next hard problem.
