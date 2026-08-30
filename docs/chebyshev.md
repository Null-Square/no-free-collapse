# Chebyshev no-free-collapse bound

The geometric de-normalization argument in `conditioning.md` is correct but not minimax-optimal. A stronger universal bound follows by choosing the best residual polynomial on the norm interval.

Let an order-`r` preparation have squared norm

\[
q(x)=\|\tilde\psi(x)\|^2\in[q_{\min},q_{\max}],
\qquad \kappa=q_{\max}/q_{\min}.
\]

For a Walsh character `chi_S` of order `k>0`, set

\[
m=\left\lfloor\frac{k-1}{2r}\right\rfloor.
\]

## Minimax residual

Among degree-`m` polynomials `R` satisfying `R(0)=1`, the smallest possible uniform norm on `[q_min,q_max]` is attained by a shifted/scaled Chebyshev polynomial:

\[
\epsilon_m(\kappa)
=\frac{1}{\left|T_m\!\left(\frac{\kappa+1}{\kappa-1}\right)\right|}.
\]

Writing

\[
\rho=\frac{\sqrt\kappa-1}{\sqrt\kappa+1},
\]

this is

\[
\epsilon_m(\kappa)=\frac{2\rho^m}{1+\rho^{2m}}.
\]

Because `R(0)=1`, there is a polynomial `P` of degree at most `m-1` with

\[
R(q)=1-qP(q).
\]

## Operator cancellation theorem

Define normalized states

\[
|\phi_x\rangle=|\tilde\psi_x\rangle/\sqrt{q(x)}
\]

and the signed state moment

\[
\Delta_S=\mathbb E_x\left[\chi_S(x)|\phi_x\rangle\langle\phi_x|\right].
\]

Each entry of

\[
|\tilde\psi_x\rangle\langle\tilde\psi_x|P(q(x))
\]

has Boolean degree at most `2rm<k`, so its `S`-Fourier coefficient vanishes. Therefore

\[
\Delta_S
=\mathbb E_x\left[\chi_S(x)|\phi_x\rangle\langle\phi_x|R(q(x))\right].
\]

Since each normalized rank-one projector has trace norm one,

\[
\|\Delta_S\|_1\le \epsilon_m(\kappa).
\]

For any two-outcome Born effect `0<=M<=I`, the interaction coefficient is

\[
\widehat p(S)=\operatorname{Tr}(M\Delta_S).
\]

Also `Tr(Delta_S)=E chi_S=0`. By the Helstrom variational identity,

\[
\sup_{0\preceq M\preceq I}|\operatorname{Tr}(M\Delta_S)|
=\frac12\|\Delta_S\|_1.
\]

Hence

\[
\boxed{
|\widehat p(S)|
\le
\frac{1}{2\left|T_m\!\left(\frac{\kappa+1}{\kappa-1}\right)\right|}
=
\frac{\rho^m}{1+\rho^{2m}}
}.
\]

For fixed `kappa>1` and large interaction order,

\[
|\widehat p(S)|\lesssim
\exp\left[-\frac{k}{2r}\log\frac{\sqrt\kappa+1}{\sqrt\kappa-1}\right].
\]

Thus bounded norm conditioning forces exponential decay of interactions whose order greatly exceeds preparation order, even when the final measurement is chosen optimally.

## Novelty discipline

None of the following ingredients is claimed as new on its own:

- Chebyshev minimax approximation;
- rational approximation / postselected quantum-query complexity;
- Fourier analysis on the Boolean cube;
- Helstrom binary state discrimination.

The research object is their combination as a resource bound for collapse-based reasoning, plus the extremal capacity problem under preparation order, conditioning, and Gram-rank constraints.
