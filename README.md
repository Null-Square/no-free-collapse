# No Free Collapse

**Interaction-order limits of quantum-inspired reasoning.**

This repository develops exact, CPU-verifiable mathematics for a narrow question:

> If many reasoning fragments are encoded into a quantum-like latent state and an answer is produced by Born-style collapse, which higher-order reasoning interactions must already have been created before measurement?

The working principle is **no free collapse**: measurement can expose interference already present in the prepared state, while nonlinear state preparation and normalization must be counted explicitly as computational resources.

## What is verified now

The codebase currently establishes eight facts.

1. **Quadratic ceiling for fixed-norm polynomial preparations.** If each amplitude has Boolean interaction order at most `r` and the state has input-independent norm, a fixed Born effect has output interaction degree at most `2r`.
2. **Normalization loophole.** Input-dependent normalization is nonlinear and can create higher-order interactions. An exact order-1 three-variable example has cubic coefficient `-6/65`.
3. **Normalization explosion.** The loophole is unbounded in exact degree: an affine/order-1 unnormalized two-dimensional state can acquire a nonzero full `n`-way interaction after normalization for every `n`.
4. **Conditioning-controlled leakage.** Bounded norm variation forces normalization-induced high-order interactions to decay exponentially.
5. **Chebyshev/Helstrom bound.** If `m=floor((k-1)/(2r))` and `rho=(sqrt(kappa)-1)/(sqrt(kappa)+1)`, then every order-`r` model obeys
   `|p_hat(S)| <= rho^m/(1+rho^(2m))` for every order-`k` Walsh interaction, even after choosing the optimal Born measurement.
6. **Tight fixed-norm parity construction.** A 2D, exactly normalized state with preparation order `ceil(k/2)` computes a pure `k`-way parity interaction exactly, attaining the factor-of-two fixed-norm bound.
7. **Exact Gram reduction.** Every order-`r` model is a constrained rational Gram quotient `z*Az / z*Qz` with `0 <= A <= Q`. Conversely every such pair is realizable, with minimum latent dimension `rank(Q)`. For fixed `Q`, the best final Born measurement for any linear interaction objective is obtained exactly by eigendecomposition.
8. **Exact symmetric-class optimum.** For even `n>=4`, among all real order-1 Gram preparations invariant under coordinate permutations and global sign reversal, with norm condition number at most `kappa`, the exact optimum full-parity interaction is the mean-field construction

   `q(x)=1+(kappa-1)(sum_i x_i)^2/n^2`,

   with closed-form capacity

   `C_sym = kappa/(kappa-1) * n! / [2^n prod_{m=1}^{n/2}(m^2+n^2/(4(kappa-1)))]`.

The normalization/rational-function connection itself is established prior art: rational approximation is tightly connected to quantum query algorithms with postselection. Chebyshev approximation and Helstrom discrimination are also classical ingredients. We do **not** claim those components separately as novel. The research target is their resource-theoretic synthesis for collapse-based reasoning and the remaining unrestricted extremal problem.

## Run

```bash
python -m pip install -e '.[dev]'
pytest
python experiments/e1_linear_barrier.py
python experiments/e2_normalization_loophole.py
python experiments/e3_tight_parity.py
python experiments/e4_normalization_explosion.py
python experiments/e5_conditioned_leakage.py
python experiments/e6_optimal_collapse.py
python experiments/e7_chebyshev_vs_mean_field.py
python experiments/e8_asymmetric_search.py
```

Everything runs on CPU and small Boolean hypercubes.

## Current hard problem

For a subset `S` of size `k`, define

\[
C_{n,r,k}(\kappa,d)=\sup |\widehat p(S)|,
\]

where the supremum ranges over physical order-`r` Born models with squared-norm condition number at most `kappa` and latent Gram rank at most `d`.

The current dimension-free Chebyshev/Helstrom bound is

\[
C_{n,r,k}(\kappa,d)
\le
\frac{\rho^m}{1+\rho^{2m}},
\qquad
\rho=\frac{\sqrt\kappa-1}{\sqrt\kappa+1},
\quad
m=\left\lfloor\frac{k-1}{2r}\right\rfloor.
\]

For even full parity at `r=1`, the symmetric class is solved exactly. The remaining high-value question is now sharper:

> Can a symmetry-breaking real or complex Gram matrix beat the exact symmetric capacity?

Random asymmetric searches have not done so in the tested small cases, but the unrestricted objective is neither convex nor concave, so group averaging is not a valid proof. We therefore treat global optimality as open until there is either a proof or an explicit counterexample.

## Notes

- [`docs/math.md`](docs/math.md): exact interaction-order results and normalization witnesses.
- [`docs/conditioning.md`](docs/conditioning.md): geometric de-normalization and spectral decay.
- [`docs/chebyshev.md`](docs/chebyshev.md): minimax Chebyshev/Helstrom interaction bound.
- [`docs/gram.md`](docs/gram.md): Gram characterization and exact optimal collapse for fixed preparation.
- [`docs/symmetric.md`](docs/symmetric.md): exact mean-field construction and asymptotics.
- [`docs/symmetric_optimality.md`](docs/symmetric_optimality.md): proof of optimality in the full permutation/global-sign symmetric class.
