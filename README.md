# No Free Collapse

**Interaction-order limits of quantum-inspired reasoning.**

This repository develops exact, CPU-verifiable mathematics for a narrow question:

> If many reasoning fragments are encoded into a quantum-like latent state and an answer is produced by Born-style collapse, which higher-order reasoning interactions must already have been created before measurement?

The working principle is **no free collapse**: measurement can expose interference already present in the prepared state, while nonlinear state preparation and normalization must be counted explicitly as computational resources.

## What is verified now

The codebase currently establishes six facts.

1. **Quadratic ceiling for fixed-norm polynomial preparations.** If each amplitude has Boolean interaction order at most `r` and the state has input-independent norm, a fixed Born effect has output interaction degree at most `2r`.
2. **Normalization loophole.** Input-dependent normalization is nonlinear and can create higher-order interactions. An exact order-1 three-variable example has cubic coefficient `-6/65`.
3. **Normalization explosion.** The loophole is unbounded in exact degree: an affine/order-1 unnormalized two-dimensional state can acquire a nonzero full `n`-way interaction after normalization for every `n`.
4. **Conditioning-controlled leakage.** If the squared norm has condition number `kappa`, write `delta=(kappa-1)/(kappa+1)`. An order-`r` normalized Born output has a degree-`2rm` polynomial approximation with uniform error at most `delta^m`; therefore its Walsh interaction energy above order `2rm` is at most `delta^(2m)`.
5. **Tight fixed-norm parity construction.** A 2D, exactly normalized state with preparation order `ceil(k/2)` computes a pure `k`-way parity interaction exactly, attaining the factor-of-two fixed-norm bound.
6. **Exact Gram reduction.** Every order-`r` model is a constrained rational Gram quotient `z*Az / z*Qz` with `0 <= A <= Q`. Conversely every such pair is realizable, with minimum latent dimension `rank(Q)`. For fixed `Q`, the best final Born measurement for any linear interaction objective is obtained exactly by eigendecomposition.

The normalization/rational-function connection itself is established prior art: rational approximation is tightly connected to quantum query algorithms with postselection. We do **not** claim that observation as novel. The current research target is the conditioned and resource-constrained problem: how much interaction strength can a low-order, low-rank, well-conditioned preparation generate after optimal collapse?

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
```

Everything runs on CPU and small Boolean hypercubes.

## Current hard problem

For a subset `S` of size `k`, define the extremal interaction capacity

\[
C_{n,r,k}(\kappa,d)=\sup |\widehat p(S)|,
\]

where the supremum ranges over physical order-`r` Born models with squared-norm condition number at most `kappa` and latent Gram rank at most `d`.

The current theorem gives the dimension-free upper bound

\[
C_{n,r,k}(\kappa,d)\le\left(\frac{\kappa-1}{\kappa+1}\right)^{\lfloor (k-1)/(2r)\rfloor},
\]

but numerical experiments indicate it is loose. The next research stage is to determine the tight scaling, first computationally and then analytically.

## Notes

- [`docs/math.md`](docs/math.md): exact interaction-order results and normalization witnesses.
- [`docs/conditioning.md`](docs/conditioning.md): approximate de-normalization and spectral-decay theorem.
- [`docs/gram.md`](docs/gram.md): Gram characterization and closed-form optimal collapse for fixed preparation.
