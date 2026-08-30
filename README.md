# No Free Collapse

**Interaction-order limits of quantum-inspired reasoning.**

This repository develops exact, CPU-verifiable mathematics for a narrow question:

> If many reasoning fragments are encoded into a quantum-like latent state and an answer is produced by Born-style collapse, which higher-order reasoning interactions must already have been created before measurement?

The working principle is **no free collapse**: measurement can expose interference already present in the prepared state, but it does not make arbitrary global reasoning structure appear for free.

## What is verified now

The initial codebase establishes four facts.

1. **Quadratic ceiling for fixed-norm polynomial preparations.** If each amplitude has Boolean interaction order at most `r` and the state has input-independent norm, a fixed Born effect has output interaction degree at most `2r`.
2. **Normalization loophole.** Input-dependent normalization is nonlinear and can create higher-order interactions. An exact order-1 three-variable example has cubic coefficient `-6/65`. Any rigorous resource accounting must therefore charge normalization to preparation.
3. **Normalization explosion.** The loophole is unbounded: an affine/order-1 unnormalized two-dimensional state can acquire a nonzero full `n`-way interaction after input-dependent normalization for every `n`.
4. **Tight parity construction.** A 2D, exactly normalized state with preparation order `ceil(k/2)` computes a pure `k`-way parity interaction exactly, attaining the factor-of-two bound.

These are foundational checks, not yet a publication-level novelty claim. The established quantum polynomial method already contains the underlying degree machinery. The research target is the reasoning-specific resource theory built on top of it: approximation bounds, normalization/nonlinearity costs, sparse higher-order objectives, and comparisons with unrestricted classical continuous latent states.

## Run

```bash
python -m pip install -e '.[dev]'
pytest
python experiments/e1_linear_barrier.py
python experiments/e2_normalization_loophole.py
python experiments/e3_tight_parity.py
python experiments/e4_normalization_explosion.py
```

Everything in the initial suite runs on a CPU and small Boolean hypercubes.

## Research discipline

We will distinguish carefully between:

- **raw Born score:** `psi* M psi`, whose degree is directly bounded by amplitude degree;
- **physical Born probability for an already normalized/fixed-norm state:** same bound;
- **probability after input-dependent normalization:** a rational computation that can introduce higher-order interactions and must be counted as a preparation resource.

This distinction is central to avoiding a false theorem.

See [`docs/math.md`](docs/math.md) for the mathematical core.
