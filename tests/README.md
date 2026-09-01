# Tests

The pytest suite protects the implementation of the mathematical results in this repository.

## What tests do

Tests check:

- independent formulas for the same exact quantity;
- equality constructions and boundary cases;
- exact rational / Bernstein certificates used in proofs;
- deterministic regression cases;
- fixed-seed random projections/frames that catch indexing or algebraic mistakes.

## What tests do not do

A finite random test is not treated as proof of a continuum inequality. The analytic proof notes live in [`../docs/`](../docs/), and claim status is recorded in [`../docs/RESULTS.md`](../docs/RESULTS.md).

The intended evidence hierarchy is:

1. analytic proof;
2. exact finite certificate when a proof reduces to polynomial positivity;
3. deterministic identity cross-check;
4. seeded random regression;
5. exploratory numerical search.

Only levels 1–2 establish theorem status.

## Useful targeted groups

Core resource results:

```bash
python -m pytest tests/test_degree_bound.py tests/test_normalization_loophole.py tests/test_conditioning.py tests/test_chebyshev.py tests/test_gram.py
```

Projection theorem:

```bash
python -m pytest tests/test_projection_gradient.py tests/test_rank_three_global_gradient.py tests/test_rank_two_global_projection.py
```

Spectral frontier:

```bash
python -m pytest tests/test_spectral_polarization.py tests/test_mixed_rank_1555.py
```

See [`../docs/REPRODUCIBILITY.md`](../docs/REPRODUCIBILITY.md) for the complete reviewer verification guide.
