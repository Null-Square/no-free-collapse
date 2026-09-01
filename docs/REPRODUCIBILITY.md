# Reproducibility and verification

This repository is designed so that the mathematical artifacts can be inspected and tested on a standard CPU-only Python environment.

## 1. Environment

The package requires Python `>=3.10` and NumPy. The development test dependency is pytest.

Recommended setup:

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -e '.[dev]'
```

The GitHub Actions workflow tests the declared Python floor and the primary current runtime.

## 2. Full verification

Run the complete test suite:

```bash
python -m pytest
```

The suite covers:

- Boolean interaction-degree calculations;
- the normalization loophole construction;
- conditioning and Chebyshev bounds;
- Gram realization and optimal collapse identities;
- symmetric and matched-pair exact formulas;
- hafnian identities and four-variable bounds;
- six-variable local first/second-order certificates;
- rank-one, diagonal-completion, and projection-gradient identities;
- exact rank-three and rank-two defect reductions;
- rational Bernstein certificates used in the global rank-two proof;
- spectral-chain reconstruction for PSD contractions;
- the mixed `(1,5,5,5)` nested spectral theorem.

A passing test suite is evidence that the implementation and stored certificates agree with the written formulas. The analytic proofs remain in `docs/`.

## 3. Curated reproduction script

Run:

```bash
python experiments/reproduce_core.py
```

This executes the deterministic/representative experiment scripts used to expose the main phenomena in paper order:

1. linear/fixed-norm interaction barrier;
2. normalization loophole;
3. tight parity construction;
4. normalization-driven high-degree growth;
5. conditioning-controlled leakage;
6. optimal final collapse;
7. Chebyshev versus mean-field comparison;
8. universal hafnian bound illustration;
9. relaxed six-variable witness and PSD-completion issue;
10. rank-three defect scan.

The exploratory asymmetric-search and pair-breaking scripts can be included with:

```bash
python experiments/reproduce_core.py --include-search
```

Those exploratory scripts are useful for discovery history, but no theorem is justified by their search output.

## 4. Exact versus floating-point checks

The project deliberately uses three kinds of verification.

### Exact algebra / rational certificates

Where a proof ends in a finite polynomial certificate, the tests use exact rational data whenever practical. Examples include the rank-two Bernstein and mixed `(1,5,5,5)` certificates.

These checks protect the literal coefficients used in the proof rather than merely sampling the inequality numerically.

### Deterministic floating-point identity checks

Some exact matrix identities are evaluated numerically against independent formulas. The tolerances are chosen only for floating-point roundoff and are not treated as proof of inequality over a continuum.

Examples include:

- direct hafnian versus tensor/Pluecker formulas;
- spectral-chain reconstruction versus direct `q1/4-q2` evaluation;
- mixed-kernel closed forms versus the original kernel definition.

### Seeded random regression

Random orthogonal projections/frames are used as regression tests to catch indexing or algebraic mistakes. Seeds are fixed in the tests. These checks are diagnostics only.

## 5. What CI certifies

CI performs installation from the repository and runs the pytest suite on multiple supported Python versions. The workflow has no GPU, network service, database, or proprietary dependency.

A green CI run therefore certifies that a fresh environment can:

- install the package;
- import/compile the implementation;
- execute the complete test suite;
- verify the stored finite certificates and theorem-supporting identities.

CI does **not** claim that numerical search proves an open conjecture.

## 6. Reviewer-oriented targeted commands

### Core interaction-order results

```bash
python -m pytest \
  tests/test_degree_bound.py \
  tests/test_normalization_loophole.py \
  tests/test_conditioning.py \
  tests/test_chebyshev.py \
  tests/test_gram.py
```

### Exact solvable classes and hafnian regime

```bash
python -m pytest \
  tests/test_symmetric.py \
  tests/test_paired.py \
  tests/test_hafnian_bounds.py \
  tests/test_pair_local_optimality.py \
  tests/test_six_variable_second_order.py \
  tests/test_six_variable_rank_one.py \
  tests/test_diagonal_completion.py
```

### Completed projection theorem

```bash
python -m pytest \
  tests/test_projection_gradient.py \
  tests/test_rank_three_equal_diagonal.py \
  tests/test_rank_three_defect_identity.py \
  tests/test_rank_three_global_gradient.py \
  tests/test_rank_two_plucker.py \
  tests/test_rank_two_two_direction.py \
  tests/test_rank_two_balanced_diagonal.py \
  tests/test_rank_two_dual_high_pair.py \
  tests/test_rank_two_global_projection.py
```

### PSD-contraction reduction and mixed theorem

```bash
python -m pytest \
  tests/test_spectral_polarization.py \
  tests/test_mixed_rank_1555.py
```

## 7. Reproducibility boundary

The following distinction is important for evaluation:

- **analytic theorem**: written proof in `docs/`, with tests protecting identities/certificates;
- **exact finite certificate**: rational or symbolic coefficients checked directly;
- **numerical diagnostic**: seeded scan/search used only to guide conjectures;
- **open statement**: listed as open in [`RESULTS.md`](RESULTS.md), regardless of numerical support.

If a statement is not tagged as proved in [`RESULTS.md`](RESULTS.md), a successful experiment should not be interpreted as a proof.

## 8. Recommended archival snapshot

For a paper submission, record the exact Git commit used to generate the manuscript and figures. The paper should cite that immutable commit (or a release/tag made from it), not only the moving `main` branch.

The natural next repository step before submission is therefore to create a versioned release after the manuscript theorem set is frozen.
