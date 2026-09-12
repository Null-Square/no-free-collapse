# Reproducibility

## Full repository

```sh
python -m pip install -e '.[dev]'
python -m compileall -q src experiments tests tools
python -m pytest
python tools/submission_check.py
python experiments/e20_spectral_stability.py --samples 250
python experiments/e21_matching_dilation.py --samples 100
python experiments/e22_subhafnian_hierarchy.py
```

The older normalized-model outputs can still be reproduced with `python experiments/reproduce_core.py`. They are background experiments, not the main paper's test of novelty.

## Exact versus numerical checks

`matching_dilation.py` checks the linear dilation on symmetric basis elements and certifies rational PSD inputs, both spectral endpoint conditions, both operator slacks, and six-variable consequences. `hierarchy_certificate.py` adds Fraction-only principal hafnians and a hierarchy certificate for a supplied rational spectral interval. It rejects floats, bools, non-symmetric inputs, invalid intervals, and unsupported dimensions.

`subhafnian_hierarchy.py` uses floating-point arithmetic for diagnostics. Its principal-hafnian enumeration is capped at N=20. The rational enumerator is capped at N=14; the explicitly assembled matching operator at N=10. These are implementation caps, not theorem restrictions. Caches are read-only and tied to the matrix's off-diagonal data. They cannot silently be reused for different input entries or replaced by an incomplete dictionary.

Exact certificates establish the supplied instance. Finite basis and polynomial checks establish the corresponding finite identities in the tested dimensions. The written analytic arguments, not random samples or a green CI badge, establish the all-dimensional theorems.

## New audit coverage

The submission pass adds rational sharp hierarchy witnesses at r=2,3,4; odd ambient orders; exact normalized convolution identities for the parity-free Roos application; all 32,768 unweighted six-vertex support patterns as a finite equality-proof regression; sharp odd-N density witnesses; zero-energy cases; invalid input orders; stale-cache and resource-limit rejection.

Before integration, the local new-research suite passed 169 tests on Python 3.13.5, NumPy 2.3.5 and pytest 9.0.2. This is a subset of the full repository suite, not an extra count to add to it. Fresh integration counts must be read from the CI run for the actual code commit, not inferred by adding historical totals.

## Paper build

From `paper/`, run pdflatex twice with `-no-shell-escape -interaction=nonstopmode -halt-on-error`. Inspect the log for unresolved references and overfull boxes, then render every page and check the layouts. The author-review candidate is deliberately marked as such until the approval record is completed.

The downloaded research archive includes its exact sources, test outputs, reproducibility data and SHA-256 manifest. Older archive files are retained under their original versioned names and are not silently represented as the new result.
