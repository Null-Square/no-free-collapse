# Publication package and submission status

**12 September 2026 — complete proof drafts for independent review.**

The primary contribution is now the sharp complementary-subset dilation, not another projection stratum. It proves the full real six-variable PSD-gradient theorem, previously results-ledger E4/E5, and classifies all PSD equality cases. The separate Boolean-cube-normalized `1/216` optimum remains open.

## Primary manuscript

[`matching_dilation.tex`](matching_dilation.tex), titled **A complementary-subset dilation for sharp hafnian gradient inequalities**, is a complete six-page proof draft. The rendered PDF and source are also included in the downloadable research package supplied with this work.

For a real symmetric matrix A of order 2m, the weighted Kneser operator on (m-1)-subsets satisfies

`||T_m(A)||_op <= (m+1)/4 * (lambda_max(A)-lambda_min(A))`.

The coefficient is sharp for every m>=2. The proof gives an exact complementary-subset dilation and two explicit PSD slack decompositions. In order six it implies `q2(A)<=d(A)^2 q1(A)/4` and `|haf(A)|<=d(A)q1(A)/6`. Therefore the previously open PSD-contraction theorem holds for all real PSD contractions, including non-idempotent inputs.

For nonzero PSD A with L=lambda_max(A), the nondegenerate equality cases have A/L equal to `ss^T/6`, `I-ss^T/6`, or `(I+M)/2`, where s is a sign vector and M is a signed perfect-matching matrix. Diagonal PSD matrices give the degenerate scalar equalities. Necessity and sufficiency are proved, not only tested on examples.

The full proof note is [`docs/matching_dilation.md`](../docs/matching_dilation.md). The current theorem ledger is [`docs/RESULTS.md`](../docs/RESULTS.md), which supersedes historical open-problem status in older progress notes.

## Companion contribution

[`docs/zero_diagonal_spectral.md`](../docs/zero_diagonal_spectral.md) proves a complex zero-diagonal spectral bound through a four-term nonnegative defect identity. It includes equality classification, phase-preserving matching recovery, and the six-variable guarantee `||B-M||_F^2<=34 epsilon` when `||B||<=1`, `|haf(B)|>=1-epsilon`, and `epsilon<=1/48`. An exact rational rotation family establishes optimal linear order in squared distance; the constant 34 is not asserted optimal.

The seven-page companion manuscript, **Sharp spectral bounds and matching stability for four-subhafnian energy**, is supplied as LaTeX and PDF in the downloadable research package. Its stability statements are not a quantitative classification of all three real PSD equality families. The two drafts can be considered together during publication planning; two separate papers are not assumed necessary.

## Implemented and verified

- `src/no_free_collapse/matching_dilation.py`: direct/dilation operators, explicit PSD slacks, exact rational six-variable certificates. Explicit numerical enumeration is capped at order ten to avoid exponential allocations; the theorem has no dimension cap.
- `src/no_free_collapse/zero_diagonal.py`: direct and trace energy, nonnegative defect, rational spectral certificates, and constructive matching recovery.
- 79 new research tests, including integer checks of the full dilation map on all 122 symmetric matrix basis elements in orders 4,6,8,10, exact sparse-polynomial identities, rational PSD instances, equality cases, boundary rejection, and optimal-order stability witnesses.
- Full repository CI at commit `61e71a5f1a291adcb8596f708db19a068a15a1ad`, run `34689783115`: **200 passed on Python 3.10**; the Python 3.12 job also succeeded. The local research overlay separately passed all 79 new tests on Python 3.13.5 / NumPy 2.3.5.
- Seeded experiments: 4,500 complex/real zero-diagonal matrices, maximum identity residual `3.9968028886505635e-15`; 400 additional real symmetric matrices, maximum dilation residual `2.220446049250313e-16`; 20 exact rational PSD instance certificates.

Finite checks protect implementation consistency. The analytic proofs establish the universal theorems. Floating-point tolerances are not described as rigorous enclosures. Runtime reproduction uses only NumPy and the Python standard library; pytest is a development dependency.

## Reproduce

In the full repository:

```sh
python -m pip install -e '.[dev]'
python -m pytest
python experiments/e20_spectral_stability.py --samples 250
python experiments/e21_matching_dilation.py --samples 100
```

The standalone downloadable artifact is a new-research overlay, not a complete clone; its README uses `PYTHONPATH=src` and runs the 79 included tests. Actual logs, JSON results, PDFs, and integrity hashes accompany that artifact.

## Corrections and claim boundaries

The rank-one equality statement now includes coordinate projections with `q1=q2=0`. More importantly, the old statement that the PSD operator bound was false is withdrawn. The conference involution has `||T(K)||=3/sqrt(5)>1`, while its PSD projection `(I+K)/2` has half that matching-operator norm. The counterexample concerned the differently scaled involution bound. The new dilation proves the PSD statement, and a regression test preserves this distinction.

The spectral theorem bypasses, rather than proves, positivity of every mixed spectral coefficient. The cube-normalized `1/216` optimum, the sharper global rank-two ratio problem, and stability around all three PSD equality families are not claimed solved.

## Before external submission

The proposed contribution is explicit and has complete proofs, code, and reproducibility evidence. It has not received independent mathematical review or exhaustive priority certification. The manuscripts compare specific results of Roos and Kummer; the full Lindell–Staples zeon-norm paper was unavailable during this pass, so a full-text specialist comparison remains necessary.

The owner must confirm human authors, affiliations, author order, corresponding-author details, assistance disclosure, and applicable declarations, obtain independent proof review, complete the priority comparison, and apply the chosen venue's submission requirements. No authorship, funding, conflict statement, journal acceptance, or submission is invented here. The PR remains a draft and has not been merged.
