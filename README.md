# No Free Collapse

**Complementary-subset dilations and a spectral hierarchy for hafnian energies.**

[![CI](https://github.com/Null-Square/no-free-collapse/actions/workflows/ci.yml/badge.svg)](https://github.com/Null-Square/no-free-collapse/actions/workflows/ci.yml)

This repository contains analytic proofs, theorem-supporting code, exact rational certificates, and reproducible diagnostics. The current paper is a focused matrix-inequalities contribution. The earlier normalized-model research is retained as background, not as a claim of quantum computational advantage.

**Publication status: author-review candidate. Not submitted.** Technical checks and human submission approval are separate. See the [submission checklist](submission/README.md) and run `python tools/submission_check.py`.

## Main results

For real symmetric A of order 2m, m >= 2, the specified weighted Kneser operator obeys the sharp bound

\[
\|T_m(A)\|_{op}\le\frac{m+1}{4}\bigl(\lambda_{max}(A)-\lambda_{min}(A)\bigr).
\]

Writing `d(A)=lambda_max(A)-lambda_min(A)` and

\[
\overline H_{2r}(A)=\binom N{2r}^{-1}\sum_{|S|=2r}|\operatorname{haf}(A[S])|^2,
\]

we prove

\[
\boxed{\overline H_{2r}(A)\le\tfrac14 d(A)^2\overline H_{2r-2}(A)}
\]

for **integer r >= 2 and N >= 4r-2**. The one-step constant is sharp at every order across the admissible dimensions. This is not a claim of the best constant at every fixed larger N, nor a sharp iterated product.

| Result | Status | Where to read |
|---|---|---|
| Sharp all-even-dimensional weighted-operator bound | Analytic proof | [Dilation](docs/matching_dilation.md) |
| Principal subhafnian hierarchy and Roos hybrid | Analytic proof, including odd ambient order | [Hierarchy](docs/subhafnian_hierarchy.md), [audit](docs/SUBMISSION_AUDIT.md) |
| Six-variable PSD gradient bound, formerly E4/E5 | **Proved globally**, with PSD equality classification | [Manuscript](paper/matching_dilation.tex) |
| Sharp fixed-N spectral edge-density bound | Analytic proof, including the odd-N improvement | [Audit](docs/SUBMISSION_AUDIT.md) |
| Real, pure, zero-displacement collision-free GBS consequence | Corollary with explicit restrictions | [Manuscript](paper/matching_dilation.tex) |
| Complex zero-diagonal matching stability | Companion theorem, separate hypotheses | [Companion](docs/zero_diagonal_spectral.md) |
| Cube-normalized six-variable 1/216 optimum | **Open**; not a prerequisite for this paper | [Ledger](docs/RESULTS.md) |

## Reviewer quick start

```sh
git clone https://github.com/Null-Square/no-free-collapse.git
cd no-free-collapse
python -m venv .venv
# Activate the virtual environment for your shell.
python -m pip install -e '.[dev]'
python -m compileall -q src experiments tests tools
python -m pytest
python tools/submission_check.py
python experiments/e20_spectral_stability.py --samples 250
python experiments/e21_matching_dilation.py --samples 100
python experiments/e22_subhafnian_hierarchy.py
```

The first commands run the full repository suite. Exact instance certificates verify their spectral hypotheses; random tests are diagnostics, never universal proofs. NumPy is the only numerical runtime dependency. Explicit enumerators are capped to prevent accidental exponential work.

Build the paper with a LaTeX installation:

```sh
cd paper
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error matching_dilation.tex
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error matching_dilation.tex
```

## Reading order and repository map

Start with the [paper overview](paper/README.md), then the [current theorem ledger](docs/RESULTS.md), [proof audit](docs/SUBMISSION_AUDIT.md), and [reproducibility guide](docs/REPRODUCIBILITY.md). [Paper materials](docs/PAPER_MATERIALS.md) describes the frozen submission scope. [Literature review](submission/literature_review.md) distinguishes direct comparisons from unresolved priority checks.

`src/no_free_collapse/` contains numerical diagnostics and exact certificates; `tests/` contains regression and algebraic checks; `experiments/` contains seeded reproductions; `paper/` contains manuscript sources; `submission/` contains the cover-letter draft, highlights, declarations, and approval record.

## Historical notes and claim boundaries

Earlier rank-by-rank and spectral-polarization notes are preserved as research history. Any older statement that the PSD-gradient extension is open, or that the PSD matching-operator bound is false, is superseded by the current ledger and dilation proof. The former conference counterexample concerned an involution; its PSD projection has half the matching-operator norm.

The normalization/conditioning proofs remain in [math](docs/math.md), [conditioning](docs/conditioning.md), [Gram](docs/gram.md), and [paired models](docs/paired.md). Standard degree doubling, Fourier analysis, SDP duality, and zeon algebras are not claimed as inventions of this project. High-degree probabilities alone do not demonstrate decision expressivity or computational advantage.

No exhaustive priority certification, independent human proof review, journal acceptance, or submission is claimed. The release preflight deliberately fails until author identity, declarations, proof/priority approval, and venue checks are recorded.

## License and citation

Repository code retains the [MIT license](LICENSE). See [CITATION.md](CITATION.md) for citing a version without inventing a paper DOI or confirmed author list.
