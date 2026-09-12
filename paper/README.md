# Publication package and submission status

**12 September 2026 — complete proof draft for independent review.**

The primary contribution is now broader than the original six-variable PSD-gradient problem. A sharp complementary-subset dilation yields an **all-orders spectral hierarchy for averaged principal subhafnian energies**, a quantitative hybrid with Roos's 2025 log-subadditivity theorem, and a collision-free Gaussian-boson-sampling corollary. The six-variable PSD theorem and complete equality classification remain important low-order consequences. The separate Boolean-cube-normalized `1/216` optimum remains open.

## Primary manuscript

[`matching_dilation.tex`](matching_dilation.tex), titled **Complementary-subset dilations and a spectral hierarchy for hafnian energies**, is the current main proof draft.

For a real symmetric matrix A of order 2m, the weighted Kneser operator on (m-1)-subsets satisfies

`||T_m(A)||_op <= (m+1)/4 * (lambda_max(A)-lambda_min(A))`.

The coefficient is sharp for every m>=2. The proof gives an exact complementary-subset dilation and two explicit PSD slack decompositions.

Writing

`H_{2r}(A)=sum_{|S|=2r} |haf(A[S])|^2`

and `Hbar_{2r}=H_{2r}/C(N,2r)`, the operator theorem implies, for `N>=4r-2`,

`Hbar_{2r}(A) <= d(A)^2/4 * Hbar_{2r-2}(A)`.

The factor `1/4` is sharp at every order in the uniform admissible family. At critical size `N=4r-2`, balanced signed rank-one projections, their complements, and signed matching half-rank projections attain equality.

In Roos's normalization `G_r=Hbar_{2r}/((2r-1)!!)^2`, the new recurrence is

`G_r <= d(A)^2/[4(2r-1)^2] * G_{r-1}`.

Combining it with Roos's `G_r<=G_1 G_{r-1}` gives a **one-switch hybrid bound**. The sharp scale-free parameter `theta=4G_1/d(A)^2` satisfies `theta<=1/(N-1)`; equality occurs exactly, after shift and scale, at rank-N/2 constant-diagonal 1/2 projections. At critical matching projections, the new one-step coefficient improves the Roos coefficient by `(2r-1)^2/(4r-3)`, asymptotically linear in r.

In order six the same theorem gives `q2(A)<=d(A)^2 q1(A)/4` and `|haf(A)|<=d(A)q1(A)/6`, proving the full real PSD-contraction theorem. For nonzero PSD A with `L=lambda_max(A)`, the nondegenerate equality cases have `A/L` equal to `ss^T/6`, `I-ss^T/6`, or `(I+M)/2`; diagonal PSD matrices are the degenerate scalar equality cases.

The analytic notes are [`docs/matching_dilation.md`](../docs/matching_dilation.md) and [`docs/subhafnian_hierarchy.md`](../docs/subhafnian_hierarchy.md). The current theorem ledger is [`docs/RESULTS.md`](../docs/RESULTS.md).

## Gaussian-boson-sampling consequence

For a real symmetric pure GBS sampling matrix, collision-free pattern probabilities have a common normalization times `|haf(B[S])|^2`. The hierarchy therefore bounds **unconditional average per-pattern probability** across admissible collision-free photon sectors. It is deliberately not a statement about threshold-detector Torontonians, collisionful repetitions, or monotonicity of total sector probability.

## Companion stability contribution

[`docs/zero_diagonal_spectral.md`](../docs/zero_diagonal_spectral.md) proves a complex zero-diagonal spectral bound through a four-term nonnegative defect identity. It includes equality classification, phase-preserving matching recovery, and the six-variable guarantee `||B-M||_F^2<=34 epsilon` when `||B||<=1`, `|haf(B)|>=1-epsilon`, and `epsilon<=1/48`. An exact rational rotation family establishes optimal linear order in squared distance; the constant 34 is not asserted optimal.

## Implemented and verified

- `src/no_free_collapse/matching_dilation.py`: direct/dilation operators, explicit PSD slacks, exact rational six-variable certificates.
- `src/no_free_collapse/subhafnian_hierarchy.py`: all principal hafnians for small matrices, averaged hierarchy diagnostics, Roos/spectral hybrid coefficients, one-switch closed form, and edge-density diagnostics.
- `src/no_free_collapse/zero_diagonal.py`: direct and trace energy, nonnegative defect, rational spectral certificates, and constructive matching recovery.
- The new research overlay now contains **124 targeted tests**, all passing locally after the hierarchy additions. These include all 122 symmetric matrix basis elements in orders 4,6,8,10 for the dilation map, exact sparse-polynomial identities, rational PSD certificates, direct integer checks of the hafnian recursion, exact global double-counting, critical sharp families through order eight, the one-switch hybrid, sharp `theta` equality, and stability witnesses.
- Seeded hierarchy experiment: seven `(N,r)` random regimes through `r=4`, critical sharpness formulas through `r=6`, hybrid PSD diagnostics, and real-interferometer GBS-like checks.
- Previous repository-wide CI passed 200 tests on both Python 3.10 and 3.12 before the hierarchy extension. A fresh full CI run is the final integration check for this revision.

Finite checks protect implementation consistency. The analytic proofs establish the universal theorems. Floating-point tolerances are not described as rigorous enclosures.

## Reproduce

```sh
python -m pip install -e '.[dev]'
python -m pytest
python experiments/e20_spectral_stability.py --samples 250
python experiments/e21_matching_dilation.py --samples 100
python experiments/e22_subhafnian_hierarchy.py --samples 20
```

## Prior-art boundary

The manuscript compares the new hierarchy directly with Roos's 2025 averaged squared-subhafnian inequality and distinguishes it from Kummer's ordinary Kneser construction, Lindell–Staples zeon norm inequalities, and the July 2026 Ouimet–Greaves lower bound for repeated-index PSD hafnians. The current search did not locate the matrix-weighted spectral-diameter theorem or the averaged spectral recurrence, but absence from a targeted search is not a proof of priority. A specialist full-text novelty review remains a pre-submission requirement.

## Corrections and claim boundaries

The old statement that the PSD operator bound was false is withdrawn. The conference involution has `||T(K)||=3/sqrt(5)>1`, while its PSD projection `(I+K)/2` has half that matching-operator norm. The counterexample concerned the differently scaled involution bound.

The paper does **not** claim the cube-normalized `1/216` optimum, a complete equality classification for the all-orders critical theorem when `r>=3`, best constants at every fixed `N>4r-2`, or quantitative stability around all three PSD equality families.

## Before external submission

The main mathematical contribution now has complete analytic proofs, code, sharp examples, comparison theorems, and reproducibility tests. It still requires independent proof review and a specialist priority review. Human authors, affiliations, author order, corresponding-author details, assistance disclosure, funding/conflict statements, and the target venue must be supplied by the owners before submission. The pull request remains a draft and is not merged automatically.
