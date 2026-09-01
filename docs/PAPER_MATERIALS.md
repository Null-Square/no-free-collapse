# Paper preparation materials

This file freezes the journal-independent material for the first manuscript. The target venue and LaTeX template are intentionally **not** chosen here; that should be a separate decision after the repository is reviewer-ready.

## 1. Recommended paper scope

The first paper should be complete **without** the unresolved full PSD-contraction theorem.

Recommended scope:

> A resource theory of interaction order for Born-style reasoning, with exact degree and conditioning limits, exact solvable capacity classes, a low-conditioning hafnian reduction, and a complete six-variable projection-gradient theorem.

The PSD-contraction extension belongs in the final technical section as an exact reduction plus open frontier, strengthened by the proved mixed `(1,5,5,5)` coefficient.

### Paper v1 theorem cutoff

Include as completed results:

- interaction-degree ceiling and tightness;
- normalization loophole;
- conditioning/Chebyshev leakage theorem;
- exact Gram/latent-rank characterization;
- exact symmetric class;
- exact matched-pair class and symmetry breaking;
- universal hafnian bound and exact `m=4` theorem;
- six-variable first/second-order stability;
- rank-one six-variable theorem;
- minimum-trace completion and range-only thin-shell reduction;
- complete projection-gradient theorem, all ranks `0,...,6`;
- exact spectral-chain homogenization;
- mixed nested `(1,5,5,5)` / `(1,1,1,5)` theorem.

Treat as open:

- full PSD-contraction gradient inequality;
- final sharp global six-variable PSD hafnian inequality;
- general even-`m` sharp hafnian optimum;
- sharp global rank-two constant `33/160`.

## 2. Working title options

Do not finalize until journal selection, but the strongest current options are:

1. **No Free Collapse: Interaction-Order Limits of Quantum-Inspired Reasoning**
2. **No Free Collapse: Preparation, Normalization, and Interaction Capacity under Born Readout**
3. **Interaction-Order Limits of Born-Style Reasoning**

Option 1 best captures the conceptual contribution and matches the repository identity.

## 3. One-sentence contribution

> Born-style collapse can expose interference but cannot create arbitrary higher-order interaction for free: fixed-order preparation imposes an exact degree ceiling, nonlinear normalization is a distinct computational resource whose leakage is controlled by conditioning, and the resulting low-conditioning extremal problem admits exact capacity theorems including a complete six-variable projection-gradient bound.

## 4. Provisional abstract ingredients

A final abstract should contain exactly these moves:

1. **Problem:** quantum-inspired reasoning models often attribute expressive power to interference/collapse without isolating what interaction structure was already prepared.
2. **Degree theorem:** fixed-norm order-`r` preparation produces Born interaction degree at most `2r`, sharply.
3. **Normalization:** input-dependent normalization can generate arbitrary exact degree, so it must be charged as a resource.
4. **Conditioning:** bounded norm condition number gives exponential control of high-order Walsh coefficients.
5. **Capacity geometry:** derive exact Gram representation and solve symmetric and matched-pair classes, demonstrating symmetry breaking.
6. **Hafnian regime:** identify the first full-order low-conditioning term with a hafnian and solve the four-variable extremal problem.
7. **Six-variable theorem:** prove `q2(P)<=q1(P)/4` for every real `6 x 6` orthogonal projection.
8. **Frontier:** reduce arbitrary PSD contractions to a nested spectral kernel and prove the first genuinely mixed extreme-rank coefficient.
9. **Interpretation:** collapse exposes prepared structure; higher-order interaction must be paid for in preparation or normalization/conditioning.

Avoid putting unresolved PSD-contraction claims in the abstract.

## 5. Recommended manuscript structure

### 1. Introduction

- motivation: interaction order versus interference/collapse;
- conceptual statement of no free collapse;
- contributions table;
- relation to quantum-inspired reasoning, Boolean Fourier analysis, rational/postselected computation, and hafnian/Gaussian structures;
- theorem roadmap figure.

### 2. Model and interaction order

- Boolean inputs and Walsh basis;
- order-`r` polynomial amplitudes;
- Born score and normalized probability;
- preparation order and interaction degree.

### 3. Fixed-norm degree ceiling

Main theorem:

\[
\deg p\le 2r.
\]

Then give the tight parity construction.

### 4. Normalization is a computational resource

- explicit high-degree normalization witness;
- explain why normalization cannot be hidden inside “collapse”;
- define condition number as a quantitative resource.

### 5. Conditioning-controlled leakage

State and prove the Chebyshev/Helstrom bound

\[
|\widehat p(S)|\le \frac{\rho^m}{1+\rho^{2m}},
\qquad
\rho=\frac{\sqrt\kappa-1}{\sqrt\kappa+1}.
\]

This is the main quantitative resource theorem.

### 6. Gram geometry and optimal readout

- exact `p=z^*Az/(z^*Qz)` representation;
- `0<=A<=Q`;
- minimum latent dimension `rank(Q)`;
- optimal linear readout by eigendecomposition.

### 7. Exact capacity classes and symmetry breaking

- solve permutation/global-sign symmetric class;
- construct matched-pair improvement;
- solve full matched-pair class;
- emphasize that imposing symmetry can hide capacity.

### 8. Low-conditioning expansion and hafnians

- derive full-order hafnian coefficient;
- universal minimax/Chebyshev bound;
- exact `m=4` global optimum.

### 9. Six-variable extremal structure

Organize this as a funnel rather than chronology:

1. local first/second-order pair stability;
2. rank-one global theorem;
3. minimum-trace completion reduction;
4. range-only theorem and thin shell;
5. gradient formulation.

### 10. Complete projection-gradient theorem

State prominently:

\[
\boxed{q_2(P)\le\frac14q_1(P)}
\]

for every real `6 x 6` orthogonal projection.

Recommended proof organization:

- ranks `1/5`;
- rank `3`;
- ranks `2/4`;
- synthesis theorem.

The main text should explain each geometric mechanism; long exact polynomial certificates can move to appendices.

### 11. Beyond projections: spectral contraction frontier

- normalized spectral chain;
- symmetric homogenized four-slot kernel;
- exact reconstruction of the contraction defect;
- emphasize that kernel is not matrix-four-linear;
- prove mixed `(1,5,5,5)` theorem;
- state full PSD-contraction inequality as a conjecture/open problem.

### 12. Discussion

- “collapse exposes rather than manufactures” interpretation;
- preparation versus normalization resource accounting;
- relevance and limits of the quantum analogy;
- open problems.

## 6. Main-text versus appendix split

### Main text

Keep the conceptual and strongest geometric arguments visible:

- degree ceiling + tightness;
- normalization loophole;
- conditioning theorem;
- Gram characterization;
- symmetric and matched-pair exact optima;
- hafnian identification + `m=4` theorem;
- six-variable completion/range funnel;
- projection theorem with proof architecture;
- spectral reduction and mixed `(1,5,5,5)` statement.

### Appendices / supplementary proofs

Move lengthy finite certificates and algebra here:

- detailed Chebyshev polynomial manipulations if space-constrained;
- second-order six-variable tangent expansions;
- elliptope primal/dual witness calculations;
- exact rank-three defect algebra;
- rank-two Pluecker tensor derivation;
- rank-two balanced/high-pair/middle-strip polynomial details;
- rational Bernstein coefficient tables;
- spectral homogenization expansion;
- heavy-atom certificate for `(1,5,5,5)`.

## 7. Proposed result tables

### Table 1 — Resource limits and exact classes

| Layer | Result | Type |
| --- | --- | --- |
| Preparation | `deg p <= 2r` at fixed norm | sharp theorem |
| Normalization | arbitrary exact degree possible | construction |
| Conditioning | exponential high-order leakage bound | quantitative theorem |
| Geometry | Gram representation; rank-minimal latent realization | characterization |
| Symmetric class | exact optimum | exact optimization |
| Matched pairs | exact stronger optimum | symmetry-breaking optimization |

### Table 2 — Hafnian / six-variable hierarchy

| Problem | Status |
| --- | --- |
| universal low-conditioning hafnian bound | proved |
| `m=4` global optimum | solved exactly |
| six-variable rank-one PSD stratum | proved |
| range-only inequality | proved |
| minimum-trace completion reduction | exact |
| all orthogonal projection ranks | proved |
| mixed nested `(1,5,5,5)` | proved |
| arbitrary PSD contraction | open |

### Table 3 — Projection rank proof mechanisms

| Rank | Mechanism |
| --- | --- |
| `0,6` | trivial |
| `1,5` | elementary symmetric polynomials + complementation |
| `3` | involution / exact defect / capacity dual |
| `2,4` | Pluecker geometry / diagonal dual / threshold graph + complementation |

## 8. Figure inventory

Canonical editable figures are in [`FIGURES.md`](FIGURES.md).

Recommended manuscript figures:

1. **Resource chain** — introduction.
2. **Theorem dependency map** — end of introduction or supplement.
3. **Projection-rank closure** — projection theorem section.
4. **Six-variable extremal funnel** — six-variable section.
5. **Evidence hierarchy** — probably repository/supplement only, not necessarily paper.

Potential data/plot figures to generate after journal selection:

- conditioning bound versus exact/symmetric/matched-pair capacities for representative `kappa`;
- symmetric versus matched-pair capacity gap as `n` or `kappa` varies;
- rank-two ratio landscape highlighting the `2+4` block geometry (clearly labeled diagnostic if sharp global constant remains open).

Do not include exploratory plots merely because they are visually interesting; every figure should support a theorem or conceptual claim.

## 9. Notation table for manuscript

| Symbol | Meaning |
| --- | --- |
| `x_i` | Boolean input in `{-1,+1}` |
| `x_S` | Walsh monomial `prod_{i in S} x_i` |
| `r` | preparation order |
| `\tilde psi(x)` | unnormalized latent amplitude |
| `p(x)` | normalized Born output |
| `Q` | preparation Gram / norm matrix |
| `A` | effect/readout Gram matrix, `0<=A<=Q` |
| `kappa` | condition number of squared norm |
| `rho` | `(sqrt(kappa)-1)/(sqrt(kappa)+1)` |
| `C` | zero-diagonal six-variable interaction matrix |
| `tau(C)` | minimum-trace PSD completion cost |
| `a(C),s(C)` | negative/positive Boolean range magnitudes |
| `q1(A)` | off-diagonal edge-square energy |
| `q2(A)` | complementary four-variable hafnian-square energy |
| `P` | orthogonal projection |
| `D` / `mathcal D` | homogenized nested spectral defect kernel |

Avoid recycling `A` for unrelated objects within the same section; the final manuscript should normalize notation globally.

## 10. Novelty boundary to defend carefully

The paper should not claim novelty for background machinery itself:

- Boolean Fourier/Walsh analysis;
- Chebyshev minimax approximation;
- Helstrom discrimination;
- hafnians and Gaussian moment identities;
- generic SDP duality / elliptope methods;
- zeon algebras;
- Schur-Horn projection diagonal theory;
- symmetric-polynomial half-degree principles.

The novelty claim should be attached to the **specific resource theory, reductions, capacity theorems, symmetry-breaking result, six-variable completion/range funnel, and projection-gradient theorem**.

Before submission, perform a fresh literature search around every theorem-level novelty claim and build the related-work section from primary sources.

## 11. Reviewer claims checklist

Before manuscript submission, every theorem in the main text should satisfy:

- [ ] theorem statement matches [`RESULTS.md`](RESULTS.md);
- [ ] proof note or manuscript proof exists;
- [ ] supporting identities/certificates have tests;
- [ ] equality/sharpness statement is proved, not inferred from numerics;
- [ ] prior-art boundary has been checked;
- [ ] open strengthenings are labeled as conjectures;
- [ ] notation agrees across paper and repository;
- [ ] exact repository commit/release is cited.

## 12. Material to prepare after journal selection

The next task, after this repository cleanup, should be:

1. compare appropriate journals/conferences for scope, audience, length, and open-access constraints;
2. choose one target venue;
3. download/use its current official LaTeX template;
4. perform a fresh related-work/novelty search;
5. write the full manuscript directly in the target template;
6. export the Mermaid figures to venue-compatible vector graphics;
7. create a tagged repository release matching the submitted manuscript.

No venue preference is assumed in this file.
