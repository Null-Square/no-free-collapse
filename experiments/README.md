# Experiments

The experiment scripts reproduce representative calculations and document the discovery path. They are **not** substitutes for analytic proofs. The theorem ledger is [`../docs/RESULTS.md`](../docs/RESULTS.md).

Run the curated sequence with:

```bash
python experiments/reproduce_core.py
```

Add the explicitly exploratory search scripts with:

```bash
python experiments/reproduce_core.py --include-search
```

## Script index

| Script | Role | Notes |
| --- | --- | --- |
| `e1_linear_barrier.py` | deterministic illustration | fixed-norm / low-order interaction barrier |
| `e2_normalization_loophole.py` | deterministic construction | input-dependent normalization creates higher-order interaction |
| `e3_tight_parity.py` | deterministic construction | factor-of-two degree ceiling is attained |
| `e4_normalization_explosion.py` | deterministic illustration | high exact degree produced through normalization |
| `e5_conditioned_leakage.py` | deterministic comparison | conditioning-controlled leakage |
| `e6_optimal_collapse.py` | deterministic calculation | optimal final readout / collapse geometry |
| `e7_chebyshev_vs_mean_field.py` | deterministic comparison | universal Chebyshev bound versus symmetric mean-field class |
| `e8_asymmetric_search.py` | **exploratory numerical search** | helped discover symmetry breaking; output is not a theorem |
| `e9_pair_breaking.py` | **exploratory numerical search** | pair-structure perturbation diagnostics; output is not a theorem |
| `e10_hafnian_bound.py` | deterministic illustration | low-conditioning hafnian bound |
| `e11_relaxed_six_variable_witness.py` | deterministic diagnostic | distinguishes a relaxed non-PSD witness from the PSD-completion problem |
| `e12_rank_three_defect_scan.py` | seeded regression / diagnostic | sanity scan for the rank-three defect identity; the global theorem is analytic |
| `reproduce_core.py` | runner | executes the curated set in paper order |

## Interpretation rule

A script may motivate or illustrate a statement, but a statement is considered proved only if it is marked **PROVED** or **SOLVED EXACTLY** in [`../docs/RESULTS.md`](../docs/RESULTS.md) and has an analytic proof note.

This distinction is especially important for open sharpenings such as the global rank-two `33/160` constant and the full PSD-contraction inequality.
