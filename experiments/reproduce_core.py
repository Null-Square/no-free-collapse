"""Run the curated No Free Collapse experiment sequence.

This script is for reviewer-facing reproduction of representative outputs. It
is intentionally separate from the pytest suite: theorem status is determined
by analytic proofs and exact certificates, not by experiment output.

Use ``--include-search`` to also run the exploratory asymmetric-search and
pair-breaking scripts.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
EXPERIMENTS = ROOT / "experiments"

CORE = (
    "e1_linear_barrier.py",
    "e2_normalization_loophole.py",
    "e3_tight_parity.py",
    "e4_normalization_explosion.py",
    "e5_conditioned_leakage.py",
    "e6_optimal_collapse.py",
    "e7_chebyshev_vs_mean_field.py",
    "e10_hafnian_bound.py",
    "e11_relaxed_six_variable_witness.py",
    "e12_rank_three_defect_scan.py",
)

EXPLORATORY = (
    "e8_asymmetric_search.py",
    "e9_pair_breaking.py",
)


def run_script(name: str) -> None:
    path = EXPERIMENTS / name
    if not path.is_file():
        raise FileNotFoundError(path)
    print(f"\n{'=' * 72}\n{name}\n{'=' * 72}", flush=True)
    subprocess.run([sys.executable, str(path)], cwd=ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--include-search",
        action="store_true",
        help="also run exploratory asymmetric-search / pair-breaking scripts",
    )
    args = parser.parse_args()

    selected = list(CORE)
    if args.include_search:
        selected.extend(EXPLORATORY)

    print("No Free Collapse — curated reproduction")
    print("Theorem status is defined by docs/RESULTS.md, not by these outputs.")
    for name in selected:
        run_script(name)

    print(f"\nCompleted {len(selected)} experiment scripts successfully.")


if __name__ == "__main__":
    main()
