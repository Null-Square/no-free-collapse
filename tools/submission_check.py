"""Mechanical preflight; it cannot certify mathematical review or author consent.

Default mode checks the technical source inventory. --submission additionally
requires recorded human approvals and removal of candidate-only notices.
This script never submits, approves, or changes any file.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import re

REQUIRED = (
    "README.md", "docs/RESULTS.md", "docs/REPRODUCIBILITY.md",
    "docs/SUBMISSION_AUDIT.md", "docs/matching_dilation.md",
    "docs/subhafnian_hierarchy.md", "paper/matching_dilation.tex",
    "paper/author_metadata.tex", "submission/README.md",
    "submission/cover_letter.md", "submission/highlights.txt",
    "submission/ai_disclosure.md", "submission/literature_review.md",
    "submission/review_record.json", "src/no_free_collapse/hierarchy_certificate.py",
    "tests/test_submission_audit.py",
)


def _text(value: object) -> bool:
    return isinstance(value, str) and bool(value.strip())


def check(root: Path) -> dict[str, object]:
    """Check local files without running code or trusting metadata as a review."""
    technical: list[str] = []
    blockers: list[str] = []
    for relative in REQUIRED:
        path = root / relative
        if not path.is_file() or path.stat().st_size == 0:
            technical.append(f"Missing or empty file: {relative}")
    main = root / "paper/matching_dilation.tex"
    metadata = root / "paper/author_metadata.tex"
    if main.is_file():
        tex = main.read_text(encoding="utf-8")
        # This preflight supports the supplied single-file bibliography, not arbitrary TeX.
        stripped = re.sub(r"(?<!\\)%[^\n]*", "", tex)
        labels = re.findall(r"\\label\{([^{}]+)\}", stripped)
        duplicates = [key for key, count in Counter(labels).items() if count > 1]
        technical.extend(f"Duplicate LaTeX label: {key}" for key in duplicates)
        refs = re.findall(r"\\(?:eqref|ref|autoref)\{([^{}]+)\}", stripped)
        technical.extend(f"Undefined LaTeX label: {key}" for key in sorted(set(refs)-set(labels)))
        bib = set(re.findall(r"\\bibitem(?:\[[^\]]*\])?\{([^{}]+)\}", stripped))
        cited = {key.strip() for group in re.findall(r"\\cite\w*(?:\[[^\]]*\])?\{([^{}]+)\}", stripped)
                 for key in group.split(",")}
        technical.extend(f"Undefined bibliography key: {key}" for key in sorted(cited-bib))
        pending = ("has not yet received documented human approval", "pending author confirmation")
        if any(phrase in tex for phrase in pending):
            blockers.append("Replace candidate-only approval notices in the manuscript after actual review")
    if metadata.is_file():
        meta_text = metadata.read_text(encoding="utf-8")
        if "pending confirmation" in meta_text.lower() or "remain pending" in meta_text.lower():
            blockers.append("Confirm and replace the visible manuscript author metadata")
    highlights = root / "submission/highlights.txt"
    if highlights.is_file():
        lines = [line.strip() for line in highlights.read_text(encoding="utf-8").splitlines() if line.strip()]
        if not 3 <= len(lines) <= 5:
            technical.append("Highlights must contain 3-5 nonempty lines")
        if any(len(line) > 85 for line in lines):
            technical.append("A highlight exceeds 85 characters")
    record_path = root / "submission/review_record.json"
    record: dict = {}
    if record_path.is_file():
        try:
            loaded = json.loads(record_path.read_text(encoding="utf-8"))
            if not isinstance(loaded, dict):
                raise ValueError("Expected an object")
            record = loaded
        except (ValueError, UnicodeError) as exc:
            technical.append(f"Invalid review record: {exc}")
    authors = record.get("authors")
    if not isinstance(authors, list) or not authors or any(
        not isinstance(author, dict) or not _text(author.get("name"))
        or not _text(author.get("affiliation")) for author in authors
    ):
        blockers.append("Record the human authors, order and affiliations")
    email = record.get("corresponding_author_email")
    if not _text(email) or not re.fullmatch(r"[^\s@]+@[^\s@]+\.[^\s@]+", email):
        blockers.append("Record a valid corresponding-author email")
    for key in ("target_journal", "funding_statement", "competing_interests_statement",
                "contributions_statement", "approved_ai_disclosure"):
        if not _text(record.get(key)):
            blockers.append(f"Confirm {key}")
    for key in ("proof_review", "priority_review"):
        review = record.get(key)
        if not isinstance(review, dict) or review.get("approved") is not True or not all(
            _text(review.get(field)) for field in ("reviewer", "record")
        ):
            blockers.append(f"Record an actual human {key} and sign-off")
    for key in ("all_authors_approved", "exclusive_submission_confirmed", "venue_requirements_verified"):
        if record.get(key) is not True:
            blockers.append(f"Obtain {key}")
    return {
        "technical_checks_passed": not technical,
        "technical_errors": technical,
        "recorded_submission_checks_passed": not technical and not blockers,
        "submission_blockers": blockers,
        "notice": "Mechanical checklist only: values must record genuine human review and approval; no submission is performed.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--submission", action="store_true", help="Fail until all human approval fields are complete")
    args = parser.parse_args()
    result = check(args.root.resolve())
    print(json.dumps(result, indent=2))
    key = "recorded_submission_checks_passed" if args.submission else "technical_checks_passed"
    return 0 if result[key] else 1


if __name__ == "__main__":
    raise SystemExit(main())
