"""The technical gate must not mistake passing tests for author approval."""
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("nfc_submission_check", ROOT / "tools/submission_check.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


@pytest.fixture
def candidate(tmp_path):
    for name in MODULE.REQUIRED:
        target = tmp_path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, target)
    path = tmp_path / "submission/review_record.json"
    record = json.loads(path.read_text())
    record["all_authors_approved"] = False
    path.write_text(json.dumps(record))
    return tmp_path


def test_current_candidate_passes_technical_not_submission(candidate):
    result = MODULE.check(candidate)
    assert result["technical_checks_passed"], result["technical_errors"]
    assert not result["recorded_submission_checks_passed"]
    assert result["submission_blockers"]


@pytest.mark.parametrize("defect", ["missing_file", "missing_label", "missing_citation", "long_highlight", "invalid_json"])
def test_rejects_mechanical_defects(candidate, defect):
    if defect == "missing_file":
        (candidate / "docs/SUBMISSION_AUDIT.md").unlink()
    elif defect in ("missing_label", "missing_citation"):
        with (candidate / "paper/matching_dilation.tex").open("a") as stream:
            stream.write(r"\ref{nonexistent-label}" if defect == "missing_label" else r"\cite{nonexistent-source}")
    elif defect == "long_highlight":
        (candidate / "submission/highlights.txt").write_text("x"*86+"\na\nb\n")
    else:
        (candidate / "submission/review_record.json").write_text("not json")
    assert not MODULE.check(candidate)["technical_checks_passed"]


def test_cli_modes_are_distinct(candidate):
    script = ROOT / "tools/submission_check.py"
    ordinary = subprocess.run([sys.executable, str(script), "--root", str(candidate)], capture_output=True, text=True)
    strict = subprocess.run([sys.executable, str(script), "--root", str(candidate), "--submission"], capture_output=True, text=True)
    assert ordinary.returncode == 0
    assert strict.returncode == 1
    assert json.loads(strict.stdout)["submission_blockers"]


def test_strings_are_not_boolean_approvals(candidate):
    path = candidate / "submission/review_record.json"
    record = json.loads(path.read_text())
    record["all_authors_approved"] = "true"
    path.write_text(json.dumps(record))
    assert "Obtain all_authors_approved" in MODULE.check(candidate)["submission_blockers"]
