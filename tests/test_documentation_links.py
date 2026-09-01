from __future__ import annotations

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
LINK_RE = re.compile(r"\[[^\]]+\]\(([^)]+)\)")


def _markdown_files() -> list[Path]:
    files = [ROOT / "README.md", ROOT / "experiments" / "README.md", ROOT / "tests" / "README.md"]
    files.extend(sorted((ROOT / "docs").glob("*.md")))
    return files


def test_local_markdown_links_resolve() -> None:
    missing: list[str] = []
    for source in _markdown_files():
        text = source.read_text(encoding="utf-8")
        for raw_target in LINK_RE.findall(text):
            target = raw_target.strip().split()[0]
            if not target or target.startswith(("http://", "https://", "mailto:", "#")):
                continue
            path_part = target.split("#", 1)[0]
            if not path_part:
                continue
            resolved = (source.parent / path_part).resolve()
            if not resolved.exists():
                missing.append(f"{source.relative_to(ROOT)} -> {target}")

    assert not missing, "broken local Markdown links:\n" + "\n".join(missing)
