"""The sign-off rules of docs/CLINICAL_REVIEW.md, required before protocol-v1 (ROADMAP M8): a
ticked sign-off must carry its names and date, a filled-in sign-off must be ticked, and an item
may be ticked only with the sign-offs it needs (both for items that need a clinician)."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

DOC = Path(__file__).resolve().parents[1] / "docs" / "CLINICAL_REVIEW.md"
ITEM = re.compile(r"^- \[( |x)\] \*\*(CR\d+)\.")
SIGN = re.compile(r"^  - \[( |x)\] (Literature checked|Clinician reviewed)(.*)$")
FILLED = {
    "Literature checked": re.compile(r"^: [^,]+, \d{4}-\d{2}-\d{2}\.$"),
    "Clinician reviewed": re.compile(r"^: [^,]+, [^,]+, \d{4}-\d{2}-\d{2}\.$"),
}
BLANK = {"Literature checked": " (name, date)", "Clinician reviewed": " (name, role, date)"}


def problems(text: str) -> list[str]:
    out, items, cur = [], [], None
    for line in text.splitlines():
        m = ITEM.match(line)
        if m:
            cur = {"id": m[2], "ticked": m[1] == "x", "clinician": False, "signs": {}}
            items.append(cur)
            continue
        if cur is None:
            continue
        if line.startswith("  - **Needs a clinician"):
            cur["clinician"] = True
        s = SIGN.match(line)
        if s:
            ticked, kind, rest = s[1] == "x", s[2], s[3]
            filled = bool(FILLED[kind].match(rest))
            if rest != BLANK[kind] and not filled:
                out.append(f"{cur['id']} {kind}: malformed sign-off {rest!r}")
            if ticked and not filled:
                out.append(f"{cur['id']} {kind}: ticked without its names and date")
            if filled and not ticked:
                out.append(f"{cur['id']} {kind}: filled in but not ticked")
            cur["signs"][kind] = ticked and filled
    for it in items:
        need = ["Literature checked"] + (["Clinician reviewed"] if it["clinician"] else [])
        missing = [k for k in need if not it["signs"].get(k)]
        if it["ticked"] and missing:
            out.append(f"{it['id']}: ticked without {', '.join(missing)}")
    if not items:
        out.append("no items found")
    return out


def test_clinical_review_sign_offs_are_consistent():
    assert problems(DOC.read_text()) == []


ITEM_TEXT = """- [{t}] **CR1. A claim.**
  - **Needs a clinician:** a clinical judgment.
  - [{l}] Literature checked{lt}
  - [{c}] Clinician reviewed{ct}
"""


@pytest.mark.parametrize(
    ("kw", "bad"),
    [
        (dict(t=" ", l=" ", lt=" (name, date)", c=" ", ct=" (name, role, date)"), False),
        (dict(t=" ", l="x", lt=" (name, date)", c=" ", ct=" (name, role, date)"), True),
        (dict(t=" ", l=" ", lt=": A. Reader, 2026-11-01.", c=" ", ct=" (name, role, date)"), True),
        (dict(t="x", l="x", lt=": A. Reader, 2026-11-01.", c=" ", ct=" (name, role, date)"), True),
        (
            dict(
                t="x",
                l="x",
                lt=": A. Reader, 2026-11-01.",
                c="x",
                ct=": B. Doc, endocrinologist, 2026-11-02.",
            ),
            False,
        ),  # fmt: skip
        (dict(t=" ", l="x", lt=": A. Reader, 2026-11.", c=" ", ct=" (name, role, date)"), True),
    ],
)
def test_sign_off_rules_reject_what_they_should(kw, bad):
    assert bool(problems(ITEM_TEXT.format(**kw))) == bad
