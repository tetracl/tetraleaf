#!/usr/bin/env python3
"""tools/helmholtz/ch4-4b-tables.py: the report's tables of Helmholtz CH4.4b
(docs/helmholtz/tasks/CH4.4b.md) from the outputs of tools/helmholtz/ch4-4b-runs.sh.

    tools/helmholtz/ch4-4b-tables.py <out-dir> summary|levels|m2l|worst [<name regex>]

Each <out-dir>/<name>.md is one `helmholtz_boxes` launch (one case, one ε, the library's
defaults). Tables (Markdown on stdout):
  summary  per launch: the worst box of φ and ∇φ over ε (the refined reference), the ℓ2 over the
           reference's targets, the level with the largest part of Σ|Δ∇φ|² and its share, the
           rest, and of the M2L study the largest change of the worst box and of the ℓ2 error
           (run ÷ first, the farthest from 1) with its run, and the classes beside the threshold
  levels   per launch and level: k a, p, M2L, the share of Σ|e|² φ and ∇φ, the classes by decade
  m2l      per launch and M2L run: the share, the worst box and ℓ2 ÷ first, the part at the worst
           ∇φ target
  worst    per launch: the worst ∇φ boxes with the levels' and the M2L runs' parts
Python 3 standard library only.
"""

import re
import sys
from pathlib import Path


def tables(text):
    """The Markdown tables of a file: a list of (heading paragraph, header cells, rows)."""
    out = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith("|") and i + 1 < len(lines) and lines[i + 1].startswith("| -"):
            head = [c.strip() for c in lines[i].strip("|").split("|")]
            # The paragraph before the table.
            j = i - 1
            while j >= 0 and not lines[j].strip():
                j -= 1
            k = j
            while k >= 0 and lines[k].strip():
                k -= 1
            heading = " ".join(lines[k + 1 : j + 1])
            rows = []
            i += 2
            while i < len(lines) and lines[i].startswith("|"):
                rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
                i += 1
            out.append((heading, head, rows))
        else:
            i += 1
    return out


def find(ts, prefix):
    for heading, head, rows in ts:
        if heading.startswith(prefix):
            return head, rows
    return None, []


def pair(cell):
    a, b = cell.split(",")
    return float(a), float(b)


def launches(out, pattern):
    names = sorted(p.stem for p in Path(out).glob("*.md"))
    rx = re.compile(pattern)
    for n in names:
        if rx.search(n):
            yield n, (Path(out) / f"{n}.md").read_text()


def summary(out, pattern):
    print(
        "| launch | worst box ÷ ε φ, ∇φ | ℓ2 φ, ∇φ | level (share of Σ|Δ∇φ|²) | rest φ, ∇φ "
        "| largest M2L change: worst box ∇φ (run) | ℓ2 ∇φ (run) | classes in (0.1, 1] |"
    )
    print("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for name, text in launches(out, pattern):
        ts = tables(text)
        _, runs = find(ts, "Runs:")
        if not runs:
            print(f"| {name} | (no output) | | | | | | |")
            continue
        first = runs[0]
        worst = first[5]
        l2 = first[4]
        _, shares = find(ts, "Each level's part")
        level, rest = "–", "–"
        if shares:
            body = [r for r in shares if r[0] != "rest"]
            if body:
                top = max(body, key=lambda r: float(r[6]))
                level = f"{top[0]} ({top[2]}, {float(top[6]):.2f})"
            rest = ", ".join(r for r in [shares[-1][5], shares[-1][6]])
        _, m2l = find(ts, "The M2L's parts")
        change_box, change_l2 = "–", "–"
        if m2l:
            box = max(m2l, key=lambda r: abs(pair(r[3])[1] - 1))
            ell = max(m2l, key=lambda r: abs(pair(r[4])[1] - 1))
            change_box = f"{pair(box[3])[1]:.3f} ({box[0]})"
            change_l2 = f"{pair(ell[4])[1]:.3f} ({ell[0]})"
        _, levels = find(ts, "Levels as run")
        beside = [
            f"{r[0]}: {r[5]}" for r in levels if len(r) > 5 and r[5] not in ("", "[]")
        ]
        print(
            f"| {name} | {worst} | {l2} | {level} | {rest} | {change_box} | {change_l2} "
            f"| {'; '.join(beside) or '–'} |"
        )


def levels(out, pattern):
    print("| launch | level | k a | p | M2L | share φ | share ∇φ | classes (0.1, 1] | (0.01, 0.1] | failing |")
    print("| --- | ---: | ---: | ---: | --- | ---: | ---: | --- | ---: | ---: |")
    for name, text in launches(out, pattern):
        ts = tables(text)
        _, lv = find(ts, "Levels as run")
        _, shares = find(ts, "Each level's part")
        share = {r[0]: (r[5], r[6]) for r in shares}
        for r in lv:
            s = share.get(r[0], ("–", "–"))
            print(
                f"| {name} | {r[0]} | {r[1]} | {r[2]} | {r[4]} | {s[0]} | {s[1]} | {r[5]} | {r[6]} "
                f"| {r[8]} |"
            )


def m2l(out, pattern):
    print("| launch | run | share φ | share ∇φ | worst box ÷ first φ, ∇φ | ℓ2 ÷ first φ, ∇φ | at the worst ∇φ target |")
    print("| --- | --- | ---: | ---: | --- | --- | --- |")
    for name, text in launches(out, pattern):
        _, rows = find(tables(text), "The M2L's parts")
        for r in rows:
            print(f"| {name} | " + " | ".join(r) + " |")


def worst(out, pattern):
    for name, text in launches(out, pattern):
        head, rows = find(tables(text), "Worst boxes of ∇φ")
        if not rows:
            continue
        print(f"**{name}**")
        print()
        print("| " + " | ".join(head) + " |")
        print("|" + " --- |" * len(head))
        for r in rows:
            print("| " + " | ".join(r) + " |")
        print()


def main():
    if len(sys.argv) < 3:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    out, what = sys.argv[1], sys.argv[2]
    pattern = sys.argv[3] if len(sys.argv) > 3 else "."
    {"summary": summary, "levels": levels, "m2l": m2l, "worst": worst}[what](out, pattern)


if __name__ == "__main__":
    main()
