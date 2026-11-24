#!/usr/bin/env python3
"""tools/helmholtz/ch4-8a-tables.py: the report's tables of Helmholtz CH4.8a
(docs/helmholtz/tasks/CH4.8a.md) from the outputs of tools/helmholtz/ch4-8a-runs.sh, beside
tools/helmholtz/ch4-4b-tables.py (summary, levels, m2l, worst, pairs), which reads the same files.

    tools/helmholtz/ch4-8a-tables.py <out-dir> rules [<gate-dir>]
    tools/helmholtz/ch4-8a-tables.py <out-dir> cost
    tools/helmholtz/ch4-8a-tables.py <out-dir> raised

  rules   per case and ε the first run of each `remedy-<rule>-…` launch (and the default's, from
          the `map-`, `survey-` or `all-` launch of the case, else the H4 gate's `h4gate-` launch
          in <gate-dir>): the degrees by level, §9.1's ℓ2 of φ and ∇φ over ε, the worst box of φ
          and ∇φ over ε; a cell in bold where §9.1's ℓ2 exceeds ε or the worst box 100 ε
  cost    per case and ε the build and the evaluation of each `cost-<rule>-…` launch (one
          setting per launch), and their sum over the default's
  raised  the `all-` launches' runs: every level with V pairs raised at once by 1 and 2

Each <name>.md is one `helmholtz_boxes` launch. Python 3 standard library only.
"""

import re
import sys
from pathlib import Path

CASE = re.compile(r"^(?P<stage>[a-z0-9]+)-(?P<rule>[a-z0-9]+-)?(?P<case>w[0-9](?::[0-9]+)?-kd[0-9]+-constant)-(?P<eps>1e-[0-9])$")


def runs(path):
    """The rows of the runs table of a launch: (label, p, §9.1 [φ, ∇φ], worst [φ, ∇φ], build, evaluate)."""
    out = []
    for line in path.read_text().splitlines():
        if not line.startswith("| ") or line.startswith("| run ") or line.startswith("| ---"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) != 9:
            continue
        try:
            s91 = [float(v) for v in cells[3].split(",")]
            worst = [float(v) for v in cells[5].split(",")]
            build, evaluate = [float(v) for v in cells[8].split(",")]
        except ValueError:
            continue
        out.append((cells[0], cells[1], s91, worst, build, evaluate))
    return out


def launches(out, stage):
    found = {}
    for path in sorted(Path(out).glob(f"{stage}-*.md")):
        m = CASE.match(path.stem)
        if not m:
            continue
        rule = (m["rule"] or "default-").rstrip("-")
        rows = runs(path)
        if rows:
            found[(m["case"], m["eps"], rule)] = rows
    return found


def bold(text, over):
    return f"**{text}**" if over else text


def cells(eps, row):
    e = float(eps)
    s91 = ", ".join(bold(f"{v / e:.2f}", v > e) for v in row[2])
    worst = ", ".join(bold(f"{v:.1f}", v > 100) for v in row[3])
    return row[1], s91, worst


def rules(out, gate=None):
    remedy = launches(out, "remedy")
    base = {}
    for stage in ["map", "survey", "all"]:
        for (case, eps, _), rows in launches(out, stage).items():
            base.setdefault((case, eps), rows[0])
    if gate:
        for (case, eps, _), rows in launches(gate, "h4gate").items():
            base.setdefault((case, eps), rows[0])
    keys = sorted({(c, e) for (c, e, _) in remedy} | set(base), key=lambda k: (k[0], -float(k[1])))
    names = sorted({r for (_, _, r) in remedy})
    print("| case | ε | rule | p by level | §9.1 φ, ∇φ ÷ ε | worst box φ, ∇φ ÷ ε |")
    print("| --- | --- | --- | --- | --- | --- |")
    for case, eps in keys:
        if (case, eps) in base:
            print(f"| {case} | {eps} | default | " + " | ".join(cells(eps, base[(case, eps)])) + " |")
        for name in names:
            if (case, eps, name) in remedy:
                row = remedy[(case, eps, name)][0]
                print(f"| {case} | {eps} | {name} | " + " | ".join(cells(eps, row)) + " |")


def cost(out):
    found = launches(out, "cost")
    keys = sorted({(c, e) for (c, e, _) in found}, key=lambda k: (k[0], -float(k[1])))
    print("| case | ε | rule | p by level | build s | evaluate s | build + evaluate ÷ default |")
    print("| --- | --- | --- | --- | ---: | ---: | ---: |")
    for case, eps in keys:
        names = ["default"] + sorted(r for (c, e, r) in found if (c, e) == (case, eps) and r != "default")
        ref = found.get((case, eps, "default"))
        for name in names:
            if (case, eps, name) not in found:
                continue
            row = found[(case, eps, name)][0]
            ratio = (row[4] + row[5]) / (ref[0][4] + ref[0][5]) if ref else float("nan")
            print(f"| {case} | {eps} | {name} | {row[1]} | {row[4]:.2f} | {row[5]:.2f} | {ratio:.3f} |")


def raised(out):
    print("| case | ε | run | p by level | §9.1 φ, ∇φ ÷ ε | worst box φ, ∇φ ÷ ε |")
    print("| --- | --- | --- | --- | --- | --- |")
    for (case, eps, _), rows in sorted(launches(out, "all").items(), key=lambda kv: (kv[0][0], -float(kv[0][1]))):
        for row in rows:
            label = "default" if row[0].startswith("ε =") else row[0].split(", the setting")[0]
            print(f"| {case} | {eps} | {label} | " + " | ".join(cells(eps, row)) + " |")


def main():
    if len(sys.argv) < 3:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    out, what = sys.argv[1], sys.argv[2]
    if what == "rules":
        rules(out, sys.argv[3] if len(sys.argv) > 3 else None)
    elif what == "cost":
        cost(out)
    elif what == "raised":
        raised(out)
    else:
        print(__doc__, file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
