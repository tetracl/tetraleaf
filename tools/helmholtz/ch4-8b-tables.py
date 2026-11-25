#!/usr/bin/env python3
"""tools/helmholtz/ch4-8b-tables.py: the report's tables of Helmholtz CH4.8b
(docs/helmholtz/tasks/CH4.8b.md) from the outputs of tools/helmholtz/ch4-8b-runs.sh.

    tools/helmholtz/ch4-8b-tables.py <out-dir> const
    tools/helmholtz/ch4-8b-tables.py <out-dir> random [<ch4-4a-dir>]
    tools/helmholtz/ch4-8b-tables.py <out-dir> cost
    tools/helmholtz/ch4-8b-tables.py selftest

  const   per case and ε the `const-` launch (the default) beside the `before-` launch (CH4.3c3's
          rule): the degrees by level, §9.1's ℓ2 of φ and ∇φ over ε, the worst box of φ and ∇φ over
          ε (a cell in bold where the ℓ2 exceeds ε or the worst box 100 ε); then the largest of
          each over the cases, per ε, and every miss
  random  per case and ε the `random-` launch (`helmholtz_fmm3d compare`, ours alone): the
          degrees, §9.1's ℓ2 of φ and ∇φ over ε (bold above ε), the build, the evaluation and the
          memory; with <ch4-4a-dir> (CH4.4a's helmholtz-<case>.md) CH4.4a's ℓ2 at the same ε beside
  cost    per case and ε the `cost-…-<i>-<rule>` launches (A B B A, one setting per launch): per
          rule the median of the build, the evaluation (the last of two) and their sum, the
          memory, the sum's spread (max − min) over its launches; the ratios of the default's
          medians over CH4.3c3's; the geometric mean of the sums' ratios per ε and over all
  selftest  the parsers on lines of the formats

Each <name>.md is one launch. Python 3 standard library only.
"""

import math
import re
import statistics
import sys
from pathlib import Path

CONST = re.compile(r"^(?P<stage>const|before)-(?P<case>w[0-9](?:-[a-z]+)?(?:-n[0-9]+)?-kd[0-9]+)-constant-(?P<eps>1e-[0-9])$")
RANDOM = re.compile(r"^random-(?P<case>w[0-9](?:-n[0-9]+)?-kd[0-9]+(?:-planewave)?)$")
COST = re.compile(r"^cost-(?P<case>w[0-9](?:-n[0-9]+)?-kd[0-9]+)-(?P<eps>1e-[0-9])-(?P<i>[0-9])-(?P<rule>default|ch4-3c3)$")


def table_cells(line):
    return [c.strip() for c in line.strip().strip("|").split("|")]


def boxes_run(path):
    """The first run row of a `helmholtz_boxes` launch: (p, §9.1 [φ, ∇φ], worst [φ, ∇φ], build, evaluate)."""
    for line in Path(path).read_text().splitlines():
        if not line.startswith("| ε = "):
            continue
        cells = table_cells(line)
        if len(cells) != 9:
            continue
        try:
            s91 = [float(v) for v in cells[3].split(",")]
            worst = [float(v) for v in cells[5].split(",")]
            build, evaluate = [float(v) for v in cells[8].split(",")]
        except ValueError:
            continue
        return cells[1], s91, worst, build, evaluate
    return None


def compare_rows(path):
    """Our rows of a `helmholtz_fmm3d compare` launch by ε token: (p, build, run, φ, ∇φ, memory GB),
    and the "Build" table's evaluation and first evaluation by ε."""
    rows, build = {}, {}
    in_build = False
    for line in Path(path).read_text().splitlines():
        if line.startswith("## Build"):
            in_build = True
            continue
        if not line.startswith("| ") or line.startswith("| ---"):
            continue
        cells = table_cells(line)
        if not in_build and len(cells) == 17 and cells[4].startswith("ours"):
            memory = cells[16].split(" ")[0]
            rows[cells[3]] = (cells[6], float(cells[7]), float(cells[8]), float(cells[9]),
                              float(cells[10]), float(memory) if memory not in ("–", "-") else None)
        elif in_build and cells[0].startswith("1e-"):
            build[cells[0]] = (float(cells[1]), float(cells[-2]), float(cells[-1]))
    return rows, build


def bold(text, over):
    return f"**{text}**" if over else text


def ratios(values, eps, bound):
    return ", ".join(bold(f"{v / eps:.2f}" if bound == 1 else f"{v:.1f}", (v / eps if bound == 1 else v) > bound)
                     for v in values)




def case_key(case):
    m = re.match(r"w([0-9])(?:-([a-z]+))?(?:-n([0-9]+))?-kd([0-9]+)(-planewave)?", case)
    return (int(m[1]), m[2] or "", int(m[3] or 0), int(m[4]), m[5] or "")


def const(out):
    found = {}
    for path in Path(out).glob("*.md"):
        m = CONST.match(path.stem)
        if m:
            run = boxes_run(path)
            if run:
                found[(m["case"], m["eps"], m["stage"])] = run
    keys = sorted({(c, e) for c, e, _ in found}, key=lambda k: (case_key(k[0]), -int(k[1][-1])))
    print("| case | ε | p by level (CH4.3c3 → default) | §9.1 φ, ∇φ ÷ ε (CH4.3c3 → default) | "
          "worst box φ, ∇φ ÷ ε (CH4.3c3 → default) |")
    print("| --- | --- | --- | --- | --- |")
    worst = {}
    misses = []
    for case, eps in keys:
        e = float(eps)
        after = found.get((case, eps, "const"))
        before = found.get((case, eps, "before"))
        p = " → ".join(r[0] if r else "–" for r in (before, after))
        s91 = " → ".join(ratios(r[1], e, 1) if r else "–" for r in (before, after))
        box = " → ".join(ratios([v * 1.0 for v in r[2]], 1.0, 100) if r else "–" for r in (before, after))
        print(f"| {case}-constant | {eps} | {p} | {s91} | {box} |")
        if after:
            w = worst.setdefault(eps, [0.0, 0.0, "", ""])
            l2 = max(after[1]) / e
            if l2 > w[0]:
                w[0], w[2] = l2, case
            if max(after[2]) > w[1]:
                w[1], w[3] = max(after[2]), case
            if l2 > 1 or max(after[2]) > 100:
                misses.append(f"{case} at {eps}: ℓ2 {l2:.2f} ε, worst box {max(after[2]):.1f} ε")
    print()
    print("| ε | largest §9.1 ℓ2 ÷ ε (case) | largest worst box ÷ ε (case) |")
    print("| --- | --- | --- |")
    for eps in sorted(worst, key=lambda x: -int(x[-1])):
        w = worst[eps]
        print(f"| {eps} | {w[0]:.2f} ({w[2]}) | {w[1]:.1f} ({w[3]}) |")
    print()
    print("Misses under the default:", "; ".join(misses) if misses else "none.")


def random(out, ch4_4a=None):
    found = {}
    for path in Path(out).glob("random-*.md"):
        m = RANDOM.match(path.stem)
        if m:
            found[m["case"]] = compare_rows(path)[0]
    print("| case | ε | p by level | §9.1 φ, ∇φ ÷ ε | CH4.4a's | build s | evaluate s | memory GB |")
    print("| --- | --- | --- | --- | --- | ---: | ---: | ---: |")
    largest = (0.0, "")
    for case in sorted(found, key=case_key):
        before = {}
        if ch4_4a:
            path = Path(ch4_4a) / f"helmholtz-{case}.md"
            if path.exists():
                before = compare_rows(path)[0]
        for eps, (p, build, run, phi, grad, memory) in sorted(found[case].items(), key=lambda x: -int(x[0][-1])):
            e = float(eps)
            earlier = before.get(eps)
            old = f"{earlier[3] / e:.2f}, {earlier[4] / e:.2f}" if earlier else "–"
            mem = f"{memory:.2f}" if memory is not None else "–"
            print(f"| {case} | {eps} | {p} | {ratios([phi, grad], e, 1)} | {old} | {build:.3f} | {run:.3f} | {mem} |")
            if max(phi, grad) / e > largest[0]:
                largest = (max(phi, grad) / e, f"{case} at {eps}")
    print()
    print(f"Largest §9.1 ℓ2 ÷ ε: {largest[0]:.2f} ({largest[1]}).")


def cost(out):
    found = {}
    for path in Path(out).glob("cost-*.md"):
        m = COST.match(path.stem)
        if not m:
            continue
        rows, build = compare_rows(path)
        if m["eps"] not in rows:
            continue
        p, b, run, _, _, memory = rows[m["eps"]]
        found.setdefault((m["case"], m["eps"]), {}).setdefault(m["rule"], []).append((p, b, run, memory))
    print("| case | ε | p by level (CH4.3c3 / default) | build s | evaluate s | build + evaluate s (spread) | "
          "memory GB | default ÷ CH4.3c3: build, evaluate, sum, memory |")
    print("| --- | --- | --- | --- | --- | --- | --- | --- |")
    by_eps = {}
    for case, eps in sorted(found, key=lambda k: (case_key(k[0]), -int(k[1][-1]))):
        rules = found[(case, eps)]
        if "default" not in rules or "ch4-3c3" not in rules:
            continue
        med = {}
        for rule, runs in rules.items():
            sums = [b + r for _, b, r, _ in runs]
            mems = [m for *_, m in runs if m is not None]
            med[rule] = (runs[0][0], statistics.median(b for _, b, _, _ in runs),
                         statistics.median(r for _, _, r, _ in runs), statistics.median(sums),
                         max(sums) - min(sums), statistics.median(mems) if mems else None)
        a, b = med["ch4-3c3"], med["default"]
        mem = " / ".join(f"{m[5]:.2f}" if m[5] is not None else "–" for m in (a, b))
        ratio = [b[1] / a[1], b[2] / a[2], b[3] / a[3]]
        mratio = f"{b[5] / a[5]:.3f}" if a[5] and b[5] else "–"
        by_eps.setdefault(eps, []).append(ratio[2])
        print(f"| {case} | {eps} | {a[0]} / {b[0]} | {a[1]:.3f} / {b[1]:.3f} | {a[2]:.3f} / {b[2]:.3f} | "
              f"{a[3]:.3f} ({a[4]:.3f}) / {b[3]:.3f} ({b[4]:.3f}) | {mem} | "
              f"{ratio[0]:.3f}, {ratio[1]:.3f}, {ratio[2]:.3f}, {mratio} |")
    print()
    every = [r for rs in by_eps.values() for r in rs]
    for eps in sorted(by_eps, key=lambda x: -int(x[-1])):
        rs = by_eps[eps]
        print(f"- {eps}: build + evaluate {min(rs):.3f}–{max(rs):.3f}×, geometric mean "
              f"{math.exp(sum(map(math.log, rs)) / len(rs)):.3f}")
    if every:
        print(f"- all: {min(every):.3f}–{max(every):.3f}×, geometric mean "
              f"{math.exp(sum(map(math.log, every)) / len(every)):.3f}")


def selftest():
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        boxes = Path(d) / "const-w2-kd30-constant-1e-6.md"
        boxes.write_text(
            "| run | p by level | leaves | §9.1 φ, ∇φ | all φ, ∇φ | worst box / ε | worst / ℓ2 | boxes > ε φ, ∇φ | build, evaluate s |\n"
            "| ε = 1e-6, defaults | 0×2, 17, 15×2 | 232 | 1.65e-8, 1.38e-7 | 1.74e-8, 1.70e-7 | 0.46, 7.14 | 27, 42 | 0, 16 of 1080 | 0.04, 0.05 |\n")
        assert boxes_run(boxes) == ("0×2, 17, 15×2", [1.65e-8, 1.38e-7], [0.46, 7.14], 0.04, 0.05)
        cmp_ = Path(d) / "cost-w1-kd30-1e-3-1-default.md"
        cmp_.write_text(
            "| Workload | N | kD | ε | code | threads | p by level | build s | run s | φ | ∇φ | all φ | all ∇φ | worst box φ / ε, ∇φ / ε | worst box / all (φ / ∇φ) of boxes | checksum | memory GB |\n"
            "| W1 | 7162 | 30 | 1e-3 | ours (avx512) | 1 × 40 | 0×2, 12 | 0.015 | 0.010 | 1.30e-5 | 2.26e-5 | – | – | – | – | 1673d481b3cbf8ec | 0.02 (rank ≤ 0.02) |\n"
            "\n## Build\n\n| ε | build | settings | evaluation | first evaluation |\n"
            "| 1e-3 | 0.0153 | 0.0030 | 0.0104 | 0.0110 |\n")
        rows, build = compare_rows(cmp_)
        assert rows == {"1e-3": ("0×2, 12", 0.015, 0.010, 1.30e-5, 2.26e-5, 0.02)}, rows
        assert build == {"1e-3": (0.0153, 0.0104, 0.0110)}, build
        assert CONST.match("before-w3-n1000000-kd500-constant-1e-8")
        assert CONST.match("const-w4-plummer-n100000-kd30-constant-1e-6")  # the H4 gate's W4 run
        assert RANDOM.match("random-w1-kd1000-planewave")
        assert COST.match("cost-w3-n1000000-kd500-1e-6-4-ch4-3c3")
    print("selftest: ok")


if __name__ == "__main__":
    if sys.argv[1:] == ["selftest"]:
        selftest()
    elif len(sys.argv) >= 3 and sys.argv[2] == "const":
        const(sys.argv[1])
    elif len(sys.argv) >= 3 and sys.argv[2] == "random":
        random(sys.argv[1], sys.argv[3] if len(sys.argv) > 3 else None)
    elif len(sys.argv) >= 3 and sys.argv[2] == "cost":
        cost(sys.argv[1])
    else:
        sys.exit(__doc__)
