#!/usr/bin/env python3
"""tools/helmholtz/ch4-3c2-tables.py: the tables of Helmholtz CH4.3c2's report
(`tetraleaf-validate/results/helmholtz-ch4-3c2.md`) from the outputs of
`tools/helmholtz/ch4-3c2-runs.sh`.

    tools/helmholtz/ch4-3c2-tables.py matrix <dir>        # §5: verify-matrix's outputs
    tools/helmholtz/ch4-3c2-tables.py w1 <file.md> …      # §6: helmholtz_boxes outputs

`matrix` reads `<rule>-<m2l>-<case>.md` (`helmholtz_accuracy`, rule `calibrated` or `ch3-3`,
m2l `rotation` or `defaults`), pairs each run under the default rule with the same run under
CH3.3's, prints the runs whose degrees differ (the degrees, the larger of φ's and ∇φ's
relative ℓ2 error over ε under each, the evaluation seconds) and the largest error over ε per
M2L, with every run above ε. `w1` prints one row per run of `helmholtz_boxes` (the file's
name, ε, M2L, degrees, φ and ∇φ over ε at the samples, the evaluation seconds). Python 3
standard library only; prints Markdown.
"""

import collections
import glob
import os
import sys


def table_rows(lines, header_start):
    """The rows of the first Markdown table whose header starts with `header_start`."""
    for i, line in enumerate(lines):
        if line.startswith(header_start):
            header = [c.strip() for c in line.strip().strip("|").split("|")]
            for row in lines[i + 2 :]:
                if not row.startswith("|"):
                    return
                yield dict(zip(header, [c.strip() for c in row.strip().strip("|").split("|")]))
            return


def matrix(directory):
    runs = {}
    for path in sorted(glob.glob(os.path.join(directory, "*.md"))):
        name = os.path.basename(path)[:-3]
        rule, rest = name.split("-", 1) if name.startswith("calibrated-") else ("ch3-3", name[6:])
        m2l, case = rest.split("-", 1)
        lines = open(path).read().split("\n")
        for c in table_rows(lines, "| Workload"):
            eps = float(c["degrees"].split("=")[1])
            key = (m2l, case, c["Workload"], c["N"], c["kD"], eps)
            runs[(rule,) + key] = (
                max(float(c["φ"]), float(c["∇φ"])) / eps,
                c["p by level"],
                float(c["evaluate s"]),
            )
    worst = collections.defaultdict(lambda: [0.0, 0.0])
    over = []
    print("| M2L | case | N | kD | ε | p, default | p, CH3.3 | error / ε, default | CH3.3 | evaluate s, default / CH3.3 |")
    print("| --- | --- | ---: | --- | --- | --- | --- | ---: | ---: | --- |")
    for key in sorted(k[1:] for k in runs if k[0] == "calibrated"):
        new, old = runs[("calibrated",) + key], runs.get(("ch3-3",) + key)
        worst[key[0]][0] = max(worst[key[0]][0], new[0])
        if old:
            worst[key[0]][1] = max(worst[key[0]][1], old[0])
        if new[0] > 1.0:
            over.append(key)
        if old and new[1] != old[1]:
            m2l, case, _, n, kd, eps = key
            print(
                f"| {m2l} | {case} | {n} | {kd} | {eps:.0e} | {new[1]} | {old[1]} | {new[0]:.3f} | "
                f"{old[0]:.3f} | {new[2]:.2f} / {old[2]:.2f} |"
            )
    print()
    print(f"{sum(1 for k in runs if k[0] == 'calibrated')} runs under the default; the largest error / ε:")
    for m2l, (new, old) in sorted(worst.items()):
        print(f"- {m2l}: {new:.3f} (CH3.3's rule {old:.3f})")
    print(f"- above ε: {over if over else 'none'}")


def w1(paths):
    print("| run | ε | M2L | p by level | φ / ε | ∇φ / ε | evaluate s |")
    print("| --- | --- | --- | --- | ---: | ---: | ---: |")
    for path in paths:
        lines = open(path).read().split("\n")
        for c in table_rows(lines, "| run |"):
            if not c["run"].startswith("ε = "):
                continue
            eps_text, m2l = c["run"][4:].split(", ")
            eps = float(eps_text)
            phi, grad = (float(v) / eps for v in c["§9.1 φ, ∇φ"].split(", "))
            seconds = c["build, evaluate s"].split(", ")[1]
            name = os.path.basename(path)[:-3]
            print(f"| {name} | {eps:.0e} | {m2l} | {c['p by level']} | {phi:.2f} | {grad:.2f} | {seconds} |")


def main():
    if len(sys.argv) < 3 or sys.argv[1] not in ("matrix", "w1"):
        print(__doc__)
        sys.exit(2)
    if sys.argv[1] == "matrix":
        matrix(sys.argv[2])
    else:
        w1(sys.argv[2:])


if __name__ == "__main__":
    main()
