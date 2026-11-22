#!/usr/bin/env python3
"""tools/helmholtz/ch4-3c3-tables.py: the tables of Helmholtz CH4.3c3's report
(`tetraleaf-validate/results/helmholtz-ch4-3c3.md`) from the outputs of
`tools/helmholtz/ch4-3c3-runs.sh scans` (`helmholtz_calibrate`'s `point:` lines).

    tools/helmholtz/ch4-3c3-tables.py scans <file.md> …    # one row per scan: p at ε/2 per ε
    tools/helmholtz/ch4-3c3-tables.py margin <file.md> …   # per x: the margin against the rules
    tools/helmholtz/ch4-3c3-tables.py matrix <dir>         # `verify-matrix`'s outputs

`scans` prints one row per (workload, density, kD, level) with the smallest degree meeting ε/2
at each ε of 10⁻³ … 10⁻⁸ (a launch split by ε, kD = 1000's levels 2 and 3, gives one row from
both files). `margin` takes, per x (to 0.1), the largest degree over the scans of a density other
than the random charges, and prints per ε: that degree p_c, the random charges' measured degree
at the same x (CH4.3c2's report §2 at x = 27–216; CH3.3's at x = 13.57, W4's clusters), CH4.3c2's
fit (d = 3) and CH4.3c2's rule (d = 5, the default before CH4.3c3) at the scan's x, p_c − fit (the
margin needed) and rule − p_c (the rule's spare), with a last line per x of the largest margin needed and the
smallest spare. `matrix` reads `<rule>-<m2l>-<case>.md` (`helmholtz_accuracy`, rule `calibrated`
or `ch4-3c2`, m2l `rotation` or `defaults`), pairs each run under the default rule with the same
run under CH4.3c2's, prints the runs whose degrees differ (the degrees, the larger of φ's and ∇φ's
relative ℓ2 error over ε under each, the evaluation seconds), the largest error over ε per M2L and
rule, and every run above ε. Python 3 standard library only; prints Markdown.
"""

import collections
import glob
import math
import os
import sys

EPSILONS = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]

# The random charges' measured degrees at ε/2 (the largest over the workloads at one x):
# CH4.3c2's report §2, and CH3.3's level scan of W4's clusters at x = 13.57 (CH3.3's report §3).
RANDOM = {
    13.6: [13, 16, 18, 20, 22, 24],
    27.1: [32, 34, 37, 39, 41, 44],
    54.1: [58, 62, 65, 68, 71, 74],
    108.3: [112, 117, 121, 126, 129, 132],
    216.5: [218, 225, 230, 235, 240, 244],
}


def up(v):
    return max(0, math.ceil(v - 1e-9))


def degree(x, eps, d):
    """`DegreeRule::CH4_3C2` with offset d (5 CH4.3c2's rule, 3 its fit)."""
    big_l = math.log10(1 / eps)
    ln23 = math.log(1 / eps) ** (2 / 3)
    plateau = up((0.30 * big_l + 0.25) * big_l + 2.60)
    excess = max(0.0, 2.30 * ln23 - 0.75 * big_l - 5.73)
    branch = min(x + excess * x ** (1 / 3) + d, x + 0.80 * ln23 * x ** (1 / 3))
    return max(plateau, up(branch))


def points(files):
    """(workload, density, kd, level, x, eps, p) of every `point:` line of the files."""
    out = []
    for path in files:
        with open(path) as f:
            for line in f:
                if not line.startswith("point: "):
                    continue
                kind, workload, kd, level, x, eps, half, _, _ = line.split()[1:]
                if kind != "level":
                    continue
                name, _, density = workload.partition("/")
                out.append(
                    (name, density or "random", kd, int(level), float(x), float(eps), int(half))
                )
    return out


def scans(files):
    rows = collections.OrderedDict()
    for name, density, kd, level, x, eps, p in points(files):
        key = (name, density, kd, level)
        rows.setdefault(key, [x, {}])[1][eps] = p
    print("| workload | density | kD | level | x | " + " | ".join(f"{e:.0e}" for e in EPSILONS) + " |")
    print("| --- | --- | ---: | ---: | ---: |" + " ---: |" * len(EPSILONS))
    for (name, density, kd, level), (x, by_eps) in sorted(
        rows.items(), key=lambda r: (float(r[0][2]), r[0][3])
    ):
        cells = [str(by_eps.get(e, "–")) for e in EPSILONS]
        print(f"| {name} | {density} | {kd} | {level} | {x:.3f} | " + " | ".join(cells) + " |")


def margin(files):
    by_x = collections.defaultdict(dict)  # x (0.1) -> eps -> (p, x, sources)
    for name, density, kd, level, x, eps, p in points(files):
        if density == "random":
            continue
        key = round(x, 1)
        best = by_x[key].get(eps)
        source = f"{kd}/{level}"
        if best is None or p > best[0]:
            by_x[key][eps] = (p, x, [source])
        elif p == best[0]:
            best[2].append(source)
    print(
        "| x | ε | p_c | kD/level | random | fit (d = 3) | CH4.3c2 (d = 5) | p_c − fit | "
        "CH4.3c2 − p_c |"
    )
    print("| ---: | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: |")
    summary = []
    for key in sorted(by_x):
        needed, spare = [], []
        for i, eps in enumerate(EPSILONS):
            if eps not in by_x[key]:
                continue
            p, x, sources = by_x[key][eps]
            random = next((r[i] for k, r in RANDOM.items() if abs(k - key) <= 0.15), None)
            fit, default = degree(x, eps, 3.0), degree(x, eps, 5.0)
            needed.append(p - fit)
            spare.append(default - p)
            print(
                f"| {x:.2f} | {eps:.0e} | {p} | {', '.join(sources)} | "
                f"{random if random is not None else '–'} | {fit} | {default} | "
                f"{p - fit:+d} | {default - p:+d} |"
            )
        summary.append((key, max(needed), min(spare)))
    print()
    print("| x | largest p_c − fit | smallest CH4.3c2 − p_c |")
    print("| ---: | ---: | ---: |")
    for key, n, s in summary:
        print(f"| {key} | {n:+d} | {s:+d} |")


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


def matrix(args):
    runs = {}
    for path in sorted(glob.glob(os.path.join(args[0], "*.md"))):
        name = os.path.basename(path)[:-3]
        rule = next(r for r in ("calibrated", "ch4-3c2") if name.startswith(r + "-"))
        m2l, case = name[len(rule) + 1 :].split("-", 1)
        lines = open(path).read().split("\n")
        for c in table_rows(lines, "| Workload"):
            eps = float(c["degrees"].split("=")[1])
            key = (m2l, case, c["Workload"], c["N"], c["kD"], eps)
            runs[(rule,) + key] = (
                max(float(c["φ"]), float(c["∇φ"])) / eps,
                c["p by level"],
                float(c["evaluate s"]),
            )
    worst = collections.defaultdict(lambda: [0.0, 0.0, None])
    over = []
    print(
        "| M2L | case | N | kD | ε | p, refit | p, CH4.3c2 | error / ε, refit | CH4.3c2 | "
        "evaluate s, refit / CH4.3c2 |"
    )
    print("| --- | --- | ---: | --- | --- | --- | --- | ---: | ---: | --- |")
    for key in sorted(k[1:] for k in runs if k[0] == "calibrated"):
        new, old = runs[("calibrated",) + key], runs.get(("ch4-3c2",) + key)
        w = worst[key[0]]
        if new[0] > w[0]:
            w[0], w[2] = new[0], key[1:]
        if old:
            w[1] = max(w[1], old[0])
        if new[0] > 1.0:
            over.append(key)
        if old and new[1] != old[1]:
            m2l, case, _, n, kd, eps = key
            print(
                f"| {m2l} | {case} | {n} | {kd} | {eps:.0e} | {new[1]} | {old[1]} | {new[0]:.3f} | "
                f"{old[0]:.3f} | {new[2]:.2f} / {old[2]:.2f} |"
            )
    print()
    n_new = sum(1 for k in runs if k[0] == "calibrated")
    n_old = sum(1 for k in runs if k[0] == "ch4-3c2")
    print(f"{n_new} runs under the refit, {n_old} under CH4.3c2's rule; the largest error / ε:")
    for m2l, (new, old, where) in sorted(worst.items()):
        print(f"- {m2l}: {new:.3f} at {where} (CH4.3c2's rule {old:.3f})")
    print(f"- above ε: {over if over else 'none'}")


def main():
    commands = {"scans": scans, "margin": margin, "matrix": matrix}
    if len(sys.argv) < 3 or sys.argv[1] not in commands:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    commands[sys.argv[1]](sys.argv[2:])


if __name__ == "__main__":
    main()
