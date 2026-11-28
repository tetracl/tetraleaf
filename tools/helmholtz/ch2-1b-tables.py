#!/usr/bin/env python3
"""tools/helmholtz/ch2-1b-tables.py: the tables of Helmholtz CH2.1b's report
(tetraleaf-validate/results/helmholtz-ch2-1b.md) from the launches of tools/helmholtz/ch2-1b-runs.sh.

    python3 tools/helmholtz/ch2-1b-tables.py accuracy <acc-dir>      # the acceptance's errors
    python3 tools/helmholtz/ch2-1b-tables.py determinism <acc-dir>   # checksums across splits
    python3 tools/helmholtz/ch2-1b-tables.py extra <acc-dir>         # Dense and complex charges
    python3 tools/helmholtz/ch2-1b-tables.py checksums <acc-dir>     # main against the branch
    python3 tools/helmholtz/ch2-1b-tables.py time <time-dir> [prefix] # D4: main against the branch
    python3 tools/helmholtz/ch2-1b-tables.py kinds <time-dir> [prefix] # D4: per kind, the same
    python3 tools/helmholtz/ch2-1b-tables.py selftest

Python 3.9 standard library only; the parsers are tools/helmholtz/common.py's.
"""

import re
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

KAPPAS = ["0.01", "0.1", "1", "10", "100", "1000"]
EPS = ["1e-3", "1e-6", "1e-8"]
WORKLOADS = [
    ("w1", "W1 (N = 10⁵)"),
    ("w3", "W3 (N = 10⁵)"),
    ("w4-cube", "W4 cube"),
    ("w4-sphere", "W4 sphere"),
    ("w4-plummer", "W4 Plummer"),
    ("w4-clusters", "W4 clusters"),
]


def e(x):
    return f"{x:.2e}" if x is not None else "–"


def worst(rows):
    """The larger of φ's and ∇φ's error of each row, by ε."""
    out = {}
    for r in rows:
        phi, grad = common.number(r["φ"]), common.number(r["∇φ"])
        out[r["ε"]] = (phi, grad, r["p by level"])
    return out


def accuracy(d):
    d = Path(d)
    print("| Workload | ε | " + " | ".join(f"κD = {k}" for k in KAPPAS) + " | max / ε |")
    print("| --- | --- | " + " | ".join("---:" for _ in KAPPAS) + " | ---: |")
    overall = 0.0
    for tag, name in WORKLOADS:
        found = {k: worst(common.accuracy_rows(d / f"acc-{tag}-{k}.md")) for k in KAPPAS}
        for eps in EPS:
            cells, ratio = [], 0.0
            for k in KAPPAS:
                phi, grad, _ = found[k].get(eps, (None, None, None))
                m = max(x for x in (phi, grad) if x is not None) if phi is not None else None
                ratio = max(ratio, m / float(eps)) if m is not None else ratio
                cells.append(common.bold(e(m), m is not None and m > float(eps)))
            overall = max(overall, ratio)
            print(f"| {name} | {eps} | " + " | ".join(cells) + f" | {ratio:.1e} |")
    print()
    print("| W3 (N = 10⁶), ε = 1e-6 | " + " | ".join(f"κD = {k}" for k in KAPPAS) + " |")
    print("| --- | " + " | ".join("---:" for _ in KAPPAS) + " |")
    cells = []
    for k in KAPPAS:
        rows = common.accuracy_rows(d / f"acc-w3large-{k}.md")
        phi, grad, _ = worst(rows).get("1e-6", (None, None, None))
        m = max(phi, grad) if phi is not None else None
        overall = max(overall, m / 1e-6) if m is not None else overall
        cells.append(e(m))
    print("| max(φ, ∇φ) | " + " | ".join(cells) + " |")
    print()
    print(f"Largest error over ε: {overall:.2e}")
    print()
    print("Degrees by level (the plateau, x = 0) and evaluate s at κD = 1:")
    for tag, name in WORKLOADS:
        rows = common.accuracy_rows(d / f"acc-{tag}-1.md")
        parts = [f"{r['ε']}: {r['p by level']}, {r['evaluate s']} s" for r in rows]
        print(f"- {name}: " + "; ".join(parts))


def checksum_rows(path):
    return [(r["ranks × threads"], r["checksum"]) for r in common.accuracy_rows(path)]


def determinism(d):
    d = Path(d)
    print("| case | 1 × 1 | 1 × 40 | 4 × 10 | the matrix's 1 × 40 | equal |")
    print("| --- | --- | --- | --- | --- | --- |")
    for k in ["1", "100"]:
        for tag in ["w1", "w4-plummer"]:
            one = checksum_rows(d / f"det-{tag}-{k}-1x.md")
            four = checksum_rows(d / f"det-{tag}-{k}-4x.md")
            matrix = [
                c
                for r in common.accuracy_rows(d / f"acc-{tag}-{k}.md")
                if r["ε"] == "1e-6"
                for c in [r["checksum"]]
            ]
            sums = [c for _, c in one] + [c for _, c in four] + matrix
            print(
                f"| {tag}, κD = {k}, 1e-6 | "
                + " | ".join(sums)
                + f" | {'yes' if len(set(sums)) == 1 else '**no**'} |"
            )


def extra(d):
    d = Path(d)
    print("| case | M2L | φ | ∇φ | checksum | ÷ Rotation |")
    print("| --- | --- | ---: | ---: | --- | --- |")
    for k in ["1", "30"]:
        for r in common.accuracy_rows(d / f"dense-w4-cube-{k}.md"):
            print(
                f"| W4 cube, κD = {k}, 1e-6 | {r['M2L']} | {r['φ']} | {r['∇φ']} | "
                f"{r['checksum']} | {r['÷ Rotation']} |"
            )
    for k in ["1", "100"]:
        for r in common.accuracy_rows(d / f"complex-w1-{k}.md"):
            print(
                f"| W1 (N = 10⁵), κD = {k}, 1e-6, complex charges | {r['M2L']} | {r['φ']} | "
                f"{r['∇φ']} | {r['checksum']} | – |"
            )


def checksums(d):
    d = Path(d)
    print("| case | M2L | main | branch | equal |")
    print("| --- | --- | --- | --- | --- |")
    bad = 0
    for f in sorted(d.glob("sum-*-main.md")):
        other = f.with_name(f.name.replace("-main.md", "-branch.md"))
        a = [r["checksum"] for r in common.accuracy_rows(f)]
        b = [r["checksum"] for r in common.accuracy_rows(other)]
        same = a == b and len(a) > 0
        bad += not same
        case = f.name[len("sum-") : -len("-main.md")]
        m2l = "rotation" if case.endswith("-rotation") else "defaults"
        case = case.rsplit("-", 1)[0].removesuffix("-auto")
        print(f"| {case} | {m2l} | {', '.join(a)} | {', '.join(b)} | {'yes' if same else '**no**'} |")
    print()
    print(f"Pairs that differ: {bad}")


def time(d, prefix="time"):
    d = Path(d)
    cases = sorted(
        {re.sub(r"-(main|branch)-\d[ab]\.md$", "", f.name) for f in d.glob(f"{prefix}-*.md")}
    )
    print(
        "| case | N | build s main | branch | ratio | evaluate s main | branch | ratio | "
        "evaluations | checksums equal |"
    )
    print("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |")
    ratios = []
    for c in cases:
        runs = {"main": [], "branch": []}
        sums = set()
        n = "–"
        for f in sorted(d.glob(f"{c}-*.md")):
            who = "main" if "-main-" in f.name else "branch"
            rows = common.accuracy_rows(f)
            build = common.build_rows(f)
            if not rows or not build:
                continue
            n = rows[0]["N"]
            # The "Build" table's times (to the millisecond; the first table rounds to 10 ms).
            b = build[0]
            runs[who].append((common.number(b["build s"]), common.number(b["evaluate s"])))
            sums.add(rows[0]["checksum"])
        med = {w: [statistics.median(x[i] for x in v) for i in (0, 1)] for w, v in runs.items() if v}
        if len(med) < 2:
            continue
        rb = med["branch"][0] / med["main"][0]
        re_ = med["branch"][1] / med["main"][1]
        ratios.append((rb, re_))
        print(
            f"| {c.removeprefix(prefix + '-')} | {n} | {med['main'][0]:.3f} | {med['branch'][0]:.3f} | "
            f"{rb:.3f} | {med['main'][1]:.3f} | {med['branch'][1]:.3f} | {re_:.3f} | "
            f"{len(runs['main'])} + {len(runs['branch'])} | {'yes' if len(sums) == 1 else '**no**'} |"
        )
        spread = {
            w: (min(x[1] for x in v), max(x[1] for x in v)) for w, v in runs.items() if v
        }
        print(
            f"|   evaluate s range | | | | | {spread['main'][0]:.3f}–{spread['main'][1]:.3f} | "
            f"{spread['branch'][0]:.3f}–{spread['branch'][1]:.3f} | | | |"
        )
    if ratios:
        print()
        print(
            f"Geometric mean of the ratios (branch / main): build "
            f"{common.gmean(r[0] for r in ratios):.3f}, evaluate {common.gmean(r[1] for r in ratios):.3f}"
        )


def selftest():
    assert e(1.234e-5) == "1.23e-05"
    assert e(None) == "–"
    print("selftest passed")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    {
        "accuracy": accuracy,
        "determinism": determinism,
        "extra": extra,
        "checksums": checksums,
        "time": time,
        "selftest": lambda: selftest(),
    }[sys.argv[1]](*sys.argv[2:])
