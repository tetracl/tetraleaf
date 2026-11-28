#!/usr/bin/env python3
"""tools/helmholtz/ch4-9b2-tables.py: the report tables of Helmholtz CH4.9b2, L2P with gradients by
the batched kernel (docs/helmholtz/tasks/CH4.9b2.md), from the outputs of
`tools/helmholtz/ch4-9b2-runs.sh`.

    python3 tools/helmholtz/ch4-9b2-tables.py runs <out-dir>     # per case: L2P, P2M, evaluation, errors
    python3 tools/helmholtz/ch4-9b2-tables.py levels <out-dir>   # per case and level: L2P and P2M
    python3 tools/helmholtz/ch4-9b2-tables.py worst <out-dir>    # the constant density's errors
    python3 tools/helmholtz/ch4-9b2-tables.py parts <out-dir>    # helmholtz_leaf: the L2P kernel per point
    python3 tools/helmholtz/ch4-9b2-tables.py native <out-dir>   # the baseline build against target-cpu=native
    python3 tools/helmholtz/ch4-9b2-tables.py selftest

Each case of `levels` has four launches, kernel-1, reference-1, reference-2, kernel-2 (ABBA): a
setting's time is the mean of its two launches, its spread the two launches' ratio. `kernel` is the
library's default (P2M and L2P by the batched kernels: after CH4.9a–b), `reference` both by the
per-point reference calls (before CH4.9a–b). L2P and P2M per run are the sums over the levels of
"Translations and leaf expansions"' "L2P s" and "P2M s" (the second evaluation's, the slowest
rank's per level). The errors (φ, ∇φ) are each setting's first launch's; the two launches of a
setting must give one checksum (asserted). Python 3.9 standard library only.
"""

import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

# The levels table leaves out levels whose reference L2P and P2M both take less than this (no
# occupied leaves: the level call's overhead, at the clock's resolution).
MIN_SECONDS = 1e-3

CASES = [
    ("w1-250-1e-6", "W1", "250", "1e-6"),
    ("w1-500-1e-6", "W1", "500", "1e-6"),
    ("w1-1000-1e-6", "W1", "1000", "1e-6"),
    ("w1-500-1e-3", "W1", "500", "1e-3"),
    ("w1-500-1e-8", "W1", "500", "1e-8"),
    ("w3-500-1e-6", "W3 (10⁶)", "500", "1e-6"),
]


def run(path):
    """One launch: its run row and its per-level {level: (L2P s, P2M s, p)}."""
    rows = common.accuracy_rows(path)
    if len(rows) != 1:
        raise SystemExit(f"{path}: {len(rows)} run rows, expected 1")
    t = common.find_table(common.tables(path), heading="Translations and leaf expansions")
    levels = {}
    for r in t.rows if t else []:
        if r["level"].isdigit():
            levels[int(r["level"])] = (common.number(r["L2P s"], 0.0),
                                       common.number(r["P2M s"], 0.0), r["p"])
    return rows[0], levels


def launches(out, name):
    """The case's four launches as {setting: [(row, levels), (row, levels)]}."""
    got = {}
    for setting in ("kernel", "reference"):
        pair = [run(Path(out) / f"{name}-{setting}-{i}.md") for i in (1, 2)]
        sums = {r["checksum"] for r, _ in pair}
        if len(sums) != 1:
            raise SystemExit(f"{name}, {setting}: the launches' checksums differ: {sums}")
        got[setting] = pair
    return got


def total(levels, part):
    return sum(v[part] for v in levels.values())


def mean(values):
    return sum(values) / len(values)


def runs(out):
    print("| case | ε | L2P s, kernel | L2P s, reference | L2P ÷ | P2M s, kernel | P2M s, reference "
          "| P2M ÷ | P2M + L2P share, reference | P2M + L2P share, kernel | evaluate s, kernel "
          "| evaluate s, reference | evaluate ÷ (spread) | φ ÷ | ∇φ ÷ |")
    print("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- "
          "| ---: | ---: |")
    for name, w, kd, eps in CASES:
        got = launches(out, name)
        k, r = got["kernel"], got["reference"]
        lk, lr = (mean([total(lv, 0) for _, lv in s]) for s in (k, r))
        pk, pr = (mean([total(lv, 1) for _, lv in s]) for s in (k, r))
        ek = [float(row["evaluate s"]) for row, _ in k]
        er = [float(row["evaluate s"]) for row, _ in r]
        spread = max(max(ek) / min(ek), max(er) / min(er))
        phi = float(k[0][0]["φ"]) / float(r[0][0]["φ"])
        grad = float(k[0][0]["∇φ"]) / float(r[0][0]["∇φ"])
        print(f"| {w}, kD = {kd} | {eps} | {lk:.4f} | {lr:.4f} | {lk / lr:.3f} | {pk:.4f} | "
              f"{pr:.4f} | {pk / pr:.3f} | {(lr + pr) / mean(er):.3f} | "
              f"{(lk + pk) / mean(ek):.3f} | {mean(ek):.3f} | {mean(er):.3f} | "
              f"{mean(ek) / mean(er):.3f} ({spread:.3f}) | {phi:.4f} | {grad:.4f} |")


def levels(out):
    print("| case | ε | level | p | L2P s, kernel | L2P s, reference | ÷ | P2M s, kernel "
          "| P2M s, reference | ÷ |")
    print("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for name, w, kd, eps in CASES:
        got = launches(out, name)
        k, r = got["kernel"], got["reference"]
        for level in sorted(k[0][1]):
            lk, lr = (mean([lv[level][0] for _, lv in s]) for s in (k, r))
            pk, pr = (mean([lv[level][1] for _, lv in s]) for s in (k, r))
            if lr < MIN_SECONDS and pr < MIN_SECONDS:
                continue
            ratio = lambda a, b: f"{a / b:.3f}" if b >= MIN_SECONDS else "–"  # noqa: E731
            p = k[0][1][level][2]
            print(f"| {w}, kD = {kd} | {eps} | {level} | {p} | {lk:.4f} | {lr:.4f} | "
                  f"{ratio(lk, lr)} | {pk:.4f} | {pr:.4f} | {ratio(pk, pr)} |")


def worst(out):
    print("| case | φ, kernel | φ, reference | ÷ | ∇φ, kernel | ∇φ, reference | ÷ |")
    print("| --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for w in ("w1", "w2"):
        k, _ = run(Path(out) / f"worst-{w}-500-1e-6-kernel.md")
        r, _ = run(Path(out) / f"worst-{w}-500-1e-6-reference.md")
        pk, pr, gk, gr = (float(x) for x in (k["φ"], r["φ"], k["∇φ"], r["∇φ"]))
        print(f"| {w.upper()}, kD = 500, 10⁻⁶, constant | {pk:.3e} | {pr:.3e} | {pk / pr:.4f} | "
              f"{gk:.3e} | {gr:.3e} | {gk / gr:.4f} |")


def leaf_rows(path):
    """`helmholtz_leaf`'s rows of one launch."""
    t = common.find_table(common.tables(path), first="p")
    return t.rows if t else []


def launches_of(out, prefix):
    """{p: rows} of the launches `<prefix>-p<p>.md`, by degree."""
    got = {}
    for path in Path(out).glob(f"{prefix}-p*.md"):
        rows = leaf_rows(path)
        if rows:
            got[int(path.stem.rsplit("p", 1)[1])] = rows
    return dict(sorted(got.items()))


def med(rows, key):
    return statistics.median(float(r[key]) for r in rows)


def parts(out):
    print("| precision | p | L2P ∇ µs, reference | ∇ kernel µs | ÷ reference | ÷ projected "
          "| ∇ kernel sum µs | sum share | sum ns per term | sum GFLOP/s "
          "| reference sum ∇, leaf values µs | kernel sum ÷ it | L2P (no ∇) ÷ reference "
          "| worst check |")
    print("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: "
          "| ---: |")
    for precision in ("f64", "f32"):
        for p, rows in launches_of(out, f"parts-{precision}").items():
            ker, ks = med(rows, "L2P ∇ kernel µs"), med(rows, "∇ kernel sum µs")
            ls = med(rows, "sum ∇, leaf values µs")
            check = max(float(r["L2P kernel check"]) for r in rows)
            print(f"| {precision} | {p} | {med(rows, 'L2P ∇ µs'):.3f} | {ker:.3f} | "
                  f"{med(rows, '∇ kernel ÷ L2P ∇'):.3f} | {med(rows, '∇ kernel ÷ projected'):.2f} | "
                  f"{ks:.3f} | {ks / ker:.2f} | {med(rows, '∇ kernel sum ns per term'):.3f} | "
                  f"{med(rows, '∇ kernel sum GFLOP/s'):.2f} | {ls:.3f} | {ks / ls:.2f} | "
                  f"{med(rows, 'L2P kernel ÷ L2P'):.3f} | {check:.1e} |")
    print()
    print("Medians over k̂ ∈ {0.5, 2, 4, 8, 16} (the ratios the medians of each k̂'s ratio); "
          "\"∇ kernel sum\" the kernel with gradients less its routines' parts; GFLOP/s by the "
          "sum's (p + 1) (12.5 p + 8) operations per target; \"worst check\" the largest of the "
          "kernel's largest difference from the reference over its largest value.")


def native(out):
    print("| p | ∇ kernel µs, baseline | native | ÷ | sum ns per term, baseline | native | ÷ "
          "| sum GFLOP/s, baseline | native | P2M kernel µs, baseline | native | ÷ |")
    print("| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    base, nat = launches_of(out, "parts-f64"), launches_of(out, "native-f64")
    for p in sorted(set(base) & set(nat)):
        b, n = base[p], nat[p]
        kb, kn = med(b, "L2P ∇ kernel µs"), med(n, "L2P ∇ kernel µs")
        sb, sn = med(b, "∇ kernel sum ns per term"), med(n, "∇ kernel sum ns per term")
        pb, pn = med(b, "P2M kernel µs"), med(n, "P2M kernel µs")
        print(f"| {p} | {kb:.3f} | {kn:.3f} | {kn / kb:.2f} | {sb:.3f} | {sn:.3f} | {sn / sb:.2f} "
              f"| {med(b, '∇ kernel sum GFLOP/s'):.2f} | {med(n, '∇ kernel sum GFLOP/s'):.2f} | "
              f"{pb:.3f} | {pn:.3f} | {pn / pb:.2f} |")
    print()
    print("f64, medians over k̂; \"native\" the same source built with `-C target-cpu=native`.")


def selftest():
    common.selftest()
    levels = {2: (0.5, 0.25, "8"), 3: (1.5, 0.75, "8")}
    assert total(levels, 0) == 2.0 and total(levels, 1) == 1.0
    assert mean([1.0, 3.0]) == 2.0
    assert med([{"x": "1"}, {"x": "3"}, {"x": "2"}], "x") == 2.0
    print("selftest passed")


def main():
    if len(sys.argv) == 2 and sys.argv[1] == "selftest":
        selftest()
        return
    stages = {"runs": runs, "levels": levels, "worst": worst, "parts": parts, "native": native}
    if len(sys.argv) != 3 or sys.argv[1] not in stages:
        raise SystemExit(__doc__)
    stages[sys.argv[1]](sys.argv[2])


if __name__ == "__main__":
    main()
