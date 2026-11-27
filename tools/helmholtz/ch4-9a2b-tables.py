#!/usr/bin/env python3
"""tools/helmholtz/ch4-9a2b-tables.py: the report tables of Helmholtz CH4.9a2b, P2M by the batched
kernel (docs/helmholtz/tasks/CH4.9a2b.md), from the outputs of `tools/helmholtz/ch4-9a2b-runs.sh`.

    python3 tools/helmholtz/ch4-9a2b-tables.py runs <out-dir>     # per case: P2M, evaluation, errors
    python3 tools/helmholtz/ch4-9a2b-tables.py levels <out-dir>   # per case and level: P2M
    python3 tools/helmholtz/ch4-9a2b-tables.py worst <out-dir>    # the constant density's errors
    python3 tools/helmholtz/ch4-9a2b-tables.py parts <out-dir>    # helmholtz_leaf per point
    python3 tools/helmholtz/ch4-9a2b-tables.py selftest

Each case of `levels` has four launches, kernel-1, reference-1, reference-2, kernel-2 (ABBA): a
setting's time is the mean of its two launches, its spread the two launches' ratio. P2M per run is
the sum over the levels of "Translations and leaf expansions"' "P2M s" (the second evaluation's,
the slowest rank's per level). The errors (φ, ∇φ) are each setting's first launch's; the two launches
of a setting must give one checksum (asserted). Python 3.9 standard library only.
"""

import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

# The levels table leaves out levels whose reference P2M takes less than this (no occupied
# leaves: the level call's overhead, at the clock's resolution).
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
    """One launch: its run row and its per-level P2M seconds {level: s}."""
    rows = common.accuracy_rows(path)
    if len(rows) != 1:
        raise SystemExit(f"{path}: {len(rows)} run rows, expected 1")
    found = common.tables(path)
    t = common.find_table(found, heading="Translations and leaf expansions")
    levels = {}
    for r in t.rows if t else []:
        if r["level"].isdigit():
            levels[int(r["level"])] = (common.number(r["P2M s"], 0.0), r["p"])
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


def p2m_total(levels):
    return sum(s for s, _ in levels.values())


def mean(values):
    return sum(values) / len(values)


def runs(out):
    print("| case | ε | P2M s, kernel | P2M s, reference | P2M ÷ | P2M share of the evaluation, "
          "reference | evaluate s, kernel | evaluate s, reference | evaluate ÷ (spread) | φ ÷ | ∇φ ÷ "
          "| checksums |")
    print("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | --- |")
    for name, w, kd, eps in CASES:
        got = launches(out, name)
        k, r = got["kernel"], got["reference"]
        pk = mean([p2m_total(lv) for _, lv in k])
        pr = mean([p2m_total(lv) for _, lv in r])
        ek = [float(row["evaluate s"]) for row, _ in k]
        er = [float(row["evaluate s"]) for row, _ in r]
        spread = max(max(ek) / min(ek), max(er) / min(er))
        phi = float(k[0][0]["φ"]) / float(r[0][0]["φ"])
        grad = float(k[0][0]["∇φ"]) / float(r[0][0]["∇φ"])
        same = "differ" if k[0][0]["checksum"] != r[0][0]["checksum"] else "same"
        print(f"| {w}, kD = {kd} | {eps} | {pk:.4f} | {pr:.4f} | {pk / pr:.3f} | "
              f"{pr / mean(er):.3f} | {mean(ek):.3f} | {mean(er):.3f} | "
              f"{mean(ek) / mean(er):.3f} ({spread:.3f}) | {phi:.4f} | {grad:.4f} | {same} |")


def levels(out):
    print("| case | ε | level | p | P2M s, kernel | P2M s, reference | ÷ |")
    print("| --- | --- | ---: | ---: | ---: | ---: | ---: |")
    for name, w, kd, eps in CASES:
        got = launches(out, name)
        k, r = got["kernel"], got["reference"]
        for level in sorted(k[0][1]):
            sk = mean([lv[level][0] for _, lv in k])
            sr = mean([lv[level][0] for _, lv in r])
            # Levels without occupied leaves: the level call's overhead alone.
            if sr < MIN_SECONDS:
                continue
            p = k[0][1][level][1]
            print(f"| {w}, kD = {kd} | {eps} | {level} | {p} | {sk:.4f} | {sr:.4f} | {sk / sr:.3f} |")


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


def parts(out):
    print("| precision | p | reference P2M µs | kernel µs | kernel ÷ reference | kernel ÷ projected "
          "| kernel sum µs | sum share | reference sum, leaf values µs | kernel sum ÷ it "
          "| worst check |")
    print("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for precision in ("f64", "f32"):
        for path in sorted(Path(out).glob(f"parts-{precision}-p*.md"),
                           key=lambda p: int(p.stem.rsplit("p", 1)[1])):
            rows = leaf_rows(path)
            if not rows:
                continue
            med = lambda key: statistics.median(float(r[key]) for r in rows)  # noqa: E731
            ref, ker = med("P2M µs"), med("P2M kernel µs")
            ks = med("kernel sum µs")
            ls = med("sum, leaf values µs")
            check = max(float(r["kernel check"]) for r in rows)
            print(f"| {precision} | {rows[0]['p']} | {ref:.3f} | {ker:.3f} | "
                  f"{med('kernel ÷ P2M'):.3f} | {med('kernel ÷ projected'):.2f} | {ks:.3f} | "
                  f"{ks / ker:.2f} | {ls:.3f} | {ks / ls:.2f} | {check:.1e} |")
    print()
    print("Medians over k̂ ∈ {0.5, 2, 4, 8, 16} (the ratios the medians of each k̂'s ratio); "
          "\"kernel sum\" the kernel less its routines' parts; \"worst check\" the largest of the "
          "kernel's largest difference from the reference over its largest value.")


def selftest():
    common.selftest()
    rows = [{"P2M s": "0.5"}, {"P2M s": "1.5"}]
    assert sum(common.number(r["P2M s"], 0.0) for r in rows) == 2.0
    assert mean([1.0, 3.0]) == 2.0
    print("selftest passed")


def main():
    if len(sys.argv) == 2 and sys.argv[1] == "selftest":
        selftest()
        return
    if len(sys.argv) != 3 or sys.argv[1] not in ("runs", "levels", "worst", "parts"):
        raise SystemExit(__doc__)
    {"runs": runs, "levels": levels, "worst": worst, "parts": parts}[sys.argv[1]](sys.argv[2])


if __name__ == "__main__":
    main()
