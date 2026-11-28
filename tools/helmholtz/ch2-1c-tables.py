#!/usr/bin/env python3
"""tools/helmholtz/ch2-1c-tables.py: the tables of Helmholtz CH2.1c's report
(tetraleaf-validate/results/helmholtz-ch2-1c.md), Yukawa's batched P2M, from the launches of
tools/helmholtz/ch2-1c-runs.sh.

    python3 tools/helmholtz/ch2-1c-tables.py accuracy <acc-dir>   # kernel ÷ reference, and ε
    python3 tools/helmholtz/ch2-1c-tables.py worst <acc-dir>      # the constant density
    python3 tools/helmholtz/ch2-1c-tables.py runs <time-dir>      # P2M and the evaluation, ABBA
    python3 tools/helmholtz/ch2-1c-tables.py levels <time-dir>    # P2M per level
    python3 tools/helmholtz/ch2-1c-tables.py selftest

`accuracy`: per launch pair (kernel, reference) and ε the larger of φ's and ∇φ's error ratio
kernel ÷ reference, the worst over the cases, and each error over ε. `runs` and `levels`: each
case's four launches kernel-1, reference-1, reference-2, kernel-2; a setting's time the mean of its
two launches; P2M the sum over the levels of "Translations and leaf expansions"' "P2M s" (the
second evaluation's). Python 3.9 standard library only; the parsers are tools/helmholtz/common.py's.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

KAPPAS = ["0.01", "0.1", "1", "10", "100", "1000"]
WORKLOADS = [
    ("w1", "W1 (N = 10⁵)"),
    ("w3", "W3 (N = 10⁵)"),
    ("w4-cube", "W4 cube"),
    ("w4-sphere", "W4 sphere"),
    ("w4-plummer", "W4 Plummer"),
    ("w4-clusters", "W4 clusters"),
    ("w3large", "W3 (N = 10⁶)"),
]
TIME_CASES = [
    ("w1-1", "W1 (N = 2·10⁶)", "1"),
    ("w1-100", "W1 (N = 2·10⁶)", "100"),
    ("w3-1", "W3 (N = 10⁶)", "1"),
    ("w3-100", "W3 (N = 10⁶)", "100"),
]
# Levels whose reference P2M takes less than this are left out (no occupied leaves).
MIN_SECONDS = 1e-3


def by_eps(path):
    return {r["ε"]: r for r in common.accuracy_rows(path)}


def ratio(k, r):
    """kernel ÷ reference of φ and of ∇φ (None where the reference's is 0 or missing)."""
    out = []
    for key in ("φ", "∇φ"):
        a, b = common.number(k[key]), common.number(r[key])
        out.append(a / b if a is not None and b else None)
    return out


def pairs(d, prefix):
    """(name, kernel rows by ε, reference rows by ε) of every launch pair with this prefix."""
    for f in sorted(Path(d).glob(f"{prefix}-*-kernel.md")):
        name = f.name[: -len("-kernel.md")]
        yield name, by_eps(f), by_eps(f.with_name(f"{name}-reference.md"))


def accuracy(d):
    d = Path(d)
    print("| Workload | ε | " + " | ".join(f"κD = {k}" for k in KAPPAS) + " |")
    print("| --- | --- | " + " | ".join("---:" for _ in KAPPAS) + " |")
    lo, hi, over, count, cases, same = 1e9, 0.0, 0.0, 0, 0, 0
    for tag, name in WORKLOADS:
        found = {k: (by_eps(d / f"acc-{tag}-{k}-kernel.md"), by_eps(d / f"acc-{tag}-{k}-reference.md"))
                 for k in KAPPAS}
        epsilons = sorted({e for k in KAPPAS for e in found[k][0]}, key=float, reverse=True)
        for eps in epsilons:
            cells = []
            for k in KAPPAS:
                kern, ref = found[k]
                if eps not in kern or eps not in ref:
                    cells.append("–")
                    continue
                cases += 1
                rs = [x for x in ratio(kern[eps], ref[eps]) if x is not None]
                worst = max(rs, key=lambda x: abs(x - 1.0))
                lo, hi = min(lo, min(rs)), max(hi, max(rs))
                errors = [common.number(kern[eps][key]) for key in ("φ", "∇φ")]
                over = max(over, max(e for e in errors if e is not None) / float(eps))
                same += kern[eps]["checksum"] == ref[eps]["checksum"]
                count += max(rs) > 1.01
                cells.append(common.bold(f"{worst:.4f}", max(rs) > 1.01))
            print(f"| {name} | {eps} | " + " | ".join(cells) + " |")
    print()
    print(f"{cases} cases: kernel ÷ reference {lo:.4f}–{hi:.4f} (φ and ∇φ), {count} above 1.01; "
          f"the largest error over ε {over:.3f}; checksums equal in {same}.")
    print()
    print("| Complex charges, W1 (N = 10⁵), 1e-6 | φ, kernel | φ, reference | ÷ | checksum, kernel "
          "| reference |")
    print("| --- | ---: | ---: | ---: | --- | --- |")
    for k in ("1", "100"):
        kern = by_eps(d / f"complex-w1-{k}-kernel.md")["1e-6"]
        ref = by_eps(d / f"complex-w1-{k}-reference.md")["1e-6"]
        print(f"| κD = {k} | {kern['φ']} | {ref['φ']} | {ratio(kern, ref)[0]:.4f} | "
              f"{kern['checksum']} | {ref['checksum']} |")


def worst(d):
    print("| case | φ, kernel | φ, reference | ÷ | ∇φ, kernel | ∇φ, reference | ÷ |")
    print("| --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for name, kern, ref in pairs(d, "worst"):
        k, r = kern["1e-6"], ref["1e-6"]
        rp, rg = ratio(k, r)
        print(f"| {k['Workload']}, κD = {name.rsplit('-', 1)[1]}, 10⁻⁶, constant | {k['φ']} | "
              f"{r['φ']} | {rp:.4f} | {k['∇φ']} | {r['∇φ']} | {rg:.4f} |")


def run(path):
    """One launch: its run row and its per-level (P2M s, p)."""
    rows = common.accuracy_rows(path)
    if len(rows) != 1:
        raise SystemExit(f"{path}: {len(rows)} run rows, expected 1")
    t = common.find_table(common.tables(path), heading="Translations and leaf expansions")
    levels = {}
    for r in t.rows if t else []:
        if r["level"].isdigit():
            levels[int(r["level"])] = (common.number(r["P2M s"], 0.0), r["p"])
    return rows[0], levels


def launches(d, tag):
    got = {}
    for setting in ("kernel", "reference"):
        pair = [run(Path(d) / f"levels-{tag}-{setting}-{i}.md") for i in (1, 2)]
        if len({r["checksum"] for r, _ in pair}) != 1:
            raise SystemExit(f"{tag}, {setting}: the launches' checksums differ")
        got[setting] = pair
    return got


def mean(values):
    values = list(values)
    return sum(values) / len(values)


def runs(d):
    print("| case | κD | P2M s, kernel | P2M s, reference | P2M ÷ | P2M share, reference | "
          "evaluate s, kernel | evaluate s, reference | evaluate ÷ (spread) | φ ÷ | ∇φ ÷ |")
    print("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: |")
    for tag, name, kappa in TIME_CASES:
        got = launches(d, tag)
        k, r = got["kernel"], got["reference"]
        pk = mean(sum(s for s, _ in lv.values()) for _, lv in k)
        pr = mean(sum(s for s, _ in lv.values()) for _, lv in r)
        ek = [float(row["evaluate s"]) for row, _ in k]
        er = [float(row["evaluate s"]) for row, _ in r]
        spread = max(max(ek) / min(ek), max(er) / min(er))
        rp, rg = ratio(k[0][0], r[0][0])
        print(f"| {name} | {kappa} | {pk:.4f} | {pr:.4f} | {pk / pr:.3f} | {pr / mean(er):.3f} | "
              f"{mean(ek):.3f} | {mean(er):.3f} | {mean(ek) / mean(er):.3f} ({spread:.3f}) | "
              f"{rp:.4f} | {rg:.4f} |")


def levels(d):
    print("| case | κD | level | p | P2M s, kernel | P2M s, reference | ÷ |")
    print("| --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for tag, name, kappa in TIME_CASES:
        got = launches(d, tag)
        k, r = got["kernel"], got["reference"]
        for level in sorted(k[0][1]):
            sk = mean(lv[level][0] for _, lv in k)
            sr = mean(lv[level][0] for _, lv in r)
            if sr < MIN_SECONDS:
                continue
            print(f"| {name} | {kappa} | {level} | {k[0][1][level][1]} | {sk:.4f} | {sr:.4f} | "
                  f"{sk / sr:.3f} |")


def selftest():
    common.selftest()
    assert ratio({"φ": "2e-6", "∇φ": "0"}, {"φ": "1e-6", "∇φ": "0"}) == [2.0, None]
    assert mean([1.0, 3.0]) == 2.0
    print("selftest passed")


def main():
    if len(sys.argv) == 2 and sys.argv[1] == "selftest":
        selftest()
        return
    commands = {"accuracy": accuracy, "worst": worst, "runs": runs, "levels": levels}
    if len(sys.argv) != 3 or sys.argv[1] not in commands:
        raise SystemExit(__doc__)
    commands[sys.argv[1]](sys.argv[2])


if __name__ == "__main__":
    main()
