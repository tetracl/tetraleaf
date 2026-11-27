#!/usr/bin/env python3
"""tools/helmholtz/ch4-9b1-tables.py: the tables of Helmholtz CH4.9b1's report (the leaf harmonics'
and the real radial functions' derivatives per point against a leaf point's reference L2P parts;
docs/helmholtz/tasks/CH4.9b1.md), from the launches of `tools/helmholtz/ch4-9b1-runs.sh`
(`helmholtz_leaf` outputs `parts-<precision>-p<p>.md`).

    python3 tools/helmholtz/ch4-9b1-tables.py summary <machine>=<dir> …
    python3 tools/helmholtz/ch4-9b1-tables.py cases   <machine>=<dir> …
    python3 tools/helmholtz/ch4-9b1-tables.py selftest

- `summary`: per machine, precision and p, the medians over k̂ of the reference's L2P with gradients,
  its harmonics with derivatives and its radial functions with derivatives per point, the leaf
  derivative columns and the real radial functions with derivatives per point with their ratios to
  the reference's parts, the L2P sum with gradients on the leaf routines' values, the projected L2P
  with gradients and its ratio, that sum's share of the projection, and over k̂ the largest
  derivatives check and the radial derivatives that differ beside flushed values.
- `cases`: the same per k̂ (every case).

No timing is asserted; the numbers go into tetraleaf-validate/results/helmholtz-ch4-9b1.md.
Python 3.9 standard library only.
"""

import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402


def rows(arg):
    """`machine=dir`: the machine and every row of its launches `parts-*.md`, as dicts with the
    precision added, sorted."""
    machine, where = arg.split("=", 1)
    out = []
    for path in sorted(Path(where).glob("parts-*.md")):
        precision = path.stem.split("-")[1]
        t = common.find_table(common.tables(path), first="p")
        for r in t.rows if t else []:
            r = dict(r)
            r["precision"] = precision
            out.append(r)
    out.sort(key=lambda r: (r["precision"], int(r["p"]), float(r["k̂"])))
    return machine, out


def differ(text):
    """The radial derivatives of a "radial ′ check" cell that differ: 0 for "bit for bit"."""
    return 0 if text.strip() == "bit for bit" else int(text.split()[0])


def parts(r):
    """The numbers of one row."""
    n = lambda k: common.number(r[k])  # noqa: E731
    l2p, harm, radial = n("L2P ∇ µs"), n("harmonics ′ µs"), n("radial ′ µs")
    columns, real = n("leaf columns ′ µs"), n("real radial ′ µs")
    leaf_sum, projected = n("sum ∇, leaf values µs"), n("projected L2P ∇ µs")
    return {
        "l2p": l2p, "harmonics": harm, "radial": radial, "columns": columns, "real": real,
        "harm ratio": columns / harm, "radial ratio": real / radial, "leaf sum": leaf_sum,
        "projected": projected, "projected ratio": projected / l2p,
        "sum share": leaf_sum / projected,
        "check": n("derivatives check"), "differ": differ(r["radial ′ check"]),
    }


def show(header, body):
    print("| " + " | ".join(header) + " |")
    print("| " + " | ".join("---" for _ in header) + " |")
    for row in body:
        print("| " + " | ".join(row) + " |")
    print()


def line(machine, precision, p, s, extra=()):
    return [machine, precision, str(p), *extra, f"{s['l2p']:.3f}", f"{s['harmonics']:.3f}",
            f"{s['columns']:.3f}", f"{s['harm ratio']:.3f}", f"{s['radial']:.3f}",
            f"{s['real']:.3f}", f"{s['radial ratio']:.2f}", f"{s['leaf sum']:.3f}",
            f"{s['projected']:.3f}", f"{s['projected ratio']:.2f}",
            f"{100 * s['sum share']:.0f} %", f"{s['check']:.2f}", str(s["differ"])]


HEADER = ["L2P ∇ µs", "harmonics ′ µs", "leaf ′ µs", "÷", "radial ′ µs", "real ′ µs", "÷",
          "sum ∇, leaf values µs", "projected µs", "÷ L2P ∇", "sum share", "derivatives check",
          "radial ′ differ"]


def cmd_summary(args):
    body = []
    for arg in args:
        machine, rs = rows(arg)
        groups = {}
        for r in rs:
            groups.setdefault((r["precision"], int(r["p"])), []).append(parts(r))
        for (precision, p), ss in sorted(groups.items()):
            med = {k: statistics.median(s[k] for s in ss) for k in ss[0]}
            med["check"] = max(s["check"] for s in ss)
            med["differ"] = sum(s["differ"] for s in ss)
            body.append(line(machine, precision, p, med))
    show(["machine", "precision", "p", *HEADER], body)


def cmd_cases(args):
    body = []
    for arg in args:
        machine, rs = rows(arg)
        for r in rs:
            body.append(line(machine, r["precision"], r["p"], parts(r), (r["k̂"],)))
    show(["machine", "precision", "p", "k̂", *HEADER], body)


def selftest():
    r = {"L2P ∇ µs": "4", "harmonics ′ µs": "2", "radial ′ µs": "1", "leaf columns ′ µs": "0.2",
         "real radial ′ µs": "0.25", "sum ∇, leaf values µs": "0.5",
         "projected L2P ∇ µs": "1", "derivatives check": "0.51",
         "radial ′ check": "43 beside flushed"}
    s = parts(r)
    assert s["harm ratio"] == 0.1 and s["radial ratio"] == 0.25 and s["differ"] == 43
    assert s["projected ratio"] == 0.25 and s["sum share"] == 0.5
    assert differ("bit for bit") == 0
    print("selftest: ok")


def main():
    commands = {"summary": cmd_summary, "cases": cmd_cases}
    if len(sys.argv) >= 2 and sys.argv[1] == "selftest":
        selftest()
    elif len(sys.argv) >= 3 and sys.argv[1] in commands:
        commands[sys.argv[1]](sys.argv[2:])
    else:
        print(__doc__)
        sys.exit(2)


if __name__ == "__main__":
    main()
