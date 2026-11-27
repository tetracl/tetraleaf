#!/usr/bin/env python3
"""tools/helmholtz/ch4-9a2a-tables.py: the tables of Helmholtz CH4.9a2a's report (the leaf harmonics
and the real radial functions per point against a leaf point's reference parts;
docs/helmholtz/tasks/CH4.9a2a.md), from the launches of `tools/helmholtz/ch4-9a2a-runs.sh`
(`helmholtz_leaf` outputs `parts-<precision>-p<p>.md` and `lanes-<precision>-b<B>.md`).

    python3 tools/helmholtz/ch4-9a2a-tables.py summary <machine>=<dir> …
    python3 tools/helmholtz/ch4-9a2a-tables.py cases   <machine>=<dir> …
    python3 tools/helmholtz/ch4-9a2a-tables.py lanes   <machine>=<dir> …
    python3 tools/helmholtz/ch4-9a2a-tables.py selftest

- `summary`: per machine, precision and p, the medians over k̂ of the reference's P2M, harmonics,
  radial functions and sum per point, the leaf harmonics (column by column) and the real radial
  functions per point with their ratios to the reference's parts, the sum on the flushed radial
  values, the projected P2M and its ratio, that sum's share of the projected P2M, and over k̂ the largest harmonics check and the flushed
  radial values.
- `cases`: the same per k̂ (every case).
- `lanes`: per machine, precision, p and lane count, the leaf harmonics (`values`, columns) and the
  real radial functions per point and the projected P2M.

No timing is asserted; the numbers go into tetraleaf-validate/results/helmholtz-ch4-9a2a.md.
Python 3.9 standard library only.
"""

import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402


def rows(arg, prefix):
    """`machine=dir`: the machine and every row of its launches `<prefix>-*.md`, as dicts with
    the precision and (for `lanes`) the lane count added, sorted."""
    machine, where = arg.split("=", 1)
    out = []
    for path in sorted(Path(where).glob(f"{prefix}-*.md")):
        parts = path.stem.split("-")
        t = common.find_table(common.tables(path), first="p")
        for r in t.rows if t else []:
            r = dict(r)
            r["precision"] = parts[1]
            if prefix == "lanes":
                r["lanes"] = parts[2].lstrip("b")
            out.append(r)
    out.sort(key=lambda r: (r["precision"], int(r["p"]), int(r.get("lanes", 0)),
                            float(r["k̂"])))
    return machine, out


def flushed(text):
    """The flushed radial values of a "radial check" cell: 0 for "bit for bit"."""
    return 0 if text.strip() == "bit for bit" else int(text.split()[0])


def parts(r):
    """The numbers of one row."""
    n = lambda k: common.number(r[k])  # noqa: E731
    p2m, harm, radial, total_sum = n("P2M µs"), n("harmonics µs"), n("radial µs"), n("sum µs")
    columns, real = n("leaf columns µs"), n("real radial µs")
    projected = n("projected P2M µs")
    leaf_sum = n("sum, leaf radial µs")
    return {
        "p2m": p2m, "harmonics": harm, "radial": radial, "sum": total_sum,
        "columns": columns, "values": n("leaf harmonics µs"), "real": real,
        "projected": projected, "harm ratio": columns / harm, "radial ratio": real / radial,
        "leaf sum": leaf_sum, "projected ratio": projected / p2m,
        "sum share": leaf_sum / projected,
        "check": n("harmonics check"), "flushed": flushed(r["radial check"]),
    }


def show(header, body):
    print("| " + " | ".join(header) + " |")
    print("| " + " | ".join("---" for _ in header) + " |")
    for row in body:
        print("| " + " | ".join(row) + " |")
    print()


def line(machine, precision, p, s, extra=()):
    return [machine, precision, str(p), *extra, f"{s['p2m']:.3f}", f"{s['harmonics']:.3f}",
            f"{s['columns']:.3f}", f"{s['harm ratio']:.3f}", f"{s['radial']:.3f}",
            f"{s['real']:.3f}", f"{s['radial ratio']:.2f}", f"{s['sum']:.3f}",
            f"{s['leaf sum']:.3f}", f"{s['projected']:.3f}", f"{s['projected ratio']:.2f}",
            f"{100 * s['sum share']:.0f} %", f"{s['check']:.2f}", str(s["flushed"])]


HEADER = ["P2M µs", "harmonics µs", "leaf µs", "÷", "radial µs", "real µs", "÷", "sum µs",
          "sum, leaf radial µs", "projected µs", "÷ P2M", "sum share", "harmonics check", "flushed"]


def cmd_summary(args):
    body = []
    for arg in args:
        machine, rs = rows(arg, "parts")
        groups = {}
        for r in rs:
            groups.setdefault((r["precision"], int(r["p"])), []).append(parts(r))
        for (precision, p), ss in sorted(groups.items()):
            med = {k: statistics.median(s[k] for s in ss) for k in ss[0]}
            med["check"] = max(s["check"] for s in ss)
            med["flushed"] = sum(s["flushed"] for s in ss)
            body.append(line(machine, precision, p, med))
    show(["machine", "precision", "p", *HEADER], body)


def cmd_cases(args):
    body = []
    for arg in args:
        machine, rs = rows(arg, "parts")
        for r in rs:
            body.append(line(machine, r["precision"], r["p"], parts(r), (r["k̂"],)))
    show(["machine", "precision", "p", "k̂", *HEADER], body)


def cmd_lanes(args):
    body = []
    for arg in args:
        machine, rs = rows(arg, "lanes")
        for r in rs:
            s = parts(r)
            body.append([machine, r["precision"], r["p"], r["lanes"], f"{s['values']:.3f}",
                         f"{s['columns']:.3f}", f"{s['real']:.3f}", f"{s['projected']:.3f}",
                         f"{s['projected ratio']:.2f}"])
    show(["machine", "precision", "p", "lanes", "leaf values µs", "leaf columns µs",
          "real radial µs", "projected µs", "÷ P2M"], body)


def selftest():
    r = {"P2M µs": "2", "harmonics µs": "1", "radial µs": "0.5", "sum µs": "0.2",
         "leaf columns µs": "0.1", "leaf harmonics µs": "0.12", "real radial µs": "0.25",
         "projected P2M µs": "0.85", "harmonics check": "1.02", "radial check": "12 flushed",
         "sum, leaf radial µs": "0.17"}
    s = parts(r)
    assert s["harm ratio"] == 0.1 and s["radial ratio"] == 0.5 and s["flushed"] == 12
    assert abs(s["projected ratio"] - 0.425) < 1e-12
    assert flushed("bit for bit") == 0
    assert abs(s["sum share"] - 0.2) < 1e-12
    print("selftest: ok")


def main():
    commands = {"summary": cmd_summary, "cases": cmd_cases, "lanes": cmd_lanes}
    if len(sys.argv) >= 2 and sys.argv[1] == "selftest":
        selftest()
    elif len(sys.argv) >= 3 and sys.argv[1] in commands:
        commands[sys.argv[1]](sys.argv[2:])
    else:
        print(__doc__)
        sys.exit(2)


if __name__ == "__main__":
    main()
