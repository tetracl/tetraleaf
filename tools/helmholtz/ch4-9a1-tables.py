#!/usr/bin/env python3
"""tools/helmholtz/ch4-9a1-tables.py: the tables of Helmholtz CH4.9a1's report (a Helmholtz leaf
point's P2M and L2P time by part; docs/helmholtz/tasks/CH4.9a1.md), from the launches of
`tools/helmholtz/ch4-9a1-runs.sh` (`helmholtz_leaf` outputs `parts-<precision>-p<p>.md`).

    python3 tools/helmholtz/ch4-9a1-tables.py shares  <machine>=<dir> …
    python3 tools/helmholtz/ch4-9a1-tables.py summary <machine>=<dir> …
    python3 tools/helmholtz/ch4-9a1-tables.py selftest

- `shares`: per machine, precision, p and k̂ the P2M per point and the shares of its parts (the
  radial functions, the harmonics, the sum, the exponential), and the L2P with gradients per
  point and the shares of its parts (the radial functions with derivatives, the harmonics with
  derivatives, the rest), with the check of the reassembled P2M.
- `summary`: per machine, precision and p the median over k̂ of the same shares, and the ns per
  term of P2M and of L2P with gradients.

No timing is asserted; the numbers go into tetraleaf-validate/results/helmholtz-ch4-9a1.md.
Python 3.9 standard library only.
"""

import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402


def rows(arg):
    """`machine=dir`: the machine and every row of its launches' tables, as dicts with the
    precision added, sorted by precision, p and k̂."""
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


def shares(r):
    """The parts' shares of P2M and of L2P with gradients, and the ns per term."""
    n = lambda k: common.number(r[k])  # noqa: E731
    p2m, l2p = n("P2M µs"), n("L2P ∇ µs")
    return {
        "radial": n("radial µs") / p2m,
        "harmonics": n("harmonics µs") / p2m,
        "sum": n("sum µs") / p2m,
        "exp": n("exp µs") / p2m,
        "radial′": n("radial ′ µs") / l2p,
        "harmonics′": n("harmonics ′ µs") / l2p,
        "rest": n("rest ∇ µs") / l2p,
        "p2m": p2m,
        "l2p": l2p,
        "p2m ns": n("P2M ns per term"),
        "l2p ns": n("L2P ∇ ns per term"),
    }


def pct(x):
    return f"{100 * x:.0f} %"


def show(header, body):
    print("| " + " | ".join(header) + " |")
    print("| " + " | ".join("---" for _ in header) + " |")
    for row in body:
        print("| " + " | ".join(row) + " |")
    print()


def cmd_shares(args):
    header = ["machine", "precision", "p", "k̂", "P2M µs", "radial", "harmonics", "sum", "exp",
              "check", "L2P ∇ µs", "radial ′", "harmonics ′", "rest"]
    body = []
    for arg in args:
        machine, rs = rows(arg)
        for r in rs:
            s = shares(r)
            body.append([machine, r["precision"], r["p"], r["k̂"], f"{s['p2m']:.3f}",
                         pct(s["radial"]), pct(s["harmonics"]), pct(s["sum"]), pct(s["exp"]),
                         r["check"], f"{s['l2p']:.3f}", pct(s["radial′"]), pct(s["harmonics′"]),
                         pct(s["rest"])])
    show(header, body)


def cmd_summary(args):
    header = ["machine", "precision", "p", "P2M µs", "P2M ns per term", "radial", "harmonics",
              "sum", "L2P ∇ µs", "L2P ∇ ns per term", "radial ′", "harmonics ′", "rest"]
    body = []
    for arg in args:
        machine, rs = rows(arg)
        groups = {}
        for r in rs:
            groups.setdefault((r["precision"], int(r["p"])), []).append(shares(r))
        for (precision, p), ss in sorted(groups.items()):
            med = lambda k: statistics.median(s[k] for s in ss)  # noqa: E731
            body.append([machine, precision, str(p), f"{med('p2m'):.3f}", f"{med('p2m ns'):.1f}",
                         pct(med("radial")), pct(med("harmonics")), pct(med("sum")),
                         f"{med('l2p'):.3f}", f"{med('l2p ns'):.1f}", pct(med("radial′")),
                         pct(med("harmonics′")), pct(med("rest"))])
    show(header, body)


def selftest():
    r = {"P2M µs": "2", "radial µs": "1", "harmonics µs": "0.5", "sum µs": "0.2", "exp µs": "0.1",
         "L2P ∇ µs": "4", "radial ′ µs": "1", "harmonics ′ µs": "2", "rest ∇ µs": "1",
         "P2M ns per term": "3", "L2P ∇ ns per term": "6"}
    s = shares(r)
    assert s["radial"] == 0.5 and s["harmonics"] == 0.25 and s["rest"] == 0.25
    assert pct(0.254) == "25 %"
    print("selftest: ok")


def main():
    commands = {"shares": cmd_shares, "summary": cmd_summary}
    if len(sys.argv) >= 2 and sys.argv[1] == "selftest":
        selftest()
    elif len(sys.argv) >= 3 and sys.argv[1] in commands:
        commands[sys.argv[1]](sys.argv[2:])
    else:
        print(__doc__)
        sys.exit(2)


if __name__ == "__main__":
    main()
