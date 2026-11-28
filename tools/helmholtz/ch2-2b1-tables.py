#!/usr/bin/env python3
"""tools/helmholtz/ch2-2b1-tables.py: the tables of Helmholtz CH2.2b1's report
(tetraleaf-validate/results/helmholtz-ch2-2b1.md), RCR's M2L tables only on the levels with an RCR
pair, from the launches of tools/helmholtz/ch2-2b1-runs.sh.

    python3 tools/helmholtz/ch2-2b1-tables.py runs <dir> <case…>     # per case, main against branch
    python3 tools/helmholtz/ch2-2b1-tables.py levels <dir> <case…>   # per case and level
    python3 tools/helmholtz/ch2-2b1-tables.py selftest

<case> is a case of the runs script (`smoke`, `w1-250`, `w1-500`, `w1-1000`, `w3-500`): its four
launches <dir>/<case>-main1.md, -branch1.md, -branch2.md, -main2.md (one run each). `runs`: per case
the build and its parts "RCR tables" and "diagonal tables" (seconds, the slowest rank's), each the
mean of the two launches of a side and the branch's over `main`'s, the RCR tables' MB and
thread-seconds ("RCR tables by level" summed over the levels, of the first launch of each side),
the errors, and whether the four checksums agree. `levels`: per case and level p, q, the parts
built, the thread-seconds and the MB of `main`'s first launch and of the branch's. Python 3.9
standard library only; the parsers are tools/helmholtz/common.py's.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

LEVELS = "RCR tables by level"
SIDES = {"main": ("main1", "main2"), "branch": ("branch1", "branch2")}


def launch(directory, case, run):
    """One launch's run row, "Build" row and level rows."""
    path = Path(directory) / f"{case}-{run}.md"
    rows = common.accuracy_rows(path)
    builds = common.build_rows(path)
    t = common.find_table(common.tables(path), heading=LEVELS)
    return {
        "run": rows[0] if rows else None,
        "build": builds[0] if builds else None,
        "levels": t.rows if t else [],
    }


def total(levels, key):
    return sum(common.number(r[key], 0.0) for r in levels)


def mean(values):
    values = [v for v in values if v is not None]
    return sum(values) / len(values) if values else None


def cell(x, digits=3):
    return "–" if x is None else f"{x:.{digits}f}"


def ratio(a, b):
    return f"{a / b:.2f}" if a is not None and b else "–"


def runs(directory, cases):
    print(
        "| case | build s (main) | build s (branch) | ÷ | RCR tables s (main) "
        "| RCR tables s (branch) | ÷ | diagonal tables s (main) | diagonal tables s (branch) "
        "| RCR MB (main) | RCR MB (branch) | RCR thread-s (main) | RCR thread-s (branch) | φ | ∇φ "
        "| checksums |"
    )
    print("| --- |" + " ---: |" * 14 + " --- |")
    for case in cases:
        got = {run: launch(directory, case, run) for side in SIDES.values() for run in side}
        failed = [r for r, g in got.items() if g["run"] is None or g["build"] is None]
        if failed:
            print(f"| {case} | failed launch: {', '.join(failed)} |" + " – |" * 14)
            continue
        side = {}
        for name, both in SIDES.items():
            side[name] = {
                part: mean(common.number(got[r]["build"][part]) for r in both)
                for part in ("build s", "RCR tables s", "diagonal tables s")
            }
            first = got[both[0]]["levels"]
            side[name]["MB"] = total(first, "MB")
            side[name]["thread-s"] = total(first, "thread-s")
        m, b = side["main"], side["branch"]
        sums = {got[r]["run"]["checksum"] for r in got}
        same = "same" if len(sums) == 1 else "**differ**"
        row = got["branch1"]["run"]
        print(
            f"| {case} | {cell(m['build s'])} | {cell(b['build s'])} | {ratio(b['build s'], m['build s'])} "
            f"| {cell(m['RCR tables s'])} | {cell(b['RCR tables s'])} "
            f"| {ratio(b['RCR tables s'], m['RCR tables s'])} | {cell(m['diagonal tables s'])} "
            f"| {cell(b['diagonal tables s'])} | {m['MB']:.1f} | {b['MB']:.1f} | {m['thread-s']:.3f} "
            f"| {b['thread-s']:.3f} | {row['φ']} | {row['∇φ']} | {same} |"
        )


def levels(directory, cases):
    print(
        "| case | level | k̂_l | p | q | built (main) | built (branch) | thread-s (main) "
        "| thread-s (branch) | MB (main) | MB (branch) |"
    )
    print("| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for case in cases:
        main = launch(directory, case, "main1")["levels"]
        branch = {r["level"]: r for r in launch(directory, case, "branch1")["levels"]}
        for r in main:
            rb = branch.get(r["level"], {})
            print(
                f"| {case} | {r['level']} | {r['k̂_l']} | {r['p']} | {r['q']} | {r['built']} "
                f"| {rb.get('built', '–')} | {r['thread-s']} | {rb.get('thread-s', '–')} "
                f"| {r['MB']} | {rb.get('MB', '–')} |"
            )


SAMPLE = """| Workload | N | kD | degrees | M2L | ranks × threads | build s | evaluate s | φ | ∇φ | checksum |
| --- | ---: | --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| W1 | 10 | 500 | ε = 1e-6 | auto | 1 × 2 | {build} | 0.01 | 1e-7 | 2e-7 | {sum} |

## Build

| Workload | N | kD | degrees | M2L | ranks × threads | build s | RCR tables s | diagonal tables s | evaluate s | checksum |
| --- | ---: | --- | --- | --- | --- | ---: | ---: | ---: | --- | --- |
| W1 | 10 | 500 | ε = 1e-6 | auto | 1 × 2 | {build} | {rcr} | 0.500 | 0.010 | {sum} |

## RCR tables by level

| Workload | N | kD | degrees | M2L | ranks × threads | level | k̂_l | p | q | built | reused | thread-s | MB | store MB |
| --- | ---: | --- | --- | --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| W1 | 10 | 500 | ε = 1e-6 | auto | 1 × 2 | 2 | 125.0 + 0.0i | 129 | 129 | {built} | 0 | {seconds} | {mb} | – |
| W1 | 10 | 500 | ε = 1e-6 | auto | 1 × 2 | 3 | 62.5 + 0.0i | 71 | 71 | 4 | 0 | 0.1000 | 1.00 | – |
"""


def selftest():
    import io
    import tempfile
    from contextlib import redirect_stdout

    with tempfile.TemporaryDirectory() as d:
        for run, build, rcr, built, seconds, mb in [
            ("main1", "2.000", "1.000", 19, "3.0000", "100.00"),
            ("main2", "2.200", "1.200", 19, "3.0000", "100.00"),
            ("branch1", "1.000", "0.010", 4, "0.0100", "10.00"),
            ("branch2", "1.200", "0.012", 4, "0.0100", "10.00"),
        ]:
            Path(d, f"x-{run}.md").write_text(
                SAMPLE.format(
                    build=build, rcr=rcr, built=built, seconds=seconds, mb=mb, sum="00aa"
                ),
                encoding="utf-8",
            )
        buf = io.StringIO()
        with redirect_stdout(buf):
            runs(d, ["x"])
        lines = buf.getvalue().splitlines()
        assert len(lines) == 3, lines
        row = common.cells(lines[2])
        assert row[0] == "x" and row[1] == "2.100" and row[2] == "1.100", row
        assert row[3] == "0.52" and row[4] == "1.100" and row[5] == "0.011", row
        assert row[6] == "0.01" and row[7] == "0.500", row
        assert row[9] == "101.0" and row[10] == "11.0", row
        assert row[11] == "3.100" and row[12] == "0.110", row
        assert row[15] == "same", row
        buf = io.StringIO()
        with redirect_stdout(buf):
            levels(d, ["x"])
        lines = buf.getvalue().splitlines()
        assert len(lines) == 4, lines
        first = common.cells(lines[2])
        assert first[:2] == ["x", "2"] and first[5:7] == ["19", "4"], first
        assert first[9:11] == ["100.00", "10.00"], first
    print("selftest passed")


def main(argv):
    if argv[1:2] == ["selftest"]:
        selftest()
    elif len(argv) >= 4 and argv[1] == "runs":
        runs(argv[2], argv[3:])
    elif len(argv) >= 4 and argv[1] == "levels":
        levels(argv[2], argv[3:])
    else:
        print(__doc__, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
