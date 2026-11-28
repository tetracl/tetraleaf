#!/usr/bin/env python3
"""tools/helmholtz/ch2-2a-tables.py: the tables of Helmholtz CH2.2a's report
(tetraleaf-validate/results/helmholtz-ch2-2a.md), the RCR tables per level in a store, from the
launches of tools/helmholtz/ch2-2a-runs.sh.

    python3 tools/helmholtz/ch2-2a-tables.py sweep <dir> <name>    # per kD, store against cold
    python3 tools/helmholtz/ch2-2a-tables.py levels <dir> <name>   # per kD and level
    python3 tools/helmholtz/ch2-2a-tables.py selftest

<name> is a stage of the runs script (`real`, `rotation`, `yukawa`, `smoke`): its launches
<dir>/<name>-store.md (the sweep sharing one store) and <dir>/<name>-cold.md (every case without
one). `sweep`: per kD the build and its "RCR tables" part (seconds, the slowest rank's) of both
launches, the parts built and reused and their thread-seconds from "RCR tables by level" summed
over the levels, the tables' MB, the store's MB after the case, and whether the checksums agree.
`levels`: per kD and level k̂_l, p, q, the parts built and reused with the store, the thread-seconds
with and without it, and the level's MB. Python 3.9 standard library only; the parsers are
tools/helmholtz/common.py's.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

LEVELS = "RCR tables by level"


def runs(path):
    """Per kD cell, in order: the run row, the "Build" row and the level rows."""
    found = common.tables(path)
    out = {}
    for r in common.accuracy_rows(path):
        out.setdefault(r["kD"], {"run": r, "build": None, "levels": []})
    for r in common.build_rows(path):
        out[r["kD"]]["build"] = r
    t = common.find_table(found, heading=LEVELS)
    for r in t.rows if t else []:
        out[r["kD"]]["levels"].append(r)
    return out


def total(levels, key):
    return sum(common.number(r[key], 0.0) for r in levels)


def sweep(directory, name):
    store = runs(Path(directory) / f"{name}-store.md")
    cold = runs(Path(directory) / f"{name}-cold.md")
    print(
        "| kD | levels | build s (store) | build s (cold) | RCR tables s (store) "
        "| RCR tables s (cold) | ÷ | parts built | reused | thread-s (store) | thread-s (cold) "
        "| tables MB | store MB | checksums |"
    )
    print("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |")
    for kd, s in store.items():
        c = cold.get(kd)
        if c is None:
            continue
        rs, rc = (common.number(x["build"]["RCR tables s"]) for x in (s, c))
        ratio = f"{rs / rc:.2f}" if rs is not None and rc else "–"
        store_mb = s["levels"][0]["store MB"] if s["levels"] else "–"
        same = "same" if s["run"]["checksum"] == c["run"]["checksum"] else "**differ**"
        print(
            f"| {kd} | {len(s['levels'])} | {s['build']['build s']} | {c['build']['build s']} "
            f"| {rs:.3f} | {rc:.3f} | {ratio} | {total(s['levels'], 'built'):.0f} "
            f"| {total(s['levels'], 'reused'):.0f} | {total(s['levels'], 'thread-s'):.3f} "
            f"| {total(c['levels'], 'thread-s'):.3f} | {total(s['levels'], 'MB'):.1f} | {store_mb} "
            f"| {same} |"
        )


def levels(directory, name):
    store = runs(Path(directory) / f"{name}-store.md")
    cold = runs(Path(directory) / f"{name}-cold.md")
    print("| kD | level | k̂_l | p | q | built | reused | thread-s (store) | thread-s (cold) | MB |")
    print("| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for kd, s in store.items():
        c = cold.get(kd, {"levels": []})
        cold_levels = {r["level"]: r for r in c["levels"]}
        for r in s["levels"]:
            rc = cold_levels.get(r["level"], {})
            print(
                f"| {kd} | {r['level']} | {r['k̂_l']} | {r['p']} | {r['q']} | {r['built']} "
                f"| {r['reused']} | {r['thread-s']} | {rc.get('thread-s', '–')} | {r['MB']} |"
            )


SAMPLE = """| Workload | N | kD | degrees | charges | m | P2P | M2L | ranks × threads | levels | leaf level | p by level | leaves | max leaf | build s | evaluate s | φ | ∇φ | checksum | ÷ Rotation |
| --- | ---: | --- | --- | --- | ---: | --- | --- | --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| W1 | 10 | 15 | ε = 1e-6 | complex | 1 | neon | rotation | 1 × 2 | 3 | 2 | 0×2, 9 | 8 | 3 | 0.10 | 0.01 | 1e-7 | 1e-7 | 00000000000000aa | – |
| W1 | 10 | 30 | ε = 1e-6 | complex | 1 | neon | rotation | 1 × 2 | 4 | 3 | 0×2, 9, 9 | 8 | 3 | 0.20 | 0.01 | 1e-7 | 1e-7 | 00000000000000bb | – |

## Build

| Workload | N | kD | degrees | M2L | ranks × threads | build s | RCR tables s | evaluate s | checksum |
| --- | ---: | --- | --- | --- | --- | ---: | ---: | --- | --- |
| W1 | 10 | 15 | ε = 1e-6 | rotation | 1 × 2 | 0.100 | 0.050 | 0.010 | 00000000000000aa |
| W1 | 10 | 30 | ε = 1e-6 | rotation | 1 × 2 | 0.200 | {rcr} | 0.010 | 00000000000000bb |

## RCR tables by level

| Workload | N | kD | degrees | M2L | ranks × threads | level | k̂_l | p | q | built | reused | thread-s | MB | store MB |
| --- | ---: | --- | --- | --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| W1 | 10 | 15 | ε = 1e-6 | rotation | 1 × 2 | 2 | 1.0 + 0.0i | 9 | 9 | 10 | 0 | 0.0400 | 1.00 | {store} |
| W1 | 10 | 30 | ε = 1e-6 | rotation | 1 × 2 | 2 | 2.0 + 0.0i | 9 | 9 | 10 | 0 | 0.0400 | 1.00 | {store} |
| W1 | 10 | 30 | ε = 1e-6 | rotation | 1 × 2 | 3 | 1.0 + 0.0i | 9 | 9 | {built} | {reused} | {seconds} | 1.00 | {store} |
"""


def selftest():
    import io
    import tempfile
    from contextlib import redirect_stdout

    with tempfile.TemporaryDirectory() as d:
        Path(d, "x-store.md").write_text(
            SAMPLE.format(rcr="0.010", store="2.00", built=0, reused=10, seconds="0.0000"),
            encoding="utf-8",
        )
        Path(d, "x-cold.md").write_text(
            SAMPLE.format(rcr="0.050", store="–", built=10, reused=0, seconds="0.0400"),
            encoding="utf-8",
        )
        buf = io.StringIO()
        with redirect_stdout(buf):
            sweep(d, "x")
        lines = buf.getvalue().splitlines()
        assert len(lines) == 4, lines
        row = common.cells(lines[3])
        assert row[0] == "30" and row[1] == "2", row
        assert row[6] == "0.20", row  # 0.010 / 0.050
        assert row[7] == "10" and row[8] == "10", row  # built 10 + 0, reused 0 + 10
        assert row[9] == "0.040" and row[10] == "0.080", row
        assert row[12] == "2.00" and row[13] == "same", row
        buf = io.StringIO()
        with redirect_stdout(buf):
            levels(d, "x")
        lines = buf.getvalue().splitlines()
        assert len(lines) == 5, lines
        last = common.cells(lines[-1])
        assert last[:2] == ["30", "3"] and last[5:9] == ["0", "10", "0.0000", "0.0400"], last
    print("selftest passed")


def main(argv):
    if argv[1:2] == ["selftest"]:
        selftest()
    elif len(argv) == 4 and argv[1] == "sweep":
        sweep(argv[2], argv[3])
    elif len(argv) == 4 and argv[1] == "levels":
        levels(argv[2], argv[3])
    else:
        print(__doc__, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
