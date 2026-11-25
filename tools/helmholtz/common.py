#!/usr/bin/env python3
"""tools/helmholtz/common.py: what the report scripts of tools/helmholtz/ share (Helmholtz
CHA.10a, docs/helmholtz/tasks/CHA.10a.md): the Markdown row splitting, the parsers of the harness
examples' outputs and the geometric mean.

    python3 tools/helmholtz/common.py selftest

A new report script imports it rather than defining its own `cells`, `tables` or parsers:

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import common  # noqa: E402

(the explicit path keeps the import working under `python3 -I`, which leaves the script's directory
off `sys.path`). The scripts written before CHA.10a keep their own copies, so their reports'
reproduction lines stay valid.

What it reads (one launch per file, the harness examples of tetraleaf-validate):
- every output: `tables` (each Markdown table, under its `## ` heading) and `find_table`;
- `helmholtz_accuracy`: `accuracy_rows` (its first table, one row per run, the ε added) and
  `build_rows` (its "Build" table by ε);
- `helmholtz_boxes`: `boxes_runs` and `boxes_run` (the run rows: degrees, §9.1's ℓ2 errors, the
  worst box over ε, the build and the evaluation);
- `helmholtz_fmm3d compare`: `compare_rows` (one row per ε and code, ours and FMM3D's variants)
  and `build_rows` (its "Build" table by ε).

Cells are split on every `|`: a cell that holds one (`helmholtz_boxes`' worst-box tables have
`|v| / r` in their header) shifts the columns after it in `tables`' row dicts.

Python 3.9 standard library only.
"""

import math
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, NamedTuple, Optional

# The tolerances of D1 as the harness prints them.
EPSILONS = ("1e-3", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8")

_EPS = re.compile(r"ε = ([0-9.]+(?:e[+-]?[0-9]+)?)")


def cells(line):
    """The cells of a Markdown table row, stripped: `| a | b |` gives ["a", "b"]."""
    return [c.strip() for c in line.strip().strip("|").split("|")]


def number(text, default=None):
    """`text` as a float, or `default` where it is none (`–`, `-`, empty, None)."""
    try:
        return float(text)
    except (TypeError, ValueError):
        return default


def eps_token(text):
    """The tolerance of a label such as `ε = 1e-6` or `ε = 1e-6, defaults` (`"1e-6"`), or None."""
    m = _EPS.search(text or "")
    return m.group(1) if m else None


def bold(text, over):
    """`text` in bold where `over` holds (a bound exceeded)."""
    return f"**{text}**" if over else text


def gmean(values: Iterable[float]) -> float:
    """The geometric mean of the positive finite values, NaN where there are none."""
    xs = [v for v in values if v is not None and math.isfinite(v) and v > 0]
    return math.exp(sum(map(math.log, xs)) / len(xs)) if xs else float("nan")


class Table(NamedTuple):
    """One Markdown table of an output."""

    heading: str  # the text of the last `## ` heading before it; "" before any
    header: List[str]
    rows: List[Dict[str, str]]  # each row by the header's cells


def _is_rule(row):
    text = "".join(row)
    return "-" in text and set(text) <= set("-: ")


def tables(path) -> List[Table]:
    """Every Markdown table of an output, in order. A table is a run of lines starting with `|`;
    its first line is the header, a rule line (`| --- | ---: |`) is skipped."""
    out, heading, current = [], "", None
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.startswith("|"):
            current = None
            if line.startswith("## "):
                heading = line[3:].strip()
            continue
        row = cells(line)
        if current is None:
            current = Table(heading, row, [])
            out.append(current)
        elif not _is_rule(row):
            current.rows.append(dict(zip(current.header, row)))
    return out


def find_table(found: List[Table], heading=None, first=None) -> Optional[Table]:
    """The first table under `heading` (any where None) whose header starts with the cell
    `first` (any where None), or None."""
    for t in found:
        if (heading is None or t.heading == heading) and (first is None or t.header[:1] == [first]):
            return t
    return None


def accuracy_rows(path) -> List[Dict[str, str]]:
    """`helmholtz_accuracy`'s run rows (its first table, before any heading), each with the key
    "ε" added from its `degrees` cell (None where the run names no tolerance)."""
    t = find_table(tables(path), heading="", first="Workload")
    rows = []
    for r in t.rows if t else []:
        r = dict(r)
        r["ε"] = eps_token(r.get("degrees"))
        rows.append(r)
    return rows


def build_rows(path) -> Dict[str, Dict[str, str]]:
    """The "Build" table by ε, the build by part in seconds: `helmholtz_fmm3d compare`'s (keyed
    by its `ε` cell) and `helmholtz_accuracy`'s (by its `degrees` cell)."""
    t = find_table(tables(path), heading="Build")
    out = {}
    for r in t.rows if t else []:
        eps = r["ε"] if "ε" in r else eps_token(r.get("degrees"))
        if eps:
            out[eps] = r
    return out


class BoxesRun(NamedTuple):
    """A run row of `helmholtz_boxes`."""

    label: str  # "ε = 1e-6, defaults"
    eps: Optional[str]  # "1e-6"
    p: str  # the degrees by level, as printed
    leaves: Optional[float]
    s91: List[float]  # §9.1's relative ℓ2 errors at the samples: φ, ∇φ
    all: List[Optional[float]]  # the relative ℓ2 errors over the reference's targets: φ, ∇φ
    worst: List[float]  # the worst box of φ and ∇φ over ε
    build: float  # seconds
    evaluate: float  # seconds


def boxes_runs(path) -> List[BoxesRun]:
    """`helmholtz_boxes`' run rows (`| ε = … |`, nine cells), in order; rows whose errors or
    times do not parse are left out."""
    runs = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.startswith("| ε = "):
            continue
        c = cells(line)
        if len(c) != 9:
            continue
        try:
            s91 = [float(v) for v in c[3].split(",")]
            worst = [float(v) for v in c[5].split(",")]
            build, evaluate = [float(v) for v in c[8].split(",")]
        except ValueError:
            continue
        runs.append(BoxesRun(c[0], eps_token(c[0]), c[1], number(c[2]), s91,
                             [number(v) for v in c[4].split(",")], worst, build, evaluate))
    return runs


def boxes_run(path) -> Optional[BoxesRun]:
    """`helmholtz_boxes`' first run row, or None."""
    runs = boxes_runs(path)
    return runs[0] if runs else None


class CompareRow(NamedTuple):
    """A row of `helmholtz_fmm3d compare`'s first table: ours, or one of FMM3D's variants."""

    eps: str  # "1e-6"
    code: str  # "ours (avx512)", "FMM3D default", "FMM3D fastker"
    threads: str  # "4 × 10" (ranks × threads) or "40"
    p: str  # our degrees by level, or FMM3D's "eps = …"
    build: Optional[float]  # seconds; None for FMM3D
    run: Optional[float]  # seconds
    phi: Optional[float]  # §9.1's relative ℓ2 error of φ at the samples
    grad: Optional[float]  # and of ∇φ
    memory: Optional[float]  # GB: ours the largest node's, FMM3D's its peak
    checksum: str


def compare_rows(path) -> List[CompareRow]:
    """`helmholtz_fmm3d compare`'s rows (its first table), in order."""
    t = find_table(tables(path), heading="", first="Workload")
    rows = []
    for r in t.rows if t else []:
        if "code" not in r or "ε" not in r:
            continue
        memory = (r.get("memory GB") or "").split(" ")[0]
        rows.append(CompareRow(r["ε"], r["code"], r.get("threads", ""), r.get("p by level", ""),
                               number(r.get("build s")), number(r.get("run s")), number(r.get("φ")),
                               number(r.get("∇φ")), number(memory), r.get("checksum", "")))
    return rows


def selftest():
    """The parsers on lines of each format; raises AssertionError on a failure."""
    import tempfile

    assert cells("| a | b c |  |") == ["a", "b c", ""]
    assert number("1.5e-3") == 1.5e-3 and number("–") is None and number(None, 0.0) == 0.0
    assert eps_token("ε = 1e-6, defaults") == "1e-6" and eps_token("ε = 1e-3") == "1e-3"
    assert eps_token("p = 12") is None and eps_token(None) is None
    assert bold("1.2", True) == "**1.2**" and bold("0.8", False) == "0.8"
    assert abs(gmean([2.0, 8.0]) - 4.0) < 1e-12
    assert gmean([4.0, float("nan"), 0.0, None]) == 4.0 and math.isnan(gmean([]))
    with tempfile.TemporaryDirectory() as d:
        boxes = Path(d) / "boxes.md"
        boxes.write_text(
            "# The box study\n\nRuns: the degrees by level; …\n\n"
            "| run | p by level | leaves | §9.1 φ, ∇φ | all φ, ∇φ | worst box / ε | worst / ℓ2 | boxes > ε φ, ∇φ | build, evaluate s |\n"
            "| --- | --- | ---: | --- | --- | --- | --- | --- | --- |\n"
            "| ε = 1e-6, defaults | 0×2, 17, 15×2 | 232 | 1.65e-8, 1.38e-7 | 1.74e-8, 1.70e-7 | 0.46, 7.14 | 27, 42 | 0, 16 of 1080 | 0.04, 0.05 |\n"
            "| ε = 1e-8, defaults | 0×2, 30 | 232 | –, – | –, – | –, – | –, – | – | 0.04, 0.05 |\n"
            "\n| # | cell | centre | error / ε | targets | boundary | |v| / r, levels 2–1 | level parts | M2L parts |\n"
            "| ---: | --- | --- | ---: | ---: | --- | --- | --- | --- |\n"
            "| 1 | (1, 1, 1) | (0.750, 0.750, 0.750) | 7.14 | 99 | yes |  |  |  |\n", encoding="utf-8")
        run = boxes_run(boxes)
        assert run == BoxesRun("ε = 1e-6, defaults", "1e-6", "0×2, 17, 15×2", 232.0, [1.65e-8, 1.38e-7],
                               [1.74e-8, 1.70e-7], [0.46, 7.14], 0.04, 0.05), run
        assert len(boxes_runs(boxes)) == 1  # the row without errors is left out
        found = tables(boxes)
        assert [len(t.rows) for t in found] == [2, 1] and found[1].header[:2] == ["#", "cell"]
        assert found[1].rows[0]["error / ε"] == "7.14"  # before the `|v|` cell, unshifted

        accuracy = Path(d) / "accuracy.md"
        accuracy.write_text(
            "# Helmholtz FMM against the sampled direct sum\n\nDensity: random\n\n"
            "| Workload | N | kD | degrees | build s | evaluate s | φ | ∇φ | checksum |\n"
            "| --- | ---: | --- | --- | ---: | ---: | ---: | ---: | --- |\n"
            "| W1 | 79578 | 100 + 10i | ε = 1e-3 | 0.09 | 0.11 | 3.18e-5 | 5.19e-5 | 6c8d3555ce5aad7b |\n"
            "| W1 | 79578 | 100 + 10i | ε = 1e-6 | 0.08 | 0.17 | 1.26e-7 | 3.09e-7 | 6b9d13e8d84a707b |\n"
            "\n## Empty boxes\n\n| Workload | N | degrees | rule |\n| --- | ---: | --- | --- |\n"
            "| W1 | 79578 | ε = 1e-3 | on |\n"
            "\n## Build\n\n| Workload | degrees | build | octree |\n| --- | --- | ---: | ---: |\n"
            "| W1 | ε = 1e-3 | 0.0880 | 0.0151 |\n| W1 | ε = 1e-6 | 0.0760 | 0.0150 |\n", encoding="utf-8")
        rows = accuracy_rows(accuracy)
        assert [(r["ε"], r["φ"], r["checksum"]) for r in rows] == [
            ("1e-3", "3.18e-5", "6c8d3555ce5aad7b"), ("1e-6", "1.26e-7", "6b9d13e8d84a707b")], rows
        assert find_table(tables(accuracy), heading="Empty boxes").rows == [
            {"Workload": "W1", "N": "79578", "degrees": "ε = 1e-3", "rule": "on"}]
        assert {e: r["octree"] for e, r in build_rows(accuracy).items()} == {"1e-3": "0.0151", "1e-6": "0.0150"}
        assert find_table(tables(accuracy), heading="Leaves") is None

        compare = Path(d) / "compare.md"
        compare.write_text(
            "# Helmholtz FMM beside FMM3D: W1 at kD = 30\n\nMachine: …\n\n"
            "| Workload | N | kD | ε | code | threads | p by level | build s | run s | φ | ∇φ | all φ | all ∇φ | worst box φ / ε, ∇φ / ε | worst box / all (φ / ∇φ) of boxes | checksum | memory GB |\n"
            "| --- | ---: | --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | ---: |\n"
            "| W1 | 7162 | 30 | 1e-3 | ours (avx512) | 4 × 10 | 0×2, 11 | 0.11 | 0.01 | 5.07e-5 | 8.05e-5 | 5.27e-5 | 9.64e-5 | 0.58 / 1.60 | 10.9 / 16.6 of 56 | 4c1b1fff1a08174a | 0.87 (rank ≤ 0.22) |\n"
            "| W1 | 7162 | 30 | 1e-3 | FMM3D fastker | 40 | eps = 1e-3 | – | 0.03 | 7.48e-5 | 7.37e-5 | 7.57e-5 | 7.53e-5 | 0.35 / 0.35 | 4.6 / 4.7 of 56 | – | 0.02 |\n"
            "\n## Build\n\nOurs, seconds.\n\n| ε | build | settings | evaluation | first evaluation |\n"
            "| --- | ---: | ---: | ---: | ---: |\n| 1e-3 | 0.1043 | 0.0032 | 0.0953 | 0.0935 |\n", encoding="utf-8")
        rows = compare_rows(compare)
        assert rows == [
            CompareRow("1e-3", "ours (avx512)", "4 × 10", "0×2, 11", 0.11, 0.01, 5.07e-5, 8.05e-5, 0.87,
                       "4c1b1fff1a08174a"),
            CompareRow("1e-3", "FMM3D fastker", "40", "eps = 1e-3", None, 0.03, 7.48e-5, 7.37e-5, 0.02, "–"),
        ], rows
        assert build_rows(compare)["1e-3"]["first evaluation"] == "0.0935"
        assert compare_rows(boxes) == [] and build_rows(boxes) == {}


if __name__ == "__main__":
    if sys.argv[1:] == ["selftest"]:
        selftest()
        print("selftest: ok")
    else:
        sys.exit(__doc__)
