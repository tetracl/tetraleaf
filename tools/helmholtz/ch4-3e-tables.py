#!/usr/bin/env python3
"""The tables of Helmholtz CH4.3e's report (tetraleaf-validate/results/helmholtz-ch4-3e.md).

Reads the launches of tools/helmholtz/ch4-3e-runs.sh: each case `<case>` in two files,
`fwd-<case>.md` (c_leaf 1, 2, 4, 8) and `rev-<case>.md` (8, 4, 2, 1), every run under the
library's defaults with the time rule named.

    ch4-3e-tables.py census <machine>=<dir>[,<dir>…] …   leaves on the high-frequency levels
    ch4-3e-tables.py summary <machine>=<dir>[,<dir>…] …  geometric means over c_leaf 2's
    ch4-3e-tables.py cases <machine>=<dir>[,<dir>…] …    every case, both orders
    ch4-3e-tables.py errors <machine>=<dir>[,<dir>…] …   the errors over ε and over c_leaf 2's

A level is high frequency where k a = Re(k) √3 r_l ≥ 2.8 (S6's range, the time rule's floor
`DIAGONAL_AUTO_MIN_KA`; the domain the unit cube, so r_l = 2^-(l+1)), whatever its M2L; the census
names the M2L that ran there (the column "M2L"'s `diagonal l:c,…`, else the table "Time rule"'s
"ran": `rcr`, `waves`, `batched`, or `no V` on a level without V pairs) and lists per run the
occupied leaves on such levels (the table "Leaves"), the W interactions into them (M2P at that
level, "Leaf operators") and the X interactions out of them (P2L at the next finer level), and
the phase bound (`Row::phase_bound`).
"""

import glob
import math
import os
import re
import sys
from collections import defaultdict

MIN_KA = 2.8  # DIAGONAL_AUTO_MIN_KA
REFERENCE = 2.0  # CH4.0d3's DEFAULT_C_LEAF, the reference of every ratio


def cells(line):
    return [c.strip() for c in line.strip().strip("|").split("|")]


def tables(path):
    """The Markdown tables of a file: {title: [rows as dicts]}; the first table is "runs"."""
    out = defaultdict(list)
    title = "runs"
    header = None
    with open(path) as f:
        for line in f:
            if line.startswith("## "):
                title = line[3:].strip()
                header = None
                continue
            if not line.startswith("|"):
                header = None if not line.strip() else header
                continue
            row = cells(line)
            if header is None:
                header = row
                continue
            if set(row[0]) <= set("-: "):
                continue
            out[title].append(dict(zip(header, row)))
    return out


def c_of(m2l):
    m = re.search(r"/c([0-9.]+)$", m2l)
    return float(m.group(1)) if m else None


def form_levels(m2l):
    m = re.search(r"diagonal ([0-9:,]+)", m2l)
    if not m:
        return {}
    return {int(l): int(c) for l, c in (p.split(":") for p in m.group(1).split(","))}


def by_level(text):
    """`3:0/240 4:1092/2124` -> {3: (occupied, leaves)}."""
    out = {}
    for part in text.split():
        level, counts = part.split(":")
        occ, leaves = counts.split("/")
        out[int(level)] = (int(occ), int(leaves))
    return out


def number(text):
    try:
        return float(text)
    except ValueError:
        return float("nan")


def runs_of(path):
    """Every run of a launch, in its order: a dict per run."""
    t = tables(path)
    runs = []
    for row in t["runs"]:
        c = c_of(row["M2L"])
        runs.append(
            {
                "c": c,
                "workload": row["Workload"],
                "n": int(row["N"]),
                "kd": row["kD"],
                "eps": row["degrees"].replace("ε = ", ""),
                "build": number(row["build s"]),
                "evaluate": number(row["evaluate s"]),
                "phi": number(row["φ"]),
                "grad": number(row["∇φ"]),
                "leaf level": int(row["leaf level"]),
                "leaves": int(row["leaves"]),
                "form": form_levels(row["M2L"]),
                "checksum": row["checksum"],
            }
        )
    leaves = t["Leaves"]
    assert len(leaves) == len(runs), path
    for run, row in zip(runs, leaves):
        assert float(row["c_leaf"]) == run["c"], path
        run["by level"] = by_level(row["occupied/leaves by level"])
        run["leaf size"] = int(row["leaf size"])
        run["phase"] = number(row["Re(k) r bound"])
        run["levels"] = int(row["levels"])
    # The table "Time rule": rows per level with V pairs, the runs in order.
    groups, last = [], None
    for row in t["Time rule"]:
        level = int(row["level"])
        if last is None or level <= last:
            groups.append([])
        groups[-1].append(row)
        last = level
    assert len(groups) == len(runs), path
    for run, rows in zip(runs, groups):
        run["ka"] = {int(r["level"]): number(r["k a"]) for r in rows}
        run["ran"] = {int(r["level"]): r["ran"] for r in rows}
    # "Leaf operators": W interactions (M2P into the level's leaves), X interactions (P2L into
    # the level's boxes, from leaves of the next coarser level).
    ops = defaultdict(dict)
    for row in t["Leaf operators"]:
        if row["level"] == "Σ":
            continue
        c = c_of(row["M2L"])
        ops[c][int(row["level"])] = (number(row["W interactions"]), number(row["X interactions"]))
    for run in runs:
        run["ops"] = ops.get(run["c"], {})
    return runs


def ka_of(run, level):
    """k a = Re(k) √3 r_l on the unit cube (kD = Re(k))."""
    return float(run["kd"]) * math.sqrt(3.0) * 0.5 ** (level + 1)


def m2l_on(run, level):
    if level in run["form"]:
        return f"form {run['form'][level]}"
    ran = run["ran"].get(level)
    if ran is None:
        return "no V"
    for prefix, name in (("rotation", "rcr"), ("planewave", "waves"), ("batched", "batched")):
        if ran.startswith(prefix):
            return name
    return ran


def high_levels(run):
    """The run's high-frequency levels with their M2L: {level: name}."""
    out = {}
    for level in range(run["levels"]):
        if ka_of(run, level) >= MIN_KA:
            if level in run["ka"]:
                assert abs(run["ka"][level] - ka_of(run, level)) < 0.01, (level, run["ka"])
            out[level] = m2l_on(run, level)
    return out


def census_row(run):
    """[occupied leaves on high levels, W into them, X out of them, the levels]."""
    occ = w = x = 0
    levels = []
    for level, kind in sorted(high_levels(run).items()):
        o = run["by level"].get(level, (0, 0))[0]
        if o == 0:
            continue
        occ += o
        w += run["ops"].get(level, (0, 0))[0]
        x += run["ops"].get(level + 1, (0, 0))[1]
        levels.append(f"{level} ({kind}, k a {ka_of(run, level):.1f}): {o}")
    return occ, w, x, levels


def cases(spec):
    """{case: {"fwd": runs, "rev": runs}} of `machine=dir[,dir]`."""
    machine, dirs = spec.split("=", 1)
    out = defaultdict(dict)
    for d in dirs.split(","):
        for path in sorted(glob.glob(os.path.join(d, "*.md"))):
            name = os.path.basename(path)[:-3]
            order, _, case = name.partition("-")
            if order not in ("fwd", "rev") or case == "smoke" or os.path.getsize(path) == 0:
                continue
            try:
                out[case][order] = runs_of(path)
            except (AssertionError, KeyError) as error:
                sys.exit(f"{path}: unreadable ({error!r})")
    return machine, out


def parse(case):
    """`w1-100-1e-3` -> ("w1", 100.0, "1e-3")."""
    workload, kd, eps = case.split("-", 2)
    return workload, float(kd), eps


def eps_of(case):
    return parse(case)[2]


def sort_key(case):
    workload, kd, eps = parse(case)
    return (workload, kd, -float(eps))


def by_c(runs):
    return {r["c"]: r for r in runs}


def gmean(values):
    values = [v for v in values if v > 0 and math.isfinite(v)]
    return math.exp(sum(map(math.log, values)) / len(values)) if values else float("nan")


def ratios(case_runs, measure):
    """{c: geometric mean over both orders of measure(c)/measure(2) in the same launch}."""
    out = {}
    for c in (1.0, 2.0, 4.0, 8.0):
        rs = []
        for order in ("fwd", "rev"):
            runs = by_c(case_runs.get(order, []))
            if c in runs and REFERENCE in runs:
                rs.append(measure(runs[c]) / measure(runs[REFERENCE]))
        out[c] = gmean(rs) if rs else float("nan")
    return out


def total(run):
    return run["build"] + run["evaluate"]


def fmt(x, digits=2):
    return "–" if not math.isfinite(x) else f"{x:.{digits}f}"


def census(specs):
    print("| machine | case | c_leaf | leaf levels (occupied) | high-frequency levels with leaves: "
          "(M2L, k a): occupied | occupied there | W into them | X out of them | Re(k) r bound |")
    print("| --- | --- | ---: | --- | --- | ---: | ---: | ---: | ---: |")
    totals = defaultdict(lambda: [0, 0])
    for spec in specs:
        machine, all_cases = cases(spec)
        for case in sorted(all_cases, key=sort_key):
            runs = by_c(all_cases[case].get("fwd") or all_cases[case].get("rev"))
            for c, run in sorted(runs.items()):
                occ, w, x, levels = census_row(run)
                eps = eps_of(case)
                totals[(machine, c, eps)][0] += 1
                totals[(machine, c, eps)][1] += occ > 0
                occupied = [l for l, (o, _) in sorted(run["by level"].items()) if o > 0]
                span = f"{occupied[0]}–{occupied[-1]}" if len(occupied) > 1 else str(occupied[0])
                print(
                    f"| {machine} | {case} | {c:g} | {span} | {'; '.join(levels) or 'none'} "
                    f"| {occ} | {w:.0f} | {x:.0f} | {run['phase']:.1f} |"
                )
    print()
    print("| machine | c_leaf | ε | trees | with leaves on a high-frequency level |")
    print("| --- | ---: | --- | ---: | ---: |")
    for (machine, c, eps), (n, h) in sorted(totals.items()):
        print(f"| {machine} | {c:g} | {eps} | {n} | {h} |")


def summary(specs):
    print("| machine | ε | cases | build + evaluate ÷ c_leaf 2: c_leaf 1 | 4 | 8 "
          "| evaluate ÷ c_leaf 2: c_leaf 1 | 4 | 8 | cases where 4 or 8 changes the tree |")
    print("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for spec in specs:
        machine, all_cases = cases(spec)
        groups = defaultdict(list)
        for case, runs in all_cases.items():
            eps = eps_of(case)
            groups[eps].append(runs)
            groups["all"].append(runs)
        for eps in sorted(groups, key=lambda e: (e == "all", -float(e) if e != "all" else 0)):
            group = groups[eps]
            t = [ratios(r, total) for r in group]
            e = [ratios(r, lambda run: run["evaluate"]) for r in group]
            changed = sum(
                any(
                    by_c(r.get("fwd") or r.get("rev"))[c]["checksum"]
                    != by_c(r.get("fwd") or r.get("rev"))[REFERENCE]["checksum"]
                    for c in (4.0, 8.0)
                )
                for r in group
            )
            g = lambda rs, c: fmt(gmean([x[c] for x in rs]), 3)
            print(
                f"| {machine} | {eps} | {len(group)} | {g(t, 1.0)} | {g(t, 4.0)} | {g(t, 8.0)} "
                f"| {g(e, 1.0)} | {g(e, 4.0)} | {g(e, 8.0)} | {changed} |"
            )


def per_case(specs):
    print("| machine | case | N | leaf size | leaf level at c_leaf 1, 2, 4, 8 | c_leaf 2: build s "
          "| evaluate s | build + evaluate ÷ c_leaf 2: 1 | 4 | 8 | evaluate ÷: 1 | 4 | 8 |")
    print("| --- | --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for spec in specs:
        machine, all_cases = cases(spec)
        for case in sorted(all_cases, key=sort_key):
            r = all_cases[case]
            runs = by_c(r.get("fwd") or r.get("rev"))
            ref = runs[REFERENCE]
            build = gmean([by_c(r[o])[REFERENCE]["build"] for o in r])
            evaluate = gmean([by_c(r[o])[REFERENCE]["evaluate"] for o in r])
            t = ratios(r, total)
            e = ratios(r, lambda run: run["evaluate"])
            levels = ", ".join(str(runs[c]["leaf level"]) for c in sorted(runs))
            print(
                f"| {machine} | {case} | {ref['n']} | {ref['leaf size']} | {levels} "
                f"| {build:.2f} | {evaluate:.3f} | {fmt(t[1.0])} | {fmt(t[4.0])} | {fmt(t[8.0])} "
                f"| {fmt(e[1.0])} | {fmt(e[4.0])} | {fmt(e[8.0])} |"
            )


def errors(specs):
    print("| machine | ε | runs | φ ÷ ε max | ∇φ ÷ ε max | φ ÷ c_leaf 2's: c_leaf 1 | 4 | 8 "
          "| ∇φ ÷ c_leaf 2's: 1 | 4 | 8 | every order and run bit for bit the same output per c_leaf |")
    print("| --- | --- | ---: | ---: | ---: | --- | --- | --- | --- | --- | --- | --- |")
    for spec in specs:
        machine, all_cases = cases(spec)
        groups = defaultdict(list)
        for case, r in all_cases.items():
            groups[eps_of(case)].append(r)
        for eps in sorted(groups, key=lambda e: -float(e)):
            e = float(eps)
            n = 0
            worst = [0.0, 0.0]
            rel = {c: [[], []] for c in (1.0, 4.0, 8.0)}
            same = True
            for r in groups[eps]:
                orders = [by_c(r[o]) for o in r]
                for runs in orders:
                    for c, run in runs.items():
                        n += 1
                        worst[0] = max(worst[0], run["phi"] / e)
                        worst[1] = max(worst[1], run["grad"] / e)
                for c in rel:
                    runs = orders[0]
                    if c in runs:
                        rel[c][0].append(runs[c]["phi"] / runs[REFERENCE]["phi"])
                        rel[c][1].append(runs[c]["grad"] / runs[REFERENCE]["grad"])
                if len(orders) == 2:
                    same &= all(
                        orders[0][c]["checksum"] == orders[1][c]["checksum"] for c in orders[0]
                    )
            span = lambda xs: f"{min(xs):.3f}–{max(xs):.3f}" if xs else "–"
            print(
                f"| {machine} | {eps} | {n} | {worst[0]:.2f} | {worst[1]:.2f} "
                f"| {span(rel[1.0][0])} | {span(rel[4.0][0])} | {span(rel[8.0][0])} "
                f"| {span(rel[1.0][1])} | {span(rel[4.0][1])} | {span(rel[8.0][1])} "
                f"| {'yes' if same else 'no'} |"
            )


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    mode, specs = sys.argv[1], sys.argv[2:]
    {"census": census, "summary": summary, "cases": per_case, "errors": errors}[mode](specs)


if __name__ == "__main__":
    main()
