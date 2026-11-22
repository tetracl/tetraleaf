#!/usr/bin/env python3
"""The report tables of Helmholtz CH4.4a (docs/helmholtz/tasks/CH4.4a.md) from the cases'
Markdown that `tools/helmholtz/ch4-4a-runs.sh` leaves in its output directories
(helmholtz-<case>.md, written by `helmholtz_fmm3d compare`):

    tools/helmholtz/ch4-4a-tables.py <table> <dir>…

Tables:
  accuracy  per case and ε: ours and FMM3D's (each variant) relative ℓ2 errors of φ and ∇φ at the
            samples, ours over ε, the all-target errors and the worst box over ε (ours gated at
            1 and 100); a last line with the largest ratios and the misses
  time      per case and ε: our build, evaluation and both, FMM3D's call per variant, and ours
            (build + evaluation, and evaluation alone) over each variant's
  memory    per case and ε: our resident growth summed over the ranks and the largest rank's,
            FMM3D's peak per variant (GB)
  build     per case and ε: our build and its larger parts (seconds)
Python 3.9 standard library only.
"""

import glob
import os
import re
import sys

WORKLOAD_ORDER = {"W1": 0, "W2": 1, "W3": 2}
DENSITY_ORDER = {"random": 0, "constant": 1, "planewave": 2}


def cells(line):
    return [c.strip() for c in line.strip().strip("|").split("|")]


def parse(path):
    """One case file: its density, the main table's rows and the build table's rows."""
    text = open(path, encoding="utf-8").read()
    density = "random"
    m = re.search(r"charge vector \((\w+)\)", text)
    if m:
        density = m.group(1)
    rows, build, header, build_header = [], [], None, None
    section = "main"
    for line in text.splitlines():
        if line.startswith("## Build"):
            section = "build"
            continue
        if not line.startswith("|") or line.startswith("| ---"):
            continue
        c = cells(line)
        if section == "main":
            if c[0] == "Workload":
                header = c
                continue
            if header and len(c) == len(header):
                rows.append(dict(zip(header, c)))
        else:
            if c[0] == "ε":
                build_header = c
                continue
            if build_header and len(c) == len(build_header):
                build.append(dict(zip(build_header, c)))
    return density, rows, build


def num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def cases(dirs):
    out = []
    for d in dirs:
        for path in sorted(glob.glob(os.path.join(d, "helmholtz-*.md"))):
            density, rows, build = parse(path)
            if not rows:
                continue
            first = rows[0]
            key = (
                WORKLOAD_ORDER.get(first["Workload"], 9),
                DENSITY_ORDER.get(density, 9),
                int(first["N"]) if first["Workload"] == "W3" else 0,
                float(first["kD"].split()[0]),
            )
            out.append((key, os.path.basename(path), density, rows, build))
    out.sort(key=lambda c: c[0])
    return out


def label(rows, density):
    r = rows[0]
    name = r["Workload"]
    if name == "W3":
        name += f" (N = {int(r['N']):,})"
    if density != "random":
        name += f", {density}"
    return f"{name}, kD = {r['kD']}"


def eps_of(r):
    return float(r["ε"])


def worst(r):
    """The worst box over ε, [φ, ∇φ], from 'a / b'."""
    parts = [num(p) for p in r["worst box φ / ε, ∇φ / ε"].split("/")]
    return parts if len(parts) == 2 and None not in parts else None


def accuracy(all_cases):
    print("| Case | ε | ours φ, ∇φ | ours ÷ ε | ours all φ, ∇φ | ours worst box ÷ ε | FMM3D φ, ∇φ (by variant) | FMM3D worst box ÷ ε |")
    print("| --- | --- | --- | --- | --- | --- | --- | --- |")
    largest, largest_box, misses = 0.0, 0.0, []
    for _, _, density, rows, _ in all_cases:
        name = label(rows, density)
        by_eps = {}
        for r in rows:
            by_eps.setdefault(r["ε"], []).append(r)
        for e, group in by_eps.items():
            eps = float(e)
            ours = [r for r in group if r["code"].startswith("ours")]
            theirs = [r for r in group if r["code"].startswith("FMM3D")]
            if not ours:
                continue
            o = ours[0]
            phi, grad = num(o["φ"]), num(o["∇φ"])
            ratio = max(phi, grad) / eps
            largest = max(largest, ratio)
            box = worst(o)
            if box:
                largest_box = max(largest_box, max(box))
            if ratio > 1.0 or (box and max(box) > 100.0):
                misses.append(f"{name} at {e}")
            f3 = "; ".join(f"{r['code'][6:]} {r['φ']}, {r['∇φ']}" for r in theirs) or "–"
            fbox = "; ".join(
                f"{r['code'][6:]} {r['worst box φ / ε, ∇φ / ε']}"
                for r in theirs
                if r["worst box φ / ε, ∇φ / ε"] != "–"
            ) or "–"
            print(
                f"| {name} | {e} | {o['φ']}, {o['∇φ']} | {phi / eps:.2f}, {grad / eps:.2f} | "
                f"{o['all φ']}, {o['all ∇φ']} | {o['worst box φ / ε, ∇φ / ε']} | {f3} | {fbox} |"
            )
            name = ""
    print()
    print(
        f"Largest of ours over ε: {largest:.2f}; largest worst box over ε: {largest_box:.1f}; "
        f"misses (over ε, or a worst box above 100 ε): {', '.join(misses) or 'none'}."
    )


def time(all_cases):
    variants = []
    for _, _, _, rows, _ in all_cases:
        for r in rows:
            v = r["code"][6:] if r["code"].startswith("FMM3D") else None
            if v and v not in variants:
                variants.append(v)
    head = " | ".join(f"FMM3D {v} s" for v in variants)
    ratios = " | ".join(f"ours ÷ {v} (total, run)" for v in variants)
    print(f"| Case | N | ε | ours build s | ours run s | ours total s | {head} | {ratios} |")
    print("| --- | ---: | --- | ---: | ---: | ---: |" + " ---: |" * (2 * len(variants)))
    for _, _, density, rows, _ in all_cases:
        name = label(rows, density)
        by_eps = {}
        for r in rows:
            by_eps.setdefault(r["ε"], []).append(r)
        for e, group in by_eps.items():
            ours = [r for r in group if r["code"].startswith("ours")]
            if not ours:
                continue
            o = ours[0]
            build, run = num(o["build s"]), num(o["run s"])
            total = build + run
            secs, rats = [], []
            for v in variants:
                t = [num(r["run s"]) for r in group if r["code"] == f"FMM3D {v}"]
                t = t[0] if t else None
                secs.append(f"{t:.2f}" if t else "–")
                rats.append(f"{total / t:.2f}, {run / t:.2f}" if t else "–")
            print(
                f"| {name} | {int(o['N']):,} | {e} | {build:.2f} | {run:.2f} | {total:.2f} | "
                + " | ".join(secs) + " | " + " | ".join(rats) + " |"
            )
            name = ""


def memory(all_cases):
    variants = []
    for _, _, _, rows, _ in all_cases:
        for r in rows:
            v = r["code"][6:] if r["code"].startswith("FMM3D") else None
            if v and v not in variants:
                variants.append(v)
    head = " | ".join(f"FMM3D {v} GB" for v in variants)
    print(f"| Case | N | ε | ours GB (largest rank) | {head} |")
    print("| --- | ---: | --- | ---: |" + " ---: |" * len(variants))
    for _, _, density, rows, _ in all_cases:
        name = label(rows, density)
        by_eps = {}
        for r in rows:
            by_eps.setdefault(r["ε"], []).append(r)
        for e, group in by_eps.items():
            ours = [r for r in group if r["code"].startswith("ours")]
            if not ours:
                continue
            o = ours[0]
            mem = []
            for v in variants:
                t = [r.get("memory GB", "–") for r in group if r["code"] == f"FMM3D {v}"]
                mem.append(t[0] if t else "–")
            print(
                f"| {name} | {int(o['N']):,} | {e} | {o.get('memory GB', '–')} | "
                + " | ".join(mem) + " |"
            )
            name = ""


BUILD_SHOWN = [
    "octree",
    "plan",
    "operator",
    "RCR tables",
    "matrices",
    "plane-wave tables",
    "diagonal tables",
    "scratch",
    "replan",
    "evaluator",
]


def build(all_cases):
    print("| Case | ε | build s | " + " | ".join(BUILD_SHOWN) + " | evaluation s |")
    print("| --- | --- | ---: |" + " ---: |" * (len(BUILD_SHOWN) + 1))
    for _, _, density, rows, builds in all_cases:
        name = label(rows, density)
        for b in builds:
            print(
                f"| {name} | {b['ε']} | {b['build']} | "
                + " | ".join(b.get(p, "–") for p in BUILD_SHOWN)
                + f" | {b['evaluation']} |"
            )
            name = ""


def main():
    if len(sys.argv) < 3 or sys.argv[1] not in ("accuracy", "time", "memory", "build"):
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    all_cases = cases(sys.argv[2:])
    {"accuracy": accuracy, "time": time, "memory": memory, "build": build}[sys.argv[1]](
        all_cases
    )


if __name__ == "__main__":
    main()
