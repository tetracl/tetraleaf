#!/usr/bin/env python3
"""tools/helmholtz/ch4-3b2c-tables.py: the tables of Helmholtz CH4.3b2c's report (the diagonal
form's time rule by default; docs/helmholtz/tasks/CH4.3b2c.md), from the launches of
`tools/helmholtz/ch4-3b2c-runs.sh`.

    python3 tools/helmholtz/ch4-3b2c-tables.py errors  <machine>=<dir>[:<prefix>] …
    python3 tools/helmholtz/ch4-3b2c-tables.py choices <machine>=<dir>[:<prefix>] …
    python3 tools/helmholtz/ch4-3b2c-tables.py summary <machine>=<dir>[:<prefix>] …
    python3 tools/helmholtz/ch4-3b2c-tables.py build   <machine>=<dir>[:<prefix>] …
    python3 tools/helmholtz/ch4-3b2c-tables.py repeat  <machine>=<dir>[:<prefix>] …
    python3 tools/helmholtz/ch4-3b2c-tables.py order   <machine>=<dir> …

<prefix> keeps the launches whose names start with it (`run`, `reverse`, `repeat`, `ranks`).
`errors`, `choices` and `summary` are CH4.3b2a's (`ch4-3b2a-tables.py`): per launch the time
rule's errors and evaluation over the defaults' (the diagonal form off); per level in the rule's
range the options timed at the level's degree in the same launch and the choice within 1.05× of
the fastest; their summary per ε. `build` prints per launch the build and its operator's parts,
the defaults' and the time rule's, and the evaluation. `repeat` prints per launch of
`--diagonal off,auto,off,auto` each run's evaluation and its levels' M2L seconds in launch order
(the spread: the two `off` runs; the place in the launch: each run against the one before).
`order` pairs each `run-` launch (off first) with its `reverse-` launch (auto first): the time
rule's evaluation over the defaults' in each order. No timing is asserted; the numbers go into
tetraleaf-validate/results/helmholtz-ch4-3b2c.md.
"""

import importlib.util
import math
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("b2a", HERE / "ch4-3b2a-tables.py")
b2a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b2a)


def launches(arg):
    """`machine=dir[:prefix]`: the machine and the launches' outputs, by name."""
    machine, where = arg.split("=", 1)
    prefix = ""
    if ":" in where:
        where, prefix = where.rsplit(":", 1)
    paths = sorted(p for p in pathlib.Path(where).glob("*.md") if p.stem.startswith(prefix))
    return machine, paths


b2a.launches = launches


def choices(path):
    """CH4.3b2a's `choices`, but a level that a forced run does not list (its plan converted the
    level's V pairs, or put it out of the rule's range) offers no option from that run; the
    level is left out where the machine run does not list it."""
    main, levels, ka = b2a.runs(path)
    if not {"defaults", "machine", "form", "rcr"} <= set(levels):
        return []
    out = []
    for l, (p, ran, seconds) in sorted(levels["defaults"].items()):
        if l not in ka or ka[l] < b2a.MIN_KA or l not in levels["machine"]:
            continue
        own = b2a.kind(ran)
        if own == "other":
            continue
        options = {"waves" if own == "waves" else "rcr": seconds}
        ran_f = "–"
        if l in levels["form"]:
            _, ran_f, s_f = levels["form"][l]
            if b2a.kind(ran_f) == "form":
                options["form"] = s_f
        if l in levels["rcr"]:
            _, ran_r, s_r = levels["rcr"][l]
            if b2a.kind(ran_r) == "rcr" and "rcr" not in options:
                options["rcr"] = s_r
        p_m, ran_m, s_m = levels["machine"][l]
        out.append((l, p_m, ka[l], options, own, b2a.kind(ran_m), s_m, ran_f))
    return out


b2a.choices = choices


def ordered(path):
    """The runs of a launch in its order: per run (setting, main row, evaluation s, {level: (p,
    ran, M2L s)}, build row)."""
    t = b2a.tables(path)
    mains = t.get("", [])
    builds = t.get("Build", [])
    per_level, evaluations, current = [], [], None
    for r in t.get("M2L per level", []):
        if r["level"] == "Σ":
            m = re.match(r"evaluate ([0-9.]+) s", r["ran"])
            evaluations.append(float(m.group(1)) if m else math.nan)
            per_level.append(current or {})
            current = None
            continue
        current = current or {}
        current[int(r["level"])] = (int(r["p"]), r["ran"], float(r["M2L s"]))
    out = []
    for i, main in enumerate(mains):
        out.append(
            (
                b2a.setting(main["M2L"]),
                main,
                evaluations[i] if i < len(evaluations) else math.nan,
                per_level[i] if i < len(per_level) else {},
                builds[i] if i < len(builds) else {},
            )
        )
    return out


def build_table(args):
    keys = ["build s", "operator s", "RCR tables s", "plane-wave tables s", "diagonal tables s"]
    print("| machine | launch | setting | " + " | ".join(keys) + " | evaluate s | build ÷ off "
          "| evaluate ÷ off |")
    print("| --- | --- | --- | " + " | ".join("---:" for _ in keys) + " | ---: | ---: | ---: |")
    for arg in args:
        machine, paths = launches(arg)
        for path in paths:
            runs = ordered(path)
            off = next((r for r in runs if r[0] == "defaults"), None)
            for name, _, evaluation, _, build in runs:
                if name not in ("defaults", "machine") or not build:
                    continue
                cells = " | ".join(build.get(k, "–") for k in keys)
                ratio_b = float(build["build s"]) / float(off[4]["build s"]) if off else math.nan
                ratio_e = evaluation / off[2] if off else math.nan
                label = "off" if name == "defaults" else "auto"
                print(f"| {machine} | {path.stem} | {label} | {cells} | {evaluation:.4f} "
                      f"| {ratio_b:.3f} | {ratio_e:.3f} |")


def repeat_table(args):
    print("| machine | launch | run | setting | evaluate s | ÷ first off | M2L per level (level:ran "
          "s) |")
    print("| --- | --- | ---: | --- | ---: | ---: | --- |")
    for arg in args:
        machine, paths = launches(arg)
        for path in paths:
            runs = ordered(path)
            first = next((r[2] for r in runs if r[0] == "defaults"), math.nan)
            for i, (name, _, evaluation, levels, _) in enumerate(runs, 1):
                label = "off" if name == "defaults" else "auto"
                cells = ", ".join(
                    f"{l}:{b2a.kind(ran)} {s:.3f}" for l, (_, ran, s) in sorted(levels.items())
                )
                print(f"| {machine} | {path.stem} | {i} | {label} | {evaluation:.4f} "
                      f"| {evaluation / first:.3f} | {cells} |")


def order_table(args):
    print("| machine | case | auto ÷ off, off first | auto ÷ off, auto first | mean |")
    print("| --- | --- | ---: | ---: | ---: |")
    for arg in args:
        machine, where = arg.split("=", 1)
        for path in sorted(pathlib.Path(where).glob("run-*.md")):
            case = path.stem[len("run-"):]
            other = path.with_name(f"reverse-{case}.md")
            if not other.exists():
                continue
            ratios = []
            for p in (path, other):
                runs = {r[0]: r[2] for r in ordered(p) if r[0] in ("defaults", "machine")}
                ratios.append(runs["machine"] / runs["defaults"])
            mean = math.sqrt(ratios[0] * ratios[1])
            print(f"| {machine} | {case} | {ratios[0]:.3f} | {ratios[1]:.3f} | {mean:.3f} |")


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    {
        "errors": b2a.errors_table,
        "choices": b2a.choices_table,
        "summary": b2a.summary_table,
        "build": build_table,
        "repeat": repeat_table,
        "order": order_table,
    }[sys.argv[1]](sys.argv[2:])


if __name__ == "__main__":
    main()
