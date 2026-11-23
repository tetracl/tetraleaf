#!/usr/bin/env python3
"""tools/helmholtz/ch4-5-tables.py: the eps lists and the report's tables of Helmholtz CH4.5
(docs/helmholtz/tasks/CH4.5.md) from the output directories of tools/helmholtz/ch4-5-runs.sh.

    python3 tools/helmholtz/ch4-5-tables.py eps <grid dir> <ours dir>... [--default-1000]
    python3 tools/helmholtz/ch4-5-tables.py bind <grid dir> <ours dir>...
    python3 tools/helmholtz/ch4-5-tables.py matches <grid dir> <ours dir>
    python3 tools/helmholtz/ch4-5-tables.py compare <grid dir> <dir>...
    python3 tools/helmholtz/ch4-5-tables.py summary <grid dir> <dir>...
    python3 tools/helmholtz/ch4-5-tables.py build <dir>...
    python3 tools/helmholtz/ch4-5-tables.py memory <grid dir> <dir>...
    python3 tools/helmholtz/ch4-5-tables.py cleaf <dir>...
    python3 tools/helmholtz/ch4-5-tables.py binding <dir>...
    python3 tools/helmholtz/ch4-5-tables.py checksums <dir>...
    python3 tools/helmholtz/ch4-5-tables.py selftest <CH4.0e's log dir>

Several directories of one machine are read as one. The grid (Kathleen's `grid` stage, FMM3D's
errors over eps; FMM3D's errors are the same on both machines to the digits shown, CH4.0e)
decides the matches; each machine's own `ours` decides our errors.

The two matches, per case, ε and FMM3D variant, over the grid eps that left an output:
- **closest** (the acceptance's, as C6.8): the eps minimising d = max over φ and ∇φ of
  |log₁₀(FMM3D's error / ours)|, the faster grid call on a tie; d is printed as the factor 10^d;
- **one-sided** (CH4.0e's): the fastest eps whose φ and ∇φ errors are both at most ours; where
  none is, the most accurate (the smallest worse of the two ratios), marked `>`.

`eps` prints the `time` stage's list (`<case> <variant> <eps,…>`: per case and variant the union
over ε and the machines of both matches and ε itself; with `--default-1000`, locust's, also
`default` at kD = 1000 at `fastker`'s eps); `bind` the `bind` stage's (a few large cases at the
closest eps of 10⁻⁶, every variant). `matches`: the matches of one machine. `compare`: per case
and ε ours at its best split (the smallest build plus one evaluation over the splits and launches
a and b; the evaluation is the second of the build's two) and at 1 × N, beside FMM3D at the
closest eps (its call the smaller of the `time` launches a and b, its tree apart) for each
variant, and the ratios. `summary`: the geometric means of the ratios by ε, rule, split and N; setup apart only where
FMM3D's tree, timed once apart, is below its call.
`build`: the build by part at 1 × N and at the best split, the shared steps (octree, routes, plan;
CH4.0k) named. `memory`: ours (largest rank and node) against FMM3D's peak at the matches.
`cleaf`: c_leaf 4 over 2 per split. `binding`: FMM3D bound over unbound. `checksums`: the
agreement of our checksums over splits and launches. `selftest`: the one-sided rule on CH4.0e's
three hand-checked rows. Python 3.9 standard library only.
"""

import glob
import math
import os
import re
import sys
from collections import defaultdict

EPS = ["1e-3", "1e-6", "1e-8"]
SHARED = ["octree", "routes", "plan"]
BIND_CASES = ["w1-kd250", "w1-kd500", "w1-kd1000", "w3-kd100", "w3-n1000000-kd500"]


def cells(line):
    return [c.strip() for c in line.strip().strip("|").split("|")]


def tables(path):
    """{section heading: [row dicts]}, the first table under each heading."""
    out, section, header = {}, None, None
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                section, header = line.lstrip("#").strip(), None
                continue
            if not line.startswith("|"):
                header = None if not line.strip() else header
                continue
            row = cells(line)
            if header is None:
                header = row
                out.setdefault(section, [])
                continue
            if set("".join(row)) <= set("-: "):
                continue
            out[section].append(dict(zip(header, row)))
    return out


def num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return float("nan")


def ok(x):
    return x is not None and not (isinstance(x, float) and math.isnan(x))


def gmean(xs):
    xs = [x for x in xs if ok(x) and x > 0]
    return math.exp(sum(map(math.log, xs)) / len(xs)) if xs else float("nan")


def apart(evaluate, t):
    """Our evaluation over FMM3D's call less its tree; None where the tree, timed once apart,
    is not below the call (the smallest cases, a call of milliseconds)."""
    rest = t["s"] - t["tree"]
    return evaluate / rest if ok(rest) and rest > 0 else None


def f(x, d=3):
    return f"{x:.{d}f}" if ok(x) else "–"


def e2(x):
    return f"{x:.2e}" if ok(x) else "–"


def kd_of(name):
    return float(re.search(r"-kd(\d+)", name).group(1))


def order(name):
    w = name.split("-kd")[0]
    dens = name.split("-kd")[1].split("-", 1)[1] if "-" in name.split("-kd")[1] else ""
    rank = {"w1": 0, "w2": 1, "w3": 2, "w3-n1000000": 3}.get(w, 9)
    return ({"": 0, "planewave": 1, "constant": 2}.get(dens, 9), rank, kd_of(name))


# --- our runs ---------------------------------------------------------------------------------


def our_runs(dirs, prefix="ours"):
    """{(name, eps): {split: [run dicts]}} from <prefix>-<name>-<eps>-<split>-<a|b>.md."""
    runs = defaultdict(lambda: defaultdict(list))
    pat = re.compile(re.escape(prefix) + r"-(.+)-(1e-\d+)-(\d+x\d+)-([ab])$")
    for d in dirs:
        for path in sorted(glob.glob(os.path.join(d, prefix + "*-*.md"))):
            base = os.path.basename(path)[:-3]
            m = pat.match(base)
            if not m:
                continue
            name, e, split, launch = m.groups()
            t = tables(path)
            main = next((rows for s, rows in t.items() if s and s.startswith("Helmholtz FMM beside")), [])
            build = t.get("Build", [])
            if not main or not build:
                continue
            r, b = main[0], build[0]
            run = {
                "launch": launch,
                "n": int(r["N"]),
                "build": num(b["build"]),
                "evaluate": num(b["evaluation"]),
                "first": num(b.get("first evaluation")),
                "phi": num(r["φ"]),
                "grad": num(r["∇φ"]),
                "checksum": r["checksum"],
                "p": r["p by level"],
                "memory": r["memory GB"],
                "parts": {k: num(v) for k, v in b.items() if k not in ("ε",)},
            }
            run["total"] = run["build"] + run["evaluate"]
            runs[(name, e)][split].append(run)
    return runs


def one_rank(per_split):
    return next((s for s in per_split if s.startswith("1x")), None)


def best(per_split, splits=None):
    """(split, run): the smallest total over the splits (all, or those named) and launches."""
    choice = None
    for split, rs in per_split.items():
        if splits is not None and split not in splits:
            continue
        for r in rs:
            if choice is None or r["total"] < choice[1]["total"]:
                choice = (split, r)
    return choice


def our_errors(runs, name, e):
    per = runs.get((name, e))
    if not per:
        return None
    r = next(iter(per.values()))[0]
    return r["phi"], r["grad"], r["n"]


# --- FMM3D ------------------------------------------------------------------------------------


def errors_tables(dirs, stage):
    """{name: [row dicts of <stage>-<name>.md]}."""
    out = {}
    for d in dirs:
        for path in glob.glob(os.path.join(d, stage + "-*.md")):
            name = os.path.basename(path)[len(stage) + 1 : -3]
            t = tables(path)
            rows = [r for rows in t.values() for r in rows if "run" in r]
            out.setdefault(name, []).extend(rows)
    return out


def grid(dirs):
    """{(name, variant): {eps: {"s", "phi", "grad", "tree", "memory"}}} from the grid's rows that
    left an output."""
    out = defaultdict(dict)
    for name, rows in errors_tables(dirs, "grid").items():
        for r in rows:
            if r["s"].startswith("failed"):
                continue
            out[(name, r["run"])][r["eps"]] = {
                "s": num(r["s"]),
                "phi": num(r["φ"]),
                "grad": num(r["∇φ"]),
                "tree": num(r.get("tree s")),
                "memory": num(r.get("memory GB")),
            }
    return out


def failed(dirs):
    """{(name, variant): {eps: memory GB}} of the grid's failed rows."""
    out = defaultdict(dict)
    for name, rows in errors_tables(dirs, "grid").items():
        for r in rows:
            if r["s"].startswith("failed"):
                out[(name, r["run"])][r["eps"]] = r.get("memory GB", "–")
    return out


def timed(dirs, stage="time"):
    """{(name, variant, binding): {eps: {"s": the smaller of launches a and b, "tree", "memory"
    (the larger), "phi", "grad", "launches"}}}."""
    out = defaultdict(dict)
    for name, rows in errors_tables(dirs, stage).items():
        for r in rows:
            m = re.match(r"(default|fastker)-(free|bound)-([ab])$", r["run"])
            if not m or r["s"].startswith("failed"):
                continue
            key = (name, m.group(1), m.group(2))
            cell = out[key].setdefault(r["eps"], {"s": float("inf"), "tree": float("inf"), "memory": 0.0, "launches": 0})
            cell["s"] = min(cell["s"], num(r["s"]))
            cell["tree"] = min(cell["tree"], num(r.get("tree s")))
            cell["memory"] = max(cell["memory"], num(r.get("memory GB")))
            cell["phi"], cell["grad"] = num(r["φ"]), num(r["∇φ"])
            cell["launches"] += 1
    return out


def distance(c, phi, grad):
    return max(abs(math.log10(c["phi"] / phi)), abs(math.log10(c["grad"] / grad)))


def closest(cells_, phi, grad):
    """(eps, cell, d) of the closest match, or None."""
    measured = [(distance(c, phi, grad), c["s"], e, c) for e, c in cells_.items() if c["phi"] > 0 and c["grad"] > 0]
    if not measured or phi <= 0 or grad <= 0:
        return None
    d, _, e, c = min(measured)
    return e, c, d


def one_sided(cells_, phi, grad):
    """(eps, cell, exact) of the one-sided match, or None."""
    good = [(c["s"], e, c) for e, c in cells_.items() if c["phi"] <= phi and c["grad"] <= grad]
    if good:
        _, e, c = min(good)
        return e, c, True
    measured = [(max(c["phi"] / phi, c["grad"] / grad), e, c) for e, c in cells_.items() if phi > 0 and grad > 0]
    if not measured:
        return None
    _, e, c = min(measured)
    return e, c, False


def variants_of(fm, name):
    return sorted({v for (n, v) in fm if n == name})


# --- commands ---------------------------------------------------------------------------------


def cmd_eps(args):
    default_1000 = "--default-1000" in args
    args = [a for a in args if a != "--default-1000"]
    fm = grid([args[0]])
    machines = [our_runs([d]) for d in args[1:]]
    names = sorted({n for (n, _) in fm}, key=order)
    for name in names:
        for v in variants_of(fm, name):
            chosen = set(EPS)
            for runs in machines:
                for e in EPS:
                    o = our_errors(runs, name, e)
                    if o is None:
                        continue
                    for rule in (closest, one_sided):
                        m = rule(fm[(name, v)], o[0], o[1])
                        if m:
                            chosen.add(m[0])
            listed = sorted(chosen, key=float, reverse=True)
            print(f"{name} {v} {','.join(listed)}")
            if default_1000 and kd_of(name) == 1000 and v == "fastker":
                print(f"{name} default {','.join(listed)}")


def cmd_bind(args):
    fm = grid([args[0]])
    runs = our_runs(args[1:])
    for name in BIND_CASES:
        o = our_errors(runs, name, "1e-6")
        for v in variants_of(fm, name):
            m = closest(fm[(name, v)], o[0], o[1]) if o else None
            if m:
                print(f"{name} {v} {m[0]}")


def cmd_matches(args):
    fm = grid([args[0]])
    runs = our_runs(args[1:])
    print("| Case | N | ε | ours φ, ∇φ | variant | closest eps: φ, ∇φ (factor) | one-sided eps: φ, ∇φ |")
    print("| --- | ---: | --- | --- | --- | --- | --- |")
    for name in sorted({n for (n, _) in runs}, key=order):
        for e in EPS:
            o = our_errors(runs, name, e)
            if o is None:
                continue
            for v in variants_of(fm, name):
                c = closest(fm[(name, v)], o[0], o[1])
                s = one_sided(fm[(name, v)], o[0], o[1])
                cc = f"{c[0]}: {e2(c[1]['phi'])}, {e2(c[1]['grad'])} ({10 ** c[2]:.1f})" if c else "–"
                sc = f"{'' if s[2] else '>'}{s[0]}: {e2(s[1]['phi'])}, {e2(s[1]['grad'])}" if s else "–"
                print(f"| {name} | {o[2]} | {e} | {e2(o[0])}, {e2(o[1])} | {v} | {cc} | {sc} |")


def comparison(grid_dir, dirs):
    """Per (name, ε): ours (best split, 1 × N) and per variant FMM3D at both matches."""
    fm = grid([grid_dir])
    runs = our_runs(dirs)
    tm = timed(dirs)
    out = []
    for name in sorted({n for (n, _) in runs}, key=order):
        for e in EPS:
            per = runs.get((name, e))
            if not per:
                continue
            bs, b = best(per)
            one = one_rank(per)
            o1 = best(per, {one}) if one else None
            row = {"name": name, "eps": e, "n": b["n"], "best": (bs, b), "one": o1, "fmm3d": {}}
            for v in variants_of(fm, name):
                for rule, fn in (("closest", closest), ("one-sided", one_sided)):
                    m = fn(fm[(name, v)], b["phi"], b["grad"])
                    if not m:
                        continue
                    t = tm.get((name, v, "free"), {}).get(m[0])
                    row["fmm3d"][(v, rule)] = {"eps": m[0], "grid": m[1], "how": m[2], "time": t}
            out.append(row)
    return out


def fastest(row, rule):
    """The variant with the faster timed call at the rule's match, with its cell."""
    best_ = None
    for (v, r), c in row["fmm3d"].items():
        if r != rule or not c["time"]:
            continue
        if best_ is None or c["time"]["s"] < best_[1]["time"]["s"]:
            best_ = (v, c)
    return best_


def cmd_compare(args):
    rows = comparison(args[0], args[1:])
    print(
        "| Case | N | ε | best split | build s | evaluate s | ours s | 1 × N: build s, evaluate s, s "
        "| FMM3D (closest): variant, eps (factor), call s, tree s | ours ÷ FMM3D: best, 1 × N "
        "| evaluate ÷ (call − tree) | FMM3D (one-sided): variant, eps, call s | ÷: best, 1 × N |"
    )
    print("| --- | ---: | --- | --- | ---: | ---: | ---: | --- | --- | --- | ---: | --- | --- |")
    for r in rows:
        bs, b = r["best"]
        one = r["one"]
        ones = f"{f(one[1]['build'])}, {f(one[1]['evaluate'])}, {f(one[1]['total'])}" if one else "–"
        cells_ = [r["name"], str(r["n"]), r["eps"], bs, f(b["build"]), f(b["evaluate"]), f(b["total"]), ones]
        c = fastest(r, "closest")
        if c:
            v, m = c
            t = m["time"]
            cells_ += [
                f"{v}, {m['eps']} ({10 ** m['how']:.1f}), {f(t['s'])}, {f(t['tree'])}",
                f"{f(b['total'] / t['s'], 2)}, {f(one[1]['total'] / t['s'], 2) if one else '–'}",
                f(apart(b["evaluate"], t), 2),
            ]
        else:
            cells_ += ["–", "–", "–"]
        c = fastest(r, "one-sided")
        if c:
            v, m = c
            t = m["time"]
            cells_ += [
                f"{v}, {'' if m['how'] else '>'}{m['eps']}, {f(t['s'])}",
                f"{f(b['total'] / t['s'], 2)}, {f(one[1]['total'] / t['s'], 2) if one else '–'}",
            ]
        else:
            cells_ += ["–", "–"]
        print("| " + " | ".join(cells_) + " |")


def cmd_summary(args):
    rows = comparison(args[0], args[1:])
    groups = [("all", lambda r: True), ("N ≥ 10⁵", lambda r: r["n"] >= 100000), ("N ≥ 10⁶", lambda r: r["n"] >= 1000000)]
    print("| rule | cases | ε | n | best split ÷ FMM3D: gmean (range) | 1 × N ÷ FMM3D: gmean (range) | evaluate ÷ (call − tree): gmean (range) [cases, where fewer: the tree below the call] | ours slower (best split) |")
    print("| --- | --- | --- | ---: | --- | --- | --- | --- |")
    for rule in ("closest", "one-sided"):
        for label, keep in groups:
            for e in EPS:
                xs, ys, zs, slower = [], [], [], []
                for r in rows:
                    if r["eps"] != e or not keep(r):
                        continue
                    c = fastest(r, rule)
                    if not c:
                        continue
                    t = c[1]["time"]
                    _, b = r["best"]
                    xs.append(b["total"] / t["s"])
                    if r["one"]:
                        ys.append(r["one"][1]["total"] / t["s"])
                    z = apart(b["evaluate"], t)
                    if z is not None:
                        zs.append(z)
                    if xs[-1] > 1:
                        slower.append(r["name"])
                if not xs:
                    continue
                rng = lambda v: f"{gmean(v):.2f} ({min(v):.2f}–{max(v):.2f})" if v else "–"
                zn = f" [{len(zs)}]" if len(zs) != len(xs) else ""
                print(f"| {rule} | {label} | {e} | {len(xs)} | {rng(xs)} | {rng(ys)} | {rng(zs)}{zn} | {', '.join(slower) or '–'} |")


def cmd_build(args):
    runs = our_runs(args)
    print("| Case | ε | split | build s | " + " | ".join(SHARED) + " | shared share | operator | evaluate s |")
    print("| --- | --- | --- | ---: | " + " | ".join("---:" for _ in SHARED) + " | ---: | ---: | ---: |")
    for (name, e) in sorted(runs, key=lambda k: (order(k[0]), k[1])):
        per = runs[(name, e)]
        if per and next(iter(per.values()))[0]["n"] < 100000:
            continue
        one = one_rank(per)
        for label, choice in (("1 × N", best(per, {one}) if one else None), ("best", best(per))):
            if not choice:
                continue
            split, r = choice
            shared = sum(r["parts"].get(k, 0.0) for k in SHARED)
            print(
                f"| {name} | {e} | {label} {split} | {f(r['build'])} | "
                + " | ".join(f(r["parts"].get(k)) for k in SHARED)
                + f" | {f(shared / r['build'], 2) if r['build'] > 0 else '–'} | {f(r['parts'].get('operator'))} | {f(r['evaluate'])} |"
            )


def cmd_memory(args):
    rows = comparison(args[0], args[1:])
    print("| Case | N | ε | ours best split: node GB (largest rank) | ours 1 × N GB | FMM3D closest: variant eps GB | FMM3D one-sided: variant eps GB |")
    print("| --- | ---: | --- | --- | --- | --- | --- |")
    for r in rows:
        if r["n"] < 100000:
            continue
        _, b = r["best"]
        cells_ = [r["name"], str(r["n"]), r["eps"], b["memory"], r["one"][1]["memory"] if r["one"] else "–"]
        for rule in ("closest", "one-sided"):
            c = fastest(r, rule)
            cells_.append(f"{c[0]} {c[1]['eps']} {f(c[1]['time']['memory'], 2)}" if c else "–")
        print("| " + " | ".join(cells_) + " |")


def cmd_cleaf(args):
    two, four = our_runs(args, "cleaf2"), our_runs(args, "cleaf4")
    print("| ε | split | c_leaf 2: build, evaluate, s | c_leaf 4: build, evaluate, s | 4 ÷ 2 | φ, ∇φ (4) |")
    print("| --- | --- | --- | --- | ---: | --- |")
    for key in sorted(four, key=lambda k: k[1]):
        for split in sorted(four[key], key=lambda s: (int(s.split("x")[0]), s)):
            a = best(two.get(key, {}), {split})
            b = best(four[key], {split})
            if not a or not b:
                continue
            fa, fb = a[1], b[1]
            print(
                f"| {key[1]} | {split} | {f(fa['build'])}, {f(fa['evaluate'])}, {f(fa['total'])} | "
                f"{f(fb['build'])}, {f(fb['evaluate'])}, {f(fb['total'])} | {f(fb['total'] / fa['total'], 2)} | "
                f"{e2(fb['phi'])}, {e2(fb['grad'])} |"
            )


def cmd_binding(args):
    tm = timed(args, "bind")
    print("| Case | variant | eps | unbound s (launches) | bound s (launches) | bound ÷ unbound |")
    print("| --- | --- | --- | ---: | ---: | ---: |")
    for (name, v, bnd) in sorted(tm, key=lambda k: (order(k[0]), k[1])):
        if bnd != "free":
            continue
        for e, c in tm[(name, v, "free")].items():
            d = tm.get((name, v, "bound"), {}).get(e)
            print(
                f"| {name} | {v} | {e} | {f(c['s'])} ({c['launches']}) | "
                f"{f(d['s']) if d else '–'} ({d['launches'] if d else 0}) | {f(d['s'] / c['s'], 2) if d else '–'} |"
            )


def cmd_checksums(args):
    for prefix in ("ours", "cleaf2", "cleaf4"):
        runs = our_runs(args, prefix)
        bad = []
        for key, per in runs.items():
            sums = {r["checksum"] for rs in per.values() for r in rs}
            if len(sums) != 1:
                bad.append((key, sorted(sums)))
        print(f"{prefix}: {len(runs)} case-ε pairs, {sum(len(rs) for p in runs.values() for rs in p.values())} launches; "
              f"checksums differing: {len(bad)}")
        for key, sums in bad:
            print(f"  {key}: {sums}")


def cmd_selftest(args):
    """The one-sided rule against CH4.0e's three hand-checked rows (its report §1; our errors
    from its §2.1 table, FMM3D's grid from launch a's `fmm3d-<case>.md`)."""
    d = args[0]
    checks = [
        ("w1-500", (1.17e-05, 2.10e-05), "default", "3e-4"),
        ("w3-30", (1.39e-07, 3.53e-07), "default", "3e-4"),
        ("w4-plummer-100", (3.05e-10, 1.39e-10), "default", "3e-7"),
    ]
    failures = 0
    for case, (phi, grad), v, want in checks:
        path = os.path.join(d, f"fmm3d-{case}.md")
        if not os.path.exists(path):
            print(f"{case}: no {path}")
            failures += 1
            continue
        cells_ = {}
        for rows in tables(path).values():
            for r in rows:
                if r.get("run") == v:
                    cells_[r["eps"]] = {"s": num(r["s"]), "phi": num(r["φ"]), "grad": num(r["∇φ"])}
        got = one_sided(cells_, phi, grad)
        status = "ok" if got and got[0] == want else "FAIL"
        failures += status != "ok"
        print(f"{case}: one-sided {got[0] if got else None}, CH4.0e {want}: {status}")
    sys.exit(1 if failures else 0)


COMMANDS = {
    "eps": cmd_eps,
    "bind": cmd_bind,
    "matches": cmd_matches,
    "compare": cmd_compare,
    "summary": cmd_summary,
    "build": cmd_build,
    "memory": cmd_memory,
    "cleaf": cmd_cleaf,
    "binding": cmd_binding,
    "checksums": cmd_checksums,
    "selftest": cmd_selftest,
}

if __name__ == "__main__":
    if len(sys.argv) < 3 or sys.argv[1] not in COMMANDS:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    COMMANDS[sys.argv[1]](sys.argv[2:])
