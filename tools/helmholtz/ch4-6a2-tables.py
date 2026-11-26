#!/usr/bin/env python3
"""tools/helmholtz/ch4-6a2-tables.py: the tables of Helmholtz CH4.6a2's report (the φ step's time
per level, the product against rustfft; docs/helmholtz/tasks/CH4.6a2.md), from the launches of
`tools/helmholtz/ch4-6a2-runs.sh`.

    python3 tools/helmholtz/ch4-6a2-tables.py runs   <machine>=<dir> …
    python3 tools/helmholtz/ch4-6a2-tables.py levels <machine>=<dir> …
    python3 tools/helmholtz/ch4-6a2-tables.py smooth <machine>=<dir> …
    python3 tools/helmholtz/ch4-6a2-tables.py rule   <machine>=<dir> …
    python3 tools/helmholtz/ch4-6a2-tables.py selftest

A case is a launch name without its last part (`w1-500-1e-6` of `w1-500-1e-6-product-1`); each
case has the launches `-product-1`, `-rustfft-1`, `-rustfft-2`, `-product-2` (ABBA, one setting
per launch) and `-form-product`, `-form-rustfft` (`--diagonal-costs form`).

- `runs`: per case the evaluation (the second of each launch) and the build of each launch, the
  diagonal levels' M2L seconds summed, rustfft's over the product's (the means of the two
  launches of each), the product's spread (the larger launch over the smaller), the errors'
  ratio (φ, ∇φ) and whether the checksums of the launches of one setting agree.
- `levels`: per case and level with the diagonal form, the conversions per level call (S → F,
  F → R), their thread-µs per conversion under each step and the ratio, the pairs' thread-µs per
  pair, and the level's M2L seconds and their ratio (means of two launches each).
- `smooth`: the smooth grid priced per level (the plan's "Machine runs", item 2): rustfft's
  thread-seconds as measured on the product's grid, with S → F and F → R changed by the
  conversions times spike S9's per-box rustfft difference (`ffts-rustfft-smooth.md` less
  `ffts-rustfft-product.md`, at the level's degree) and the pairs scaled by N_φ smooth over N_φ;
  the level's M2L seconds scaled by the thread-seconds' ratio; per case the sum over the levels
  against the evaluation.
- `rule`: per level with V pairs in the time rule's range (k a ≥ 2.8) of each case: what the
  default ran and its M2L seconds (the product launches), the form forced (`-form-product`,
  `-form-rustfft`), and the model's choice with the form's conversions at the measured rustfft
  ÷ product ratio of the case's conversions (`DIAGONAL_COSTS_*`' S → F and F → R units, the
  rest of the model as printed): whether the time rule needs an FFT variant of its costs.

No timing is asserted; the numbers go into tetraleaf-validate/results/helmholtz-ch4-6a2.md.
Python 3.9 standard library only.
"""

import math
import re
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

SETTINGS = ["product-1", "rustfft-1", "rustfft-2", "product-2", "form-product", "form-rustfft"]
MIN_KA = 2.8  # `DIAGONAL_AUTO_MIN_KA`
# The time rule's S → F and F → R units, a (p + 1)² + b (p + 1)³ thread-ns
# (tetraleaf-exec/src/helmholtz/operator_diagonal.rs, `DIAGONAL_COSTS_X86`, `DIAGONAL_COSTS_NEON`).
COSTS = {
    "x86": {"s2f": (383.0, 2.3), "f2r": (313.0, 5.85)},
    "neon": {"s2f": (108.0, 1.58), "f2r": (92.0, 3.26)},
}


def unit(u, p):
    q = p + 1.0
    return u[0] * q * q + u[1] * q * q * q


def mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def ratio(a, b):
    return a / b if a is not None and b else None


def fmt(x, digits=3):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "–"
    return f"{x:.{digits}f}"


def machine_class(path):
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("Machine:"):
            return "neon" if "aarch64" in line else "x86"
    return "x86"


class Launch:
    """One launch's output: its run row, the build, and per level the tables' rows."""

    def __init__(self, path):
        self.path = path
        found = common.tables(path)
        runs = common.accuracy_rows(path)
        self.run = runs[0] if runs else None
        build = common.build_rows(path)
        self.build = common.number(build[0]["build s"]) if build else None
        self.evaluate = common.number(self.run["evaluate s"]) if self.run else None
        self.errors = [common.number(self.run[k]) for k in ("φ", "∇φ")] if self.run else None
        self.checksum = self.run["checksum"] if self.run else None
        self.machine = machine_class(path)

        def by_level(heading):
            t = common.find_table(found, heading=heading, first="Workload")
            return {int(r["level"]): r for r in (t.rows if t else []) if r["level"].isdigit()}

        self.m2l = by_level("M2L per level")
        self.diagonal = by_level("Diagonal M2L")
        self.rule = by_level("Time rule")
        self.threads = None
        if self.run:
            r, t = self.run["ranks × threads"].split("×")
            self.threads = int(r) * int(t)

    def m2l_seconds(self, level):
        r = self.m2l.get(level)
        return common.number(r["M2L s"]) if r else None


def cases(arg):
    """`machine=dir`: the machine and its cases, each {setting: Launch}, in name order."""
    machine, where = arg.split("=", 1)
    found = {}
    for p in sorted(Path(where).glob("*.md")):
        for s in SETTINGS:
            if p.stem.endswith("-" + s):
                found.setdefault(p.stem[: -len(s) - 1], {})[s] = Launch(p)
    return machine, found


def pair(case, a, b):
    return [case[s] for s in (a, b) if s in case]


def diagonal_levels(case):
    first = case.get("product-1") or next(iter(case.values()))
    return sorted(first.diagonal)


def per_conversion(launch, level):
    """[S → F, pairs, F → R] thread-µs per conversion or pair at `level`, and the counts."""
    r = launch.diagonal.get(level)
    if r is None:
        return None
    s2f, f2r = int(r["S → F"]), int(r["F → R"])
    pairs = int(r["diagonal pairs"])
    ts = [common.number(r[k]) for k in ("S → F thread-s", "pairs thread-s", "F → R thread-s")]
    return {
        "s2f": 1e6 * ts[0] / s2f if s2f else None,
        "pair": 1e6 * ts[1] / pairs if pairs else None,
        "f2r": 1e6 * ts[2] / f2r if f2r else None,
        "counts": (s2f, pairs, f2r),
        "ts": ts + [common.number(r["RCR thread-s"])],
        "p": int(r["p"]),
        "classes": r["classes"],
        "phi": r.get("φ", "–"),
    }


def show(header, rows):
    print("| " + " | ".join(header) + " |")
    print("| " + " | ".join("---" for _ in header) + " |")
    for row in rows:
        print("| " + " | ".join(row) + " |")
    print()


def cmd_runs(args):
    header = ["machine", "case", "build s P1, R1, R2, P2", "evaluate s P1, R1, R2, P2",
              "evaluate R ÷ P", "P spread", "diagonal M2L s P, R", "diagonal M2L R ÷ P",
              "errors R ÷ P (φ, ∇φ)", "checksums"]
    rows = []
    for arg in args:
        machine, found = cases(arg)
        for name, case in found.items():
            ab = [case.get(s) for s in SETTINGS[:4]]
            if any(l is None for l in ab):
                continue
            p, r = [ab[0], ab[3]], [ab[1], ab[2]]
            levels = diagonal_levels(case)
            m2l = [mean([sum(l.m2l_seconds(lv) or 0.0 for lv in levels) for l in s]) for s in (p, r)]
            ev = [mean([l.evaluate for l in s]) for s in (p, r)]
            spread = max(l.evaluate for l in p) / min(l.evaluate for l in p)
            errs = [ratio(r[0].errors[i], p[0].errors[i]) for i in (0, 1)]
            same = all(l.checksum == p[0].checksum for l in p) and all(
                l.checksum == r[0].checksum for l in r)
            rows.append([
                machine, name,
                ", ".join(fmt(l.build, 2) for l in ab),
                ", ".join(fmt(l.evaluate, 3) for l in ab),
                fmt(ratio(ev[1], ev[0])), fmt(spread), f"{fmt(m2l[0])}, {fmt(m2l[1])}",
                fmt(ratio(m2l[1], m2l[0])), ", ".join(fmt(e, 4) for e in errs),
                "agree" if same else "**differ**",
            ])
    show(header, rows)


def level_means(case, level):
    """The two settings' per-conversion data at `level`, means over their two launches."""
    out = {}
    for step, launches in (("product", pair(case, "product-1", "product-2")),
                           ("rustfft", pair(case, "rustfft-1", "rustfft-2"))):
        data = [per_conversion(l, level) for l in launches]
        data = [d for d in data if d]
        if not data:
            continue
        out[step] = {k: mean([d[k] for d in data]) for k in ("s2f", "pair", "f2r")}
        out[step]["ts"] = [mean([d["ts"][i] for d in data]) for i in range(4)]
        out[step]["counts"] = data[0]["counts"]
        out[step]["p"] = data[0]["p"]
        out[step]["classes"] = data[0]["classes"]
        out[step]["phi"] = {d["phi"] for d in data}
        out[step]["m2l"] = mean([l.m2l_seconds(level) for l in launches])
    return out


def cmd_levels(args):
    header = ["machine", "case", "level", "p", "classes", "S → F, F → R", "φ ran",
              "S → F µs P, R", "S → F R ÷ P", "F → R µs P, R", "F → R R ÷ P",
              "pair µs P, R", "M2L s P, R", "M2L R ÷ P"]
    rows = []
    for arg in args:
        machine, found = cases(arg)
        for name, case in found.items():
            for level in diagonal_levels(case):
                m = level_means(case, level)
                if "product" not in m or "rustfft" not in m:
                    continue
                P, R = m["product"], m["rustfft"]
                rows.append([
                    machine, name, str(level), str(P["p"]), P["classes"],
                    f"{P['counts'][0]}, {P['counts'][2]}",
                    "/".join(sorted(R["phi"])),
                    f"{fmt(P['s2f'], 1)}, {fmt(R['s2f'], 1)}", fmt(ratio(R["s2f"], P["s2f"])),
                    f"{fmt(P['f2r'], 1)}, {fmt(R['f2r'], 1)}", fmt(ratio(R["f2r"], P["f2r"])),
                    f"{fmt(P['pair'], 2)}, {fmt(R['pair'], 2)}",
                    f"{fmt(P['m2l'], 4)}, {fmt(R['m2l'], 4)}", fmt(ratio(R["m2l"], P["m2l"])),
                ])
    show(header, rows)


def ffts(directory):
    """Spike S9's per-box times of `ffts-<variant>-<grid>.md`: {(variant, grid): {p: (N_φ, S → F
    µs, F → R µs)}}."""
    out = {}
    for path in Path(directory).glob("ffts-*.md"):
        _, variant, grid = path.stem.split("-", 2)
        t = common.find_table(common.tables(path), first="p")
        if t is None:
            continue
        out[(variant, grid)] = {
            int(r["p"]): (int(r["N_φ"]), common.number(r["S → F µs"]), common.number(r["F → R µs"]))
            for r in t.rows
        }
    return out


def cmd_smooth(args):
    header = ["machine", "case", "level", "p", "N_φ, smooth", "rustfft per box S → F + F → R µs, "
              "smooth", "thread-s measured", "thread-s smooth", "smooth ÷ measured",
              "M2L s measured, smooth", "Δ s"]
    rows, totals = [], []
    for arg in args:
        machine, found = cases(arg)
        where = arg.split("=", 1)[1]
        times = ffts(where)
        prod, smooth = times.get(("rustfft", "product"), {}), times.get(("rustfft", "smooth"), {})
        for name, case in found.items():
            delta_case = 0.0
            for level in diagonal_levels(case):
                m = level_means(case, level).get("rustfft")
                if not m:
                    continue
                p = m["p"]
                if p not in prod or p not in smooth:
                    rows.append([machine, name, str(level), str(p), "no S9 time at this p"] +
                                ["–"] * 6)
                    continue
                n_p, ps2f, pf2r = prod[p]
                n_s, ss2f, sf2r = smooth[p]
                s2f_n, _, f2r_n = m["counts"]
                ts = m["ts"]
                priced = [
                    ts[0] + s2f_n * (ss2f - ps2f) * 1e-6,
                    ts[1] * n_s / n_p,
                    ts[2] + f2r_n * (sf2r - pf2r) * 1e-6,
                    ts[3],
                ]
                r = sum(priced) / sum(ts)
                wall = m["m2l"]
                delta = wall * (r - 1.0)
                delta_case += delta
                rows.append([
                    machine, name, str(level), str(p), f"{n_p}, {n_s}",
                    f"{fmt(ps2f + pf2r, 1)}, {fmt(ss2f + sf2r, 1)}", fmt(sum(ts), 4),
                    fmt(sum(priced), 4), fmt(r), f"{fmt(wall, 4)}, {fmt(wall * r, 4)}",
                    fmt(delta, 4),
                ])
            ev = mean([case[s].evaluate for s in ("rustfft-1", "rustfft-2") if s in case])
            totals.append([machine, name, fmt(ev), fmt(delta_case, 4),
                           fmt(ratio(ev + delta_case, ev) if ev else None)])
    show(header, rows)
    show(["machine", "case", "evaluate s (rustfft)", "Δ s of the smooth grid",
          "evaluate smooth ÷ product's grid"], totals)


def model_choice(printed, bound, diag, rcr, waves):
    """The time rule's choice on a level with the model's diagonal level at `diag` (the example's
    "model's choice (no budget)", restated): on a level with Auto's waves, or one the rule took
    from them, the waves unless the form beats them, then the form against RCR; elsewhere the form
    against the level's own M2L ("other"). Whether the level has waves is read from the choice
    `printed` with the model's own diagonal level: "waves" and "rcr" only on such a level, "other"
    only elsewhere; a lower `diag` keeps a "diagonal" choice in both cases."""
    if printed == "diagonal":
        return "diagonal"
    waves_level = printed in ("waves", "rcr")
    if waves_level and not (bound > 0 and diag < waves):
        return "waves"
    if bound > 0 and diag < rcr:
        return "diagonal"
    return "rcr" if waves_level else "other"


def cmd_rule(args):
    header = ["machine", "case", "level", "p", "k a", "bound", "default ran", "default M2L s",
              "form P M2L s", "form R M2L s", "form R ÷ default", "model diag, RCR, waves s",
              "conversions' share of diag", "model diag with R", "choice P → R"]
    rows = []
    for arg in args:
        machine, found = cases(arg)
        for name, case in found.items():
            first = case.get("product-1")
            if first is None:
                continue
            costs = COSTS[first.machine]
            # The case's measured rustfft ÷ product ratio of the conversions' thread time.
            num = den = 0.0
            for level in diagonal_levels(case):
                m = level_means(case, level)
                if "product" in m and "rustfft" in m:
                    num += m["rustfft"]["ts"][0] + m["rustfft"]["ts"][2]
                    den += m["product"]["ts"][0] + m["product"]["ts"][2]
            rho = num / den if den else None
            for level, r in sorted(first.rule.items()):
                ka = common.number(r["k a"])
                if ka is None or ka < MIN_KA:
                    continue
                p = int(r["p"])
                bound = int(common.number(r["bound classes"], 0))
                diag, rcr, waves = (common.number(r[k]) for k in
                                    ("model diagonal s", "model RCR s", "model waves s"))
                boxes = common.number(first.m2l[level]["boxes"]) if level in first.m2l else None
                conv = None
                if boxes is not None and first.threads:
                    conv = boxes * (unit(costs["s2f"], p) + unit(costs["f2r"], p)) * 1e-9 / \
                        first.threads
                diag_r = diag - (1.0 - rho) * conv if None not in (diag, conv, rho) else None
                before = r["model's choice (no budget)"]
                after = model_choice(before, bound, diag_r, rcr, waves) if diag_r is not None \
                    else before
                default = mean([case[s].m2l_seconds(level) for s in ("product-1", "product-2")
                                if s in case])
                fp = case["form-product"].m2l_seconds(level) if "form-product" in case else None
                fr = case["form-rustfft"].m2l_seconds(level) if "form-rustfft" in case else None
                ran_fr = case["form-rustfft"].rule.get(level, {}).get("ran", "") \
                    if "form-rustfft" in case else ""
                rows.append([
                    machine, name, str(level), str(p), fmt(ka, 2), str(bound), r["ran"],
                    fmt(default, 4), fmt(fp, 4),
                    fmt(fr, 4) + ("" if ran_fr.startswith("diagonal") else f" ({ran_fr})"),
                    fmt(ratio(fr, default)),
                    f"{fmt(diag, 4)}, {fmt(rcr, 4)}, {fmt(waves, 4)}",
                    fmt(ratio(conv, diag)), fmt(diag_r, 4),
                    before if before == after else f"**{before} → {after}**",
                ])
            print(f"{machine} {name}: rustfft ÷ product of the conversions' thread time "
                  f"{fmt(rho)}")
    print()
    show(header, rows)


def selftest():
    assert abs(unit((1.0, 1.0), 1) - 12.0) < 1e-12
    assert model_choice("other", 34, 1.0, 2.0, 9.0) == "diagonal"
    assert model_choice("other", 34, 3.0, 2.0, 9.0) == "other"
    assert model_choice("waves", 34, 1.0, 2.0, 1.5) == "diagonal"
    assert model_choice("waves", 34, 1.6, 2.0, 1.5) == "waves"
    assert model_choice("waves", 0, 1.0, 2.0, 1.5) == "waves"
    assert model_choice("rcr", 34, 1.4, 1.2, 1.5) == "rcr"
    assert model_choice("diagonal", 34, 9.0, 2.0, 5.0) == "diagonal"
    assert re.match(r"\d", "2")
    assert statistics.median([1, 3]) == 2
    print("selftest: ok")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    mode, args = sys.argv[1], sys.argv[2:]
    commands = {"runs": cmd_runs, "levels": cmd_levels, "smooth": cmd_smooth, "rule": cmd_rule}
    if mode == "selftest":
        selftest()
    elif mode in commands:
        commands[mode](args)
    else:
        print(__doc__)
        sys.exit(2)


if __name__ == "__main__":
    main()
