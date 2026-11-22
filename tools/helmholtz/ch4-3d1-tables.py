#!/usr/bin/env python3
"""tools/helmholtz/ch4-3d1-tables.py: the tables of Helmholtz CH4.3d1's report (the transfer
function's length, per pair; docs/helmholtz/tasks/CH4.3d1.md), from the CSV of S6's spike's
`lengths` example (`spikes/diagonal/examples/lengths.rs`).

    python3 tools/helmholtz/ch4-3d1-tables.py <lengths.csv> [section …]

Sections (all by default): `lmin` (the shortest admissible length per (kind, k a, ε) for the
nearest class and over the classes, against 2p, ⌈3p/2⌉ and the excess-bandwidth length),
`named` (the named lengths ⌈3p/2⌉ and 2ka + 1.8 (log₁₀ 1/ε)^{2/3} (2ka)^{1/3}: admissible on
how many classes), `model` (T's truncation against the model's sum: the constant and the
spread), `round` (the rounding against E at length L), `gain` (per level of W1 and W3: the
grid a length rule allows and the classes it takes from RCR), `rules` (three
length rules on the classes the bound takes: the measured lengths, the model's and Song and
Chew's, their grids and failures), `trunc` (T's truncation per length on the nearest class).

Every number is the worst over the three placements of max(potential, gradient), S7's relative
ℓ2 norm over the targets. A length L is admissible for a class at (k a, ε) when the f64
chain's total error at L is at most max(1.1 × RCR's, RCR's + 0.1 ε) (the plan's "within RCR's
error"). A unit (kind, k a, p) serves every ε whose rule gives p (the `eps` column). No timing.
"""

import collections
import csv
import math
import sys

U = 2.0**-53
U_DD = 2.0**-104
C_T, C_0, THETA = 0.15, 200.0, 0.5
A = math.sqrt(3.0) / 2.0  # the circumradius of a unit box
EPSILONS = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
# The class sizes (offsets per class, of 316): S6's `classes()`, by representative.
SIZES = {}


def worst(r, a, b):
    return max(float(r[a]), float(r[b]))


def load(path):
    """units[(kind, ka, p)] = {"eps": [...], "classes": {rep: {"d2", "trunc", "L": {L: {...}}}}}"""
    units = {}
    for r in csv.DictReader(open(path)):
        key = (r["kind"], float(r["ka"]), int(r["p"]))
        u = units.setdefault(key, {"eps": [float(e) for e in r["eps"].split(";")], "classes": {}})
        c = u["classes"].setdefault(r["rep"], {"d2": int(r["d2"]), "trunc": 0.0, "L": {}})
        c["trunc"] = max(c["trunc"], worst(r, "trunc_pot", "trunc_grad"))
        L = int(r["l"])
        e = c["L"].setdefault(
            L,
            {
                "dd": 0.0,
                "dev": 0.0,
                "round": 0.0,
                "total": 0.0,
                "term": float(r["log10_term"]),
                "tail": float(r["log10_tail"]),
            },
        )
        for k in ("dd", "dev", "round", "total"):
            v = worst(r, k + "_pot", k + "_grad")
            if math.isnan(v):
                e[k] = float("nan")
            elif not math.isnan(e[k]):
                e[k] = max(e[k], v)
    return units


def damping(kind, ka, d2):
    """e^{Im(k)(|t| + 2a)} in box widths, k = κ (1 + iη), κ = ka / a."""
    if kind == "real":
        return 1.0
    eta = float(kind.replace("complex", ""))
    return math.exp(eta * ka / A * (math.sqrt(d2) + 2 * A))


def bound(kind, ka, d2, term):
    """E at a length, from log₁₀ of its term (proposal §3.6)."""
    if term > 300:
        return float("inf")
    return C_T * U * 10.0**term * damping(kind, ka, d2) + C_0 * U


def admissible(c, L, eps):
    e = c["L"][L]
    return e["total"] <= max(1.1 * c["trunc"], c["trunc"] + 0.1 * eps)


def lmin(c, eps):
    """The shortest admissible scanned length (None if none is): not required to start an
    admissible run up to 2p, since where rounding dominates 2p itself may fail."""
    for L in sorted(c["L"]):
        if admissible(c, L, eps):
            return L
    return None


def excess_bandwidth(ka, eps):
    x = 2 * ka
    return math.ceil(x + 1.8 * math.log10(1 / eps) ** (2 / 3) * x ** (1 / 3) - 1e-9)


def pad4(n):
    return -(-n // 4) * 4


def samples(p, L):
    """The grid exact for degree 2p + L: Gauss–Legendre rings and azimuths padded to 4."""
    return -(-(2 * p + L + 1) // 2) * pad4(2 * p + L + 1)


def rows_by_eps(units):
    """(kind, ka, eps) → (p, unit)."""
    out = {}
    for (kind, ka, p), u in units.items():
        for eps in u["eps"]:
            out[(kind, ka, eps)] = (p, u)
    return out


NEAREST = "0:0:2"


def section_lmin(units):
    print("## The shortest admissible length\n")
    print(
        "| kind | k a | ε | p | 2p | ⌈3p/2⌉ | 2ka + 1.8… | nearest L_min | max L_min | "
        "max L_min − 2p | classes at 2p only |"
    )
    print("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for (kind, ka, eps), (p, u) in sorted(rows_by_eps(units).items()):
        cs = u["classes"]
        ls = {rep: lmin(c, eps) for rep, c in cs.items()}
        known = [v for v in ls.values() if v is not None]
        near = ls.get(NEAREST)
        top = max(known) if known else None
        at_full = sum(1 for v in ls.values() if v == 2 * p)
        none = sum(1 for v in ls.values() if v is None)
        print(
            f"| {kind} | {ka:.1f} | {eps:.0e} | {p} | {2 * p} | {(3 * p + 1) // 2} | "
            f"{excess_bandwidth(ka, eps)} | {near} | {top} | "
            f"{'' if top is None else top - 2 * p} | {at_full}{f' (+{none} none)' if none else ''} |"
        )
    print()


def section_named(units):
    print("## The named lengths: classes admissible (of the classes run)\n")
    print("| kind | k a | ε | p | ⌈3p/2⌉ | admissible | 2ka + 1.8… | admissible |")
    print("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for (kind, ka, eps), (p, u) in sorted(rows_by_eps(units).items()):
        cs = u["classes"]
        out = []
        for L in ((3 * p + 1) // 2, excess_bandwidth(ka, eps)):
            if L in next(iter(cs.values()))["L"]:
                n = sum(1 for c in cs.values() if admissible(c, L, eps))
                out += [str(L), f"{n}/{len(cs)}"]
            elif L >= 2 * p:
                out += [f"{L} (≥ 2p)", "–"]
            else:
                out += [f"{L} (< scan)", "–"]
        print(f"| {kind} | {ka:.1f} | {eps:.0e} | {p} | " + " | ".join(out) + " |")
    print()


def section_model(units):
    print("## T's truncation against the model\n")
    ratios = []
    per_ka = collections.defaultdict(list)
    for (kind, ka, p), u in units.items():
        for rep, c in u["classes"].items():
            for L, e in c["L"].items():
                if L >= 2 * p or math.isnan(e["dd"]) or not math.isfinite(e["tail"]):
                    continue
                if e["dd"] < 1e-20:
                    continue
                # The double-double chain must resolve T's truncation: its own rounding,
                # C_T u_dd 10^term (u_dd = 2⁻¹⁰⁴), at most 1 % of it.
                if e["dd"] < 100 * U_DD * 10.0 ** min(e["term"], 300) * damping(
                    kind, ka, c["d2"]
                ):
                    continue
                r = math.log10(e["dd"]) - e["tail"]
                ratios.append((r, kind, ka, p, rep, L))
                per_ka[(kind, ka)].append(r)
    if not ratios:
        print("no rows\n")
        return
    ratios.sort()
    top = ratios[-1]
    med = ratios[len(ratios) // 2][0]
    c_tau = 10 ** top[0] / 2
    print(
        f"{len(ratios)} (unit, class, length) rows with the double-double chain (p ≤ 100), T's "
        f"truncation above 10⁻²⁰ and at least 100× the double-double chain's own rounding. log₁₀(measured/τ-sum): max {top[0]:.2f} at {top[1:]}, median "
        f"{med:.2f}, min {ratios[0][0]:.2f}. C_τ = {c_tau:.3g} puts the largest at 2× (the median "
        f"over-prediction then {10 ** (math.log10(c_tau) - med):.3g}×).\n"
    )
    print("| kind | k a | rows | max | median | min |")
    print("| --- | --- | --- | --- | --- | --- |")
    for key in sorted(per_ka):
        v = sorted(per_ka[key])
        print(
            f"| {key[0]} | {key[1]:.1f} | {len(v)} | {v[-1]:.2f} | {v[len(v) // 2]:.2f} | "
            f"{v[0]:.2f} |"
        )
    print()


def section_round(units):
    print("## The rounding against E at the length\n")
    worst_r = []
    for (kind, ka, p), u in units.items():
        for rep, c in u["classes"].items():
            for L, e in c["L"].items():
                if math.isnan(e["round"]):
                    continue
                E = bound(kind, ka, c["d2"], e["term"])
                worst_r.append((e["round"] / E, kind, ka, p, rep, L, L == 2 * p))
    worst_r.sort()
    full = [w for w in worst_r if w[6]]
    short = [w for w in worst_r if not w[6]]
    for name, v in (("L = 2p", full), ("L < 2p", short)):
        if v:
            print(
                f"- {name}: {len(v)} rows, rounding/E max {v[-1][0]:.3f} at {v[-1][1:6]}, "
                f"median {v[len(v) // 2][0]:.3g}; above 2: {sum(1 for w in v if w[0] > 2)}"
            )
    print()


# The diagonal levels of W1: κa = (√3/2) kD / 2^l for the levels on the sweep's k a.
W1_LEVELS = {
    250: [54.127, 27.063, 13.532],
    500: [108.253, 54.127, 27.063, 13.532],
    1000: [216.506, 108.253, 54.127, 27.063, 13.532],
}


def section_gain(units):
    print("## The gain of a length rule on W1's levels (real k)\n")
    print(
        "Per level: the classes the rule takes now (E at 2p ≤ θε), the longest L_min among "
        "them (sets the grid), the grid's samples at that length over 2p's, the classes "
        "on RCR now that some admissible length puts within θε (their offsets of 316).\n"
    )
    print("| kD | k a | ε | p | diagonal now | longest L_min | samples ratio | from RCR |")
    print("| --- | --- | --- | --- | --- | --- | --- | --- |")
    by = rows_by_eps(units)
    for kd, kas in W1_LEVELS.items():
        for ka in kas:
            for eps in EPSILONS:
                hit = [k for k in by if k[0] == "real" and abs(k[1] - ka) < 2e-3 and k[2] == eps]
                if not hit:
                    continue
                p, u = by[hit[0]]
                diag, longest, moved, moved_sz = 0, None, 0, 0
                for rep, c in u["classes"].items():
                    e2p = bound("real", ka, c["d2"], c["L"][2 * p]["term"])
                    lm = lmin(c, eps)
                    if e2p <= THETA * eps:
                        diag += 1
                        if lm is not None:
                            longest = lm if longest is None else max(longest, lm)
                    else:
                        ok = [
                            L
                            for L in sorted(c["L"])
                            if admissible(c, L, eps)
                            and bound("real", ka, c["d2"], c["L"][L]["term"]) <= THETA * eps
                        ]
                        if ok:
                            moved += 1
                            moved_sz += SIZES.get(rep, 0)
                ratio = "" if longest is None else f"{samples(p, longest) / samples(p, 2 * p):.3f}"
                print(
                    f"| {kd} | {ka:.1f} | {eps:.0e} | {p} | {diag}/{len(u['classes'])} | "
                    f"{longest} | {ratio} | {moved}{f' ({moved_sz})' if moved else ''} |"
                )
    print()


def c_tau(units):
    """The model's constant: the largest measured/τ-sum over the resolved rows, halved (the
    largest at 2×, as S6's C_T)."""
    top = -math.inf
    for (kind, ka, p), u in units.items():
        for c in u["classes"].values():
            for L, e in c["L"].items():
                if L >= 2 * p or math.isnan(e["dd"]) or not math.isfinite(e["tail"]):
                    continue
                if e["dd"] < 1e-20 or e["dd"] < 100 * U_DD * 10.0 ** min(
                    e["term"], 300
                ) * damping(kind, ka, c["d2"]):
                    continue
                top = max(top, math.log10(e["dd"]) - e["tail"])
    return 10**top / 2


def section_rules(units):
    print("## Three length rules on the classes the bound takes (E at 2p ≤ θε)\n")
    ct = c_tau(units)
    print(
        "Per (kind, k a, ε): the classes the bound takes at 2p; the grid's samples at the "
        "longest length a rule gives them over 2p's: **measured**, the shortest admissible "
        f"lengths (an oracle); **model**, the shortest scanned L with 2 C_τ τ_L ≤ 0.1 ε "
        f"(C_τ = {ct:.3g}, the model's bound on T's truncation at a tenth of ε); **1.8**, "
        "Song and Chew's length (capped at 2p), and the classes on which it is not "
        "admissible. Then the classes on which the model's length is not admissible.\n"
    )
    print("| kind | k a | ε | p | classes | measured | model | 1.8 | 1.8 fails | model fails |")
    print("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    totals = collections.Counter()
    for (kind, ka, eps), (p, u) in sorted(rows_by_eps(units).items()):
        if ka < 11:
            continue
        meas, model, sc, sc_fail, model_fail, n = None, None, None, 0, 0, 0
        lsc = min(excess_bandwidth(ka, eps), 2 * p)
        for c in u["classes"].values():
            if bound(kind, ka, c["d2"], c["L"][2 * p]["term"]) > THETA * eps:
                continue
            n += 1
            ls = sorted(c["L"])
            lm = lmin(c, eps)
            lmod = next((L for L in ls if 2 * ct * 10 ** c["L"][L]["tail"] <= 0.1 * eps), 2 * p)
            meas = max(meas or 0, lm if lm is not None else 2 * p)
            model = max(model or 0, lmod)
            if not admissible(c, lmod, eps):
                model_fail += 1
            l18 = next(L for L in ls if L >= lsc)
            if not admissible(c, l18, eps):
                sc_fail += 1
        if n == 0:
            continue
        totals["levels"] += 1
        totals["model_fail"] += model_fail
        totals["classes"] += n
        r = [samples(p, x) / samples(p, 2 * p) for x in (meas, model, lsc)]
        print(
            f"| {kind} | {ka:.1f} | {eps:.0e} | {p} | {n} | {r[0]:.3f} | {r[1]:.3f} | "
            f"{r[2]:.3f} | {sc_fail} | {model_fail} |"
        )
    print(f"\n{totals['levels']} rows, {totals['classes']} classes, the model failing on "
          f"{totals['model_fail']}.\n")


def section_trunc(units):
    print("## T's truncation on the nearest class (double-double chain against RCR)\n")
    for (kind, ka, p), u in sorted(units.items()):
        if kind != "real" or p > 100 or NEAREST not in u["classes"]:
            continue
        c = u["classes"][NEAREST]
        cells = " ".join(
            f"{L}:{e['dd']:.1e}/{e['total']:.1e}" for L, e in sorted(c["L"].items())[-6:]
        )
        print(f"- k a {ka:.1f}, p {p}, RCR {c['trunc']:.1e}: {cells}")
    print()


def main():
    path = sys.argv[1]
    want = sys.argv[2:] or ["lmin", "named", "model", "round", "gain", "rules", "trunc"]
    units = load(path)
    for name in want:
        globals()["section_" + name](units)


def sizes():
    """The class sizes of S6's 34 classes: offsets o with max |oᵢ| ≤ 3 and not all |oᵢ| ≤ 1,
    by representative (min(|a|, |b|), max(|a|, |b|), |c|)."""
    for a in range(-3, 4):
        for b in range(-3, 4):
            for c in range(-3, 4):
                if max(abs(a), abs(b), abs(c)) <= 1:
                    continue
                rep = f"{min(abs(a), abs(b))}:{max(abs(a), abs(b))}:{abs(c)}"
                SIZES[rep] = SIZES.get(rep, 0) + 1
    assert sum(SIZES.values()) == 316 and len(SIZES) == 34, (sum(SIZES.values()), len(SIZES))


if __name__ == "__main__":
    sizes()
    main()
