#!/usr/bin/env python3
"""tools/helmholtz/ch4-3c2-pairs.py: the pair evidence of Helmholtz CH4.3c2
(docs/helmholtz/tasks/CH4.3c2.md, "Acceptance", pair evidence), from spike S7's per-pair
curves (`spikes/truncation/results/c1_curves.csv`: the one-level chain in extended
precision, per real k a on the half-octave grid and degree p, the worst error over the 316
offsets and three placements, φ and ∇φ, relative ℓ2 per pair).

    tools/helmholtz/ch4-3c2-pairs.py <points> [--rule name=gamma,alpha,beta,c,a,b,d[,cap] …]

<points> is a file of lines `x eps p` (the FMM's measured degree at ε/2 at x =
Re(k) a_l; other lines are skipped) or `point:` lines of `helmholtz_calibrate`'s scans
(`point: level <workload> <kd> <level> <x> <eps> <p at eps/2> <p at eps> <S7>`).

Per FMM point with x ≥ 43 the worst pair's error at the FMM's degree over ε (r), read from
the curves at the point's excess e = (p − x) / x^{1/3}: at each of the two grid k a around x
the degree k a + e (k a)^{1/3}, its log error interpolated linearly in p, then linearly in
log k a. Per rule and grid k a from 45 to 256 and per ε of D1 the worst pair's error at the
rule's degree over ε, beside the largest r of the FMM's points at that ε: the guard of the
plan is that the rule's ratio stays within that largest r (a degree no less resolved per pair
than one the FMM confirmed). Prints Markdown. Python 3 standard library only.
"""

import csv
import math
import sys
from pathlib import Path

CURVES = Path(__file__).resolve().parents[2] / "spikes/truncation/results/c1_curves.csv"
EPSILONS = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
MIN_X = 43.0


def up(v):
    return max(0, math.ceil(v - 1e-9))


class Rule:
    def __init__(self, name, gamma, alpha, beta, c, a=0.0, b=0.0, d=0.0, cap=None):
        self.name = name
        self.k = (gamma, alpha, beta, c, a, b, d)
        self.cap = cap

    def degree(self, x, eps):
        """As `DegreeRule::degree`: the excess at least 0, the branch at most the cap's."""
        gamma, alpha, beta, c, a, b, d = self.k
        big_l = math.log10(1 / eps)
        plateau = up((gamma * big_l + alpha) * big_l + beta)
        excess = max(0.0, c * math.log(1 / eps) ** (2 / 3) + a * big_l + b)
        branch = x + excess * x ** (1 / 3) + d
        if self.cap is not None:
            branch = min(branch, x + self.cap * math.log(1 / eps) ** (2 / 3) * x ** (1 / 3))
        return max(plateau, up(branch))


def load_curves():
    curves = {}
    with open(CURVES) as f:
        for row in csv.DictReader(f):
            if row["kind"] != "real":
                continue
            curves.setdefault(float(row["ka"]), {})[int(row["p"])] = float(row["worst"])
    return curves


def log_error(curve, p):
    """log₁₀ of the worst error at a real degree p, linear between the integers."""
    lo = math.floor(p)
    if lo not in curve or lo + 1 not in curve:
        return None
    t = p - lo
    return (1 - t) * math.log10(curve[lo]) + t * math.log10(curve[lo + 1])


def pair_error(curves, x, p):
    """The worst pair's error at x and degree p, by the excess (module documentation)."""
    grid = sorted(curves)
    if x in curves:
        v = log_error(curves[x], p)
        return None if v is None else 10**v
    below = [g for g in grid if g < x]
    above = [g for g in grid if g > x]
    if not below or not above:
        return None
    g1, g2 = below[-1], above[0]
    e = (p - x) / x ** (1 / 3)
    v1 = log_error(curves[g1], g1 + e * g1 ** (1 / 3))
    v2 = log_error(curves[g2], g2 + e * g2 ** (1 / 3))
    if v1 is None or v2 is None:
        return None
    t = (math.log(x) - math.log(g1)) / (math.log(g2) - math.log(g1))
    return 10 ** ((1 - t) * v1 + t * v2)


def read_points(path):
    points = []
    for line in Path(path).read_text().splitlines():
        f = line.split()
        try:
            if f and f[0] == "point:" and f[1] == "level":
                points.append((float(f[5]), float(f[6]), int(f[7])))
            elif len(f) == 3:
                points.append((float(f[0]), float(f[1]), int(f[2])))
        except (ValueError, IndexError):
            continue
    return points


def main():
    args = sys.argv[1:]
    if not args:
        print(__doc__)
        sys.exit(2)
    points = read_points(args[0])
    rules = [Rule("CH3.3", 0.30, 0.25, 2.60, 0.80)]
    rest = args[1:]
    while rest:
        if rest[0] != "--rule" or len(rest) < 2:
            sys.exit(f"unexpected {rest[0]}")
        name, values = rest[1].split("=")
        rules.append(Rule(name, *[float(v) for v in values.split(",")]))
        rest = rest[2:]
    curves = load_curves()

    print("## The FMM's points: the worst pair at the FMM's degree")
    print()
    print("| x | ε | p (FMM, ε/2) | worst pair / ε |")
    print("| ---: | --- | ---: | ---: |")
    largest = {}
    for x, eps, p in sorted(points):
        if x < MIN_X:
            continue
        err = pair_error(curves, x, p)
        if err is None:
            print(f"| {x:.2f} | {eps:.0e} | {p} | – |")
            continue
        r = err / eps
        largest[eps] = max(largest.get(eps, 0.0), r)
        print(f"| {x:.2f} | {eps:.0e} | {p} | {r:.2f} |")
    print()
    print("Largest per ε: " + ", ".join(f"{e:.0e}: {r:.2f}" for e, r in sorted(largest.items())))
    print()
    print("## The rules at S7's grid: the worst pair at the rule's degree / ε")
    print()
    header = "| k a | rule |" + "".join(f" {e:.0e} |" for e in EPSILONS)
    print(header)
    print("| ---: | --- |" + " ---: |" * len(EPSILONS))
    failures = []
    for ka in sorted(curves):
        if ka < 45:
            continue
        for rule in rules:
            cells = []
            for eps in EPSILONS:
                p = rule.degree(ka, eps)
                v = log_error(curves[ka], p)
                if v is None:
                    cells.append(f"{p}: –")
                    continue
                r = 10**v / eps
                mark = ""
                if eps in largest and r > largest[eps] * (1 + 1e-9):
                    mark = " ✗"
                    failures.append((rule.name, ka, eps, p, r, largest[eps]))
                cells.append(f"{p}: {r:.2f}{mark}")
            print(f"| {ka:.1f} | {rule.name} |" + "".join(f" {c} |" for c in cells))
    print()
    print("Each cell: the rule's degree, then the worst pair's error over ε; ✗ above the "
          "largest of the FMM's points at that ε.")
    print()
    if failures:
        print("Above the FMM's largest ratio:")
        for name, ka, eps, p, r, big in failures:
            print(f"- {name} at k a = {ka:.1f}, ε = {eps:.0e}: p = {p}, {r:.2f} > {big:.2f}")
    else:
        print("Every rule's ratio within the FMM's largest at every grid k a and ε.")


if __name__ == "__main__":
    main()
