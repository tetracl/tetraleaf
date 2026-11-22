#!/bin/bash
# tools/helmholtz/ch4-4a-runs.sh: the machine runs of Helmholtz CH4.4a, the kD sweep beside FMM3D
# (docs/helmholtz/tasks/CH4.4a.md, "Machine runs"), from the repository root:
#
#     tools/helmholtz/ch4-4a-runs.sh <out-dir> <stage…>
#
# On Kathleen inside a job (tools/kathleen/README.md), one node, for example
#
#     tools/kathleen/submit.sh --time=06:00:00 -- run.sbatch \
#         'tools/helmholtz/ch4-4a-runs.sh /scratch/scratch/ucahtbe/fmm/logs/ch4-4a low'
#
# Each stage is one call of `tools/fmm3d/run.sh helmholtz` (README.md there) on its cases, with
# our side under the library's defaults (`H_SETTINGS=defaults`) at ε = 10⁻³, 10⁻⁶, 10⁻⁸, FMM3D at
# the same eps (the variants of VARIANTS, by default `default fastker`; `default` alone on the
# density cases) and `default` at 10⁻¹² for the full-target crosscheck, no one-thread cases, and
# under GNU time FMM3D's peak memory per eps. `tools/helmholtz/ch4-4a-tables.py` makes the
# report's tables from <out-dir>'s helmholtz-<case>.md. Stages:
#   smoke    W1 at kD = 10, and W2 at kD = 10 with the constant density (a check of the script)
#   low      W1 and W2 at kD = 10, 30, 100, 250; W3 resolved at kD = 10, 30
#   const    the constant density: W1 and W2 at kD = 10, 30, 100; W3 resolved at kD = 10, 30
#   mid      W1 and W2 at kD = 500
#   pw500    W1 at kD = 500 with the plane wave's trace
#   w1-1000  W1 at kD = 1000
#   pw1000   W1 at kD = 1000 with the plane wave's trace
#   w3       W3 resolved at kD = 100; W3 at N = 10⁶ for kD = 200 and 500
# DRY_RUN (print the cases only) and TIMEOUT (per launch, run.sh's, default 7200 s) pass through.
set -eu
out=$1
shift
mkdir -p "$out"
export OUT=$out H_SETTINGS=defaults H_EPS=${H_EPS:-1e-3,1e-6,1e-8} H_T1_CASES=
both=${VARIANTS:-default fastker}

# `cases <variants> <case…>`: one call of run.sh's `helmholtz` step.
cases() {
    local variants=$1
    shift
    echo "=== $(date '+%Y-%m-%d %H:%M:%S') cases: $* (FMM3D: $variants)"
    VARIANTS=$variants H_CASES="$*" tools/fmm3d/run.sh helmholtz
}

for stage in "$@"; do
    case "$stage" in
    smoke) cases "$both" w1@10 && cases default w2@10/constant ;;
    low) cases "$both" w1@10 w1@30 w1@100 w1@250 w2@10 w2@30 w2@100 w2@250 w3@10 w3@30 ;;
    const)
        cases default w1@10/constant w1@30/constant w1@100/constant w2@10/constant \
            w2@30/constant w2@100/constant w3@10/constant w3@30/constant
        ;;
    mid) cases "$both" w1@500 w2@500 ;;
    pw500) cases default w1@500/planewave ;;
    w1-1000) cases "$both" w1@1000 ;;
    pw1000) cases default w1@1000/planewave ;;
    w3) cases "$both" w3@100 w3:1000000@200 w3:1000000@500 ;;
    *)
        echo "ch4-4a-runs.sh: unknown stage $stage (smoke, low, const, mid, pw500, w1-1000, pw1000, w3)" >&2
        exit 2
        ;;
    esac
done
