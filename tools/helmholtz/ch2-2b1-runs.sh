#!/bin/bash
# tools/helmholtz/ch2-2b1-runs.sh: the machine runs of Helmholtz CH2.2b1, RCR's M2L tables only on
# the levels with an RCR pair (docs/helmholtz/tasks/CH2.2b1.md, "Machine runs"), from the
# repository root after
#
#     cargo build --release -p tetraleaf-validate --example helmholtz_accuracy
#     MAIN=<dir> tools/helmholtz/ch2-2b1-runs.sh <out-dir> <stage…>
#
# MAIN is a directory holding `main`'s `helmholtz_accuracy` (another tree's
# target/release/examples); without it the branch's binary runs in both places. Every case runs in
# four launches, main, branch, branch, main (<name>-main1, -branch1, -branch2, -main2): one setting
# per launch (review 2026-11-24, B12), at one rank with every core (1 x 40 on a Kathleen node),
# ε = 10⁻⁶, gradients, §9.1's samples, the library's defaults named (CH2.2a's `DEFAULTS`),
# `--levels` (the tables "Build" and "RCR tables by level"). Stages:
#   smoke  W1 at 20,000 points, kD = 30
#   w1     W1's sphere at 10⁶ points (`w1:1000000`, CH2.2a's tree), kD = KDS (default 250 500 1000),
#          one case per kD
#   w3     W3 at 10⁶ points (`w3:1000000`, CH4.3b2c's) at kD = 500
# TIMEOUT per launch (default 3600 s), ONLY (a regular expression on the launch names), DRY_RUN
# (print the launches only). Each launch writes <out-dir>/<name>.md and .err and a line in
# <out-dir>/runs.log. Times are read on Kathleen's exclusive nodes only, never asserted; no run on
# the M6 without the user's go-ahead. `tools/helmholtz/ch2-2b1-tables.py` makes the report's
# tables.
set -u
out=$1
shift
mkdir -p "$out"
out=$(cd "$out" && pwd)
target=${CARGO_TARGET_DIR:-target}/release/examples
B=$target/helmholtz_accuracy
A=${MAIN:+$MAIN/helmholtz_accuracy}
A=${A:-$B}
if [ "$(uname)" = Darwin ]; then
    export HWLOC_COMPONENTS=${HWLOC_COMPONENTS--darwin}
    cores=$(sysctl -n hw.perflevel0.physicalcpu 2> /dev/null || echo 6)
else
    # GNU nproc honours OMP_NUM_THREADS: count the cores without it.
    cores=$(env -u OMP_NUM_THREADS nproc)
fi
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 \
    BLIS_NUM_THREADS=1
TIMEOUT=${TIMEOUT:-3600}
KDS=${KDS:-250 500 1000}
DEFAULTS=(--near-v cost --near-wx cost --m2l auto-pw --diagonal auto --diagonal-build interpolated)

# `mpi`: the launcher of one rank of every core.
mpi() {
    if [ "$(uname)" = Darwin ]; then
        echo "mpirun --mca btl_tcp_if_include lo0 --mca oob_tcp_if_include lo0 --bind-to none --host localhost:$cores --map-by slot -n 1"
    else
        echo "mpirun --report-bindings --map-by slot:PE=$cores --bind-to core -n 1"
    fi
}

# `record <name> <command…>`: one launch under the timeout, its output in <name>.md and .err.
record() {
    local n=$1
    shift
    [[ "$n" =~ ${ONLY:-.} ]] || return 0
    if [ -n "${DRY_RUN:-}" ]; then
        echo "$n: $*"
        return 0
    fi
    local start=$SECONDS
    {
        echo "# $(date '+%Y-%m-%d %H:%M:%S') $(hostname) job ${SLURM_JOB_ID:-none}"
        echo "# load: $(uptime)"
        echo "# $*"
    } > "$out/$n.err"
    timeout "$TIMEOUT" "$@" < /dev/null > "$out/$n.md" 2>> "$out/$n.err"
    local status=$?
    echo "$n: exit $status, $((SECONDS - start)) s, $(date +%T)" | tee -a "$out/runs.log"
}

# `case <name> <workload> <arguments…>`: the case in four launches, main, branch, branch, main.
case_abba() {
    local n=$1 w=$2
    shift 2
    local run
    for run in main1:"$A" branch1:"$B" branch2:"$B" main2:"$A"; do
        # shellcheck disable=SC2046
        record "$n-${run%%:*}" $(mpi) "${run#*:}" --threads "$cores" --workload "$w" --eps 1e-6 \
            --levels "$@"
    done
}

for stage in "$@"; do
    case $stage in
    smoke)
        case_abba smoke w1:20000 --kd 30 "${DEFAULTS[@]}"
        ;;
    w1)
        for kd in $KDS; do
            case_abba "w1-$kd" w1:1000000 --kd "$kd" "${DEFAULTS[@]}"
        done
        ;;
    w3)
        case_abba w3-500 w3:1000000 --kd 500 "${DEFAULTS[@]}"
        ;;
    *)
        echo "unknown stage $stage" >&2
        exit 2
        ;;
    esac
done
