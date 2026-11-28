#!/bin/bash
# tools/helmholtz/ch2-2a-runs.sh: the machine runs of Helmholtz CH2.2a, the RCR tables per level in
# a store (docs/helmholtz/tasks/CH2.2a.md, "Machine runs"), from the repository root after
#
#     cargo build --release -p tetraleaf-validate --example helmholtz_accuracy
#     tools/helmholtz/ch2-2a-runs.sh <out-dir> <stage…>
#
# Every launch runs `helmholtz_accuracy` at one rank with every core (1 x 40 on a Kathleen node) on
# W1's sphere at a fixed N (`w1:<n>`: the points, and so the FMM's domain, the same at every kD, so
# that k̂₀ doubles exactly with kD), a sweep of kD in octaves in one launch, ε = 10⁻⁶, gradients,
# §9.1's samples, `--levels` (the table "RCR tables by level"). Each sweep runs twice, in two
# launches: with `--table-store` (one store per launch, shared by the sweep's cases in order: warm
# from the second kD on) and without (every case cold); the tables script compares the two
# launches per kD (parts, thread-seconds, MB, checksums). Times are reported, not decided on, and
# never compared within a launch (review 2026-11-24, B12). Stages:
#   smoke     W1 at 20,000 points, kD = 15, 30, Rotation
#   real      W1 at 10⁶ points, kD = KDS (default 62.5 125 250 500 1000) under the library's
#             defaults named (`--m2l auto-pw --diagonal auto --diagonal-build interpolated`, the V and
#             W/X rules of the cost model)
#   rotation  the same sweep under Rotation (RCR's M2L on every level: every level's M2L tables
#             in the store)
#   yukawa    W1 at 10⁶ points, Yukawa at κD = KAPPAS (default 15.625 31.25 62.5 125 250 500 1000;
#             `--kd 0 --kim …`) under the defaults named (Rotation on Yukawa)
#   repeat    W1 at 10⁶ points, kD = 1000 twice under the defaults named: the second run's tables
#             all from the store (the same k again, the bound of what a store saves)
# TIMEOUT per launch (default 3600 s), STORE_MB (default 65536: no eviction on a Kathleen node),
# ONLY (a regular expression on the launch names), DRY_RUN (print the launches only). Each launch
# writes <out-dir>/<name>.md and .err and a line in <out-dir>/runs.log. Times are read on Kathleen's
# exclusive nodes only, never asserted; no run on the M6 without the user's go-ahead.
# `tools/helmholtz/ch2-2a-tables.py` makes the report's tables.
set -u
out=$1
shift
mkdir -p "$out"
out=$(cd "$out" && pwd)
target=${CARGO_TARGET_DIR:-target}/release/examples
A=$target/helmholtz_accuracy
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
STORE_MB=${STORE_MB:-65536}
KDS=${KDS:-62.5 125 250 500 1000}
KAPPAS=${KAPPAS:-15.625 31.25 62.5 125 250 500 1000}
DEFAULTS=(--near-v cost --near-wx cost --m2l auto-pw --diagonal auto --diagonal-build interpolated)
N=1000000

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

# `sweep <name> <workload> <arguments…>`: the sweep with a store (<name>-store), then without
# (<name>-cold), each a launch at every core.
sweep() {
    local n=$1 w=$2
    shift 2
    # shellcheck disable=SC2046
    record "$n-store" $(mpi) "$A" --threads "$cores" --workload "$w" --eps 1e-6 --levels \
        --table-store "$STORE_MB" "$@"
    # shellcheck disable=SC2046
    record "$n-cold" $(mpi) "$A" --threads "$cores" --workload "$w" --eps 1e-6 --levels "$@"
}

list() {
    echo "$*" | tr ' ' ','
}

for stage in "$@"; do
    case $stage in
    smoke)
        sweep smoke w1:20000 --kd 15,30 --m2l rotation
        ;;
    real)
        sweep real "w1:$N" --kd "$(list $KDS)" "${DEFAULTS[@]}"
        ;;
    rotation)
        sweep rotation "w1:$N" --kd "$(list $KDS)" --m2l rotation
        ;;
    repeat)
        sweep repeat "w1:$N" --kd 1000,1000 "${DEFAULTS[@]}"
        ;;
    yukawa)
        sweep yukawa "w1:$N" --kd 0 --kim "$(list $KAPPAS)" "${DEFAULTS[@]}"
        ;;
    *)
        echo "unknown stage $stage" >&2
        exit 2
        ;;
    esac
done
