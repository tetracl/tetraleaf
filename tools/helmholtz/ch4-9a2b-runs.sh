#!/bin/bash
# tools/helmholtz/ch4-9a2b-runs.sh: the machine runs of Helmholtz CH4.9a2b, P2M by the batched
# kernel (docs/helmholtz/tasks/CH4.9a2b.md, "Machine runs"), from the repository root after
#
#     cargo build --release -p tetraleaf-validate --example helmholtz_accuracy --example helmholtz_leaf
#     tools/helmholtz/ch4-9a2b-runs.sh <out-dir> <stage…>
#
# Every `helmholtz_accuracy` launch runs at one rank with every core, `--levels`, §9.1's complex
# charges and samples, gradients, two evaluations (the rows the second's), under the library's
# defaults named (`--m2l auto-pw --diagonal auto --diagonal-build interpolated`, the V and W/X
# rules of the cost model), and one P2M per launch (review 2026-11-24, B12: on Kathleen a second
# setting in one launch runs 1.09–1.20× slower): `--leaves kernel` (the batched kernel, the
# library's default) or `--leaves reference` (the per-point reference calls; P2P stays the
# vectorised kernel). `helmholtz_leaf` runs on one thread pinned to one core where `taskset`
# exists (CORE, default 2), one precision and degree per launch, 16 lanes.
# `tools/helmholtz/ch4-9a2b-tables.py` makes the report's tables. Stages:
#   smoke   `helmholtz_leaf` at p = 8 and 17, k̂ = 2, 64 points; W1 at kD = 30, 10⁻⁶, the kernel
#           then the reference (a check of the script and the tables)
#   levels  W1 at kD = 250, 500, 1000 at 10⁻⁶, W1 at kD = 500 at 10⁻³ and 10⁻⁸, W3 at N = 10⁶,
#           kD = 500 at 10⁻⁶ (CH4.9's cases): per case kernel, reference, reference, kernel (ABBA)
#   worst   W1's and W2's constant density at kD = 500, 10⁻⁶ (the worst density the record
#           names, CH4.3c3 and CH4.8): the kernel, then the reference (errors)
#   parts   `helmholtz_leaf` at PS (default 8 12 17 26 32 45) × KHAT (default 0.5,2,4,8,16),
#           f64 then f32, 256 points as one leaf
# TIMEOUT per launch (default 3600 s), ONLY (a regular expression on the launch names),
# DRY_RUN (print the launches only). Each launch writes <out-dir>/<name>.md and .err and a line
# in <out-dir>/runs.log.
set -u
out=$1
shift
mkdir -p "$out"
out=$(cd "$out" && pwd)
target=${CARGO_TARGET_DIR:-target}/release/examples
A=$target/helmholtz_accuracy
L=$target/helmholtz_leaf
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
PS=${PS:-8 12 17 26 32 45}
KHAT=${KHAT:-0.5,2,4,8,16}
COMMON=(--near-v cost --near-wx cost --levels --m2l auto-pw --evaluations 2
    --diagonal auto --diagonal-build interpolated)
pin=()
command -v taskset > /dev/null && pin=(taskset -c "${CORE:-2}")

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

# `launch <name> <arguments…>`: one launch of the harness at one rank with every core.
launch() {
    local n=$1
    shift
    # shellcheck disable=SC2046
    record "$n" $(mpi) "$A" "${COMMON[@]}" --threads "$cores" "$@"
}

# `leaf <name> <arguments…>`: one launch of `helmholtz_leaf` on one core.
leaf() {
    local n=$1
    shift
    record "$n" ${pin[@]+"${pin[@]}"} "$L" --lanes 16 "$@"
}

# `case_runs <name> <arguments…>`: one case's four launches, ABBA (the stage `levels`).
case_runs() {
    local n=$1
    shift
    launch "$n-kernel-1" --leaves kernel "$@"
    launch "$n-reference-1" --leaves reference "$@"
    launch "$n-reference-2" --leaves reference "$@"
    launch "$n-kernel-2" --leaves kernel "$@"
}

for stage in "$@"; do
    case $stage in
    smoke)
        leaf smoke-leaf --ps 8,17 --khat 2 --points 64
        launch smoke-kernel --leaves kernel --workload w1 --kd 30 --eps 1e-6
        launch smoke-reference --leaves reference --workload w1 --kd 30 --eps 1e-6
        ;;
    levels)
        for kd in 250 500 1000; do
            case_runs "w1-$kd-1e-6" --workload w1 --kd "$kd" --eps 1e-6
        done
        for eps in 1e-3 1e-8; do
            case_runs "w1-500-$eps" --workload w1 --kd 500 --eps "$eps"
        done
        case_runs "w3-500-1e-6" --workload w3:1000000 --kd 500 --eps 1e-6
        ;;
    worst)
        for w in w1 w2; do
            for leaves in kernel reference; do
                launch "worst-$w-500-1e-6-$leaves" --leaves "$leaves" --workload "$w" --kd 500 \
                    --eps 1e-6 --density constant
            done
        done
        ;;
    parts)
        for precision in f64 f32; do
            for p in $PS; do
                leaf "parts-$precision-p$p" --ps "$p" --khat "$KHAT" --precision "$precision"
            done
        done
        ;;
    *)
        echo "ch4-9a2b-runs.sh: unknown stage $stage" >&2
        exit 2
        ;;
    esac
done
