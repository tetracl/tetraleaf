#!/bin/bash
# tools/helmholtz/ch4-3e-runs.sh: the machine runs of Helmholtz CH4.3e, c_leaf re-decided with the
# diagonal form's time rule the default (docs/helmholtz/tasks/CH4.3e.md, "Machine runs"), from the
# repository root after
#
#     cargo build --release -p tetraleaf-validate --example helmholtz_accuracy
#     tools/helmholtz/ch4-3e-runs.sh <out-dir> <stage…>
#
# Every launch runs `helmholtz_accuracy` at one rank with every core, `--levels`, §9.1's complex
# charges and samples, two evaluations (the rows the second's), under the library's defaults named
# (`--m2l auto-pw`, the time rule `--diagonal auto` with the interpolated build, the V, W and X
# rules of the cost model, CH4.3b2d's `PerBox`, CH4.3c2's degrees, the leaf-size rule), at c_leaf
# 1, 2, 4 and 8 in one order (`fwd-*`) and in the other (`rev-*`): Kathleen's second run in a
# launch is slower whatever its setting (CH4.3b2c's report §6). `tools/helmholtz/ch4-3e-tables.py`
# makes the report's tables. Stages:
#   w1      W1 at kD = 100, 200, 500 at ε 10⁻³, 10⁻⁶, 10⁻⁷, 10⁻⁸
#   w2      W2 at kD = 100, 200, 500 at the same ε
#   w3      W3 at N = 10⁶, kD = 200 and 500 at the same ε
#   large   W1 at kD = 1000 at 10⁻³ and 10⁻⁶; W3 resolved at kD = 100 (N = 4·10⁶) at 10⁻³, 10⁻⁶,
#           10⁻⁸ (`ONLY` picks one workload)
#   smoke   W2 at kD = 60, 10⁻⁶ (a check of the script)
# C_LEAF (default "1,2,4,8"; `rev-*` runs it reversed), TIMEOUT per launch (default 3600 s), ONLY
# (a regular expression on the launch names), DRY_RUN (print the launches only).
# Each launch writes <out-dir>/<name>.md and .err and a line in <out-dir>/runs.log.
set -u
out=$1
shift
mkdir -p "$out"
target=${CARGO_TARGET_DIR:-target}/release/examples
A=$target/helmholtz_accuracy
if [ "$(uname)" = Darwin ]; then
    export HWLOC_COMPONENTS=${HWLOC_COMPONENTS--darwin}
    cores=$(sysctl -n hw.perflevel0.physicalcpu 2> /dev/null || echo 6)
else
    # GNU nproc honours OMP_NUM_THREADS: count the cores without it.
    cores=$(env -u OMP_NUM_THREADS nproc)
fi
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
TIMEOUT=${TIMEOUT:-3600}
C_LEAF=${C_LEAF:-1,2,4,8}
REVERSED=$(echo "$C_LEAF" | awk -F, '{ for (i = NF; i > 1; i--) printf "%s,", $i; print $1 }')
COMMON=(--near-v cost --near-wx cost --levels --m2l auto-pw --evaluations 2 --diagonal auto
    --diagonal-build interpolated)

# `mpi`: the launcher of one rank with every core.
mpi() {
    if [ "$(uname)" = Darwin ]; then
        echo "mpirun --mca btl_tcp_if_include lo0 --mca oob_tcp_if_include lo0 --bind-to none --host localhost:$cores --map-by slot -n 1"
    else
        echo "mpirun --report-bindings --map-by slot:PE=$cores --bind-to core -n 1"
    fi
}

# `launch <name> <arguments…>`: one launch of the example under the timeout.
launch() {
    local n=$1
    shift
    [[ "$n" =~ ${ONLY:-.} ]] || return 0
    local launcher
    launcher=$(mpi)
    if [ -n "${DRY_RUN:-}" ]; then
        echo "$n: $launcher $A ${COMMON[*]} --threads $cores $*"
        return 0
    fi
    local start=$SECONDS
    {
        echo "# $(date '+%Y-%m-%d %H:%M:%S') $(hostname) job ${SLURM_JOB_ID:-none}"
        echo "# $launcher $A ${COMMON[*]} --threads $cores $*"
    } > "$out/$n.err"
    # shellcheck disable=SC2086
    timeout "$TIMEOUT" $launcher "$A" "${COMMON[@]}" --threads "$cores" "$@" < /dev/null \
        > "$out/$n.md" 2>> "$out/$n.err"
    local status=$?
    echo "$n: exit $status, $((SECONDS - start)) s, $(date +%T)" | tee -a "$out/runs.log"
}

# `both <case> <arguments…>`: the case in both orders of c_leaf.
both() {
    local c=$1
    shift
    launch "fwd-$c" "$@" --c-leaf "$C_LEAF"
    launch "rev-$c" "$@" --c-leaf "$REVERSED"
}

EPS="1e-3 1e-6 1e-7 1e-8"
# `sweep <tag> <workload> <kD…>`: the case at every ε of EPS and each kD.
sweep() {
    local tag=$1 workload=$2 kd eps
    shift 2
    for kd in "$@"; do
        for eps in $EPS; do
            both "$tag-$kd-$eps" --workload "$workload" --kd "$kd" --eps "$eps"
        done
    done
}

for stage in "$@"; do
    case $stage in
    smoke) both smoke --workload w2 --kd 60 --eps 1e-6 ;;
    w1) sweep w1 w1 100 200 500 ;;
    w2) sweep w2 w2 100 200 500 ;;
    w3) sweep w3 w3:1000000 200 500 ;;
    large)
        for eps in 1e-3 1e-6; do
            both "w1-1000-$eps" --workload w1 --kd 1000 --eps "$eps"
        done
        for eps in 1e-3 1e-6 1e-8; do
            both "w3r-100-$eps" --workload w3 --kd 100 --eps "$eps"
        done
        ;;
    *)
        echo "ch4-3e-runs.sh: unknown stage $stage" >&2
        exit 2
        ;;
    esac
done
