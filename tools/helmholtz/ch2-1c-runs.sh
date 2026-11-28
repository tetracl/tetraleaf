#!/bin/bash
# tools/helmholtz/ch2-1c-runs.sh: the machine runs of Helmholtz CH2.1c, Yukawa's batched P2M
# (docs/helmholtz/tasks/CH2.1c.md, "Machine runs"), from the repository root after
#
#     cargo build --release -p tetraleaf-validate --example helmholtz_accuracy
#     tools/helmholtz/ch2-1c-runs.sh <out-dir> <stage…>
#
# Every launch runs `helmholtz_accuracy` at one rank with every core (1 x 40 on a Kathleen node),
# Yukawa (`--kd 0 --kim κ`, D = 1, so κD is κ), with gradients and §9.1's samples, under the
# library's defaults named (`--m2l auto-pw --diagonal auto --diagonal-build interpolated`, the V and
# W/X rules of the cost model; on Yukawa Auto runs Rotation), and one P2M per launch:
# `--leaves kernel` (the batched kernel, the library's default) or `--leaves reference` (the
# per-point reference calls; P2P stays the vectorised kernel). Accuracy may share a launch, times
# may not (review 2026-11-24, B12). `tools/helmholtz/ch2-1c-tables.py` makes the report's tables.
# Stages:
#   smoke     W1's sphere at 20,000 points, κD = 1, 10⁻⁶, the kernel then the reference
#   accuracy  CH2.1b's cases (`tools/helmholtz/ch2-1b-runs.sh accuracy`): W1 at 10⁵ points, W3 at
#             10⁵ and W4's four distributions at 10⁵ at κD = KAPPAS (default 0.01 0.1 1 10 100
#             1000) and ε = 10⁻³, 10⁻⁶, 10⁻⁸ in one launch per (workload, κD, P2M); W3 at 10⁶ at
#             10⁻⁶; real charges; and W1 at 10⁵ at κD = 1 and 100, 10⁻⁶, complex charges
#   worst     the constant density (the worst the record names, CH4.3c3 and CH4.8): W1 and W2 at
#             κD = 1 and 100, 10⁻⁶, real charges, the kernel and the reference
#   levels    P2M's time: W1 (its default N) and W3 at 10⁶, κD = 1 and 100, 10⁻⁶, `--levels`, two
#             evaluations (the rows the second's), per case kernel, reference, reference, kernel
# TIMEOUT per launch (default 3600 s), ONLY (a regular expression on the launch names), DRY_RUN
# (print the launches only). Each launch writes <out-dir>/<name>.md and .err and a line in
# <out-dir>/runs.log. Times are read on Kathleen's exclusive nodes only, never asserted; no run on
# the M6 without the user's go-ahead.
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
KAPPAS=${KAPPAS:-0.01 0.1 1 10 100 1000}
COMMON=(--near-v cost --near-wx cost --m2l auto-pw --diagonal auto --diagonal-build interpolated)
N=100000
W4S="cube sphere plummer clusters"

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

# `yukawa <name> <leaves> <workload> <κ> <arguments…>`: a Yukawa launch at every core.
yukawa() {
    local n=$1 leaves=$2 w=$3 kappa=$4
    shift 4
    # shellcheck disable=SC2046
    record "$n-$leaves" $(mpi) "$A" "${COMMON[@]}" --threads "$cores" --leaves "$leaves" \
        --workload "$w" --kd 0 --kim "$kappa" "$@"
}

# `both <name> <arguments…>`: the kernel, then the reference.
both() {
    local n=$1
    shift
    yukawa "$n" kernel "$@"
    yukawa "$n" reference "$@"
}

for stage in "$@"; do
    case $stage in
    smoke)
        both smoke w1:20000 1 --eps 1e-6 --charges real
        ;;
    accuracy)
        for kappa in $KAPPAS; do
            both "acc-w1-$kappa" "w1:$N" "$kappa" --eps 1e-3,1e-6,1e-8 --charges real
            both "acc-w3-$kappa" "w3:$N" "$kappa" --eps 1e-3,1e-6,1e-8 --charges real
            for d in $W4S; do
                both "acc-w4-$d-$kappa" "w4-$d:$N" "$kappa" --eps 1e-3,1e-6,1e-8 --charges real
            done
            both "acc-w3large-$kappa" w3:1000000 "$kappa" --eps 1e-6 --charges real
        done
        for kappa in 1 100; do
            both "complex-w1-$kappa" "w1:$N" "$kappa" --eps 1e-6 --charges complex
        done
        ;;
    worst)
        for w in w1 w2; do
            for kappa in 1 100; do
                both "worst-$w-$kappa" "$w" "$kappa" --eps 1e-6 --charges real --density constant
            done
        done
        ;;
    levels)
        for c in w1@1 w1@100 w3:1000000@1 w3:1000000@100; do
            w=${c%@*} kappa=${c#*@}
            tag=${w%%:*}-$kappa
            for o in kernel:1 reference:1 reference:2 kernel:2; do
                IFS=: read -r leaves i <<< "$o"
                # shellcheck disable=SC2046
                record "levels-$tag-$leaves-$i" $(mpi) "$A" "${COMMON[@]}" --threads "$cores" \
                    --leaves "$leaves" --workload "$w" --kd 0 --kim "$kappa" --eps 1e-6 \
                    --charges real --levels --evaluations 2
            done
        done
        ;;
    *)
        echo "ch2-1c-runs.sh: unknown stage $stage" >&2
        exit 2
        ;;
    esac
done
