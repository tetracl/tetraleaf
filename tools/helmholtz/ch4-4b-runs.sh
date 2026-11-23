#!/bin/bash
# tools/helmholtz/ch4-4b-runs.sh: the machine runs of Helmholtz CH4.4b, the crosscheck mapped to
# levels and classes (docs/helmholtz/tasks/CH4.4b.md, "Machine runs"), from the repository root
# after
#
#     cargo build --release -p tetraleaf-validate --example helmholtz_boxes
#     tools/helmholtz/ch4-4b-runs.sh <out-dir> <stage…>
#
# On Kathleen inside a job (tools/kathleen/README.md), one node, the slow nodes excluded, e.g.
#
#     tools/kathleen/submit.sh --time=06:00:00 --exclude=node-c11b-004,node-c11b-013 -- \
#         run.sbatch 'cargo build --release -p tetraleaf-validate --example helmholtz_boxes' \
#         'tools/helmholtz/ch4-4b-runs.sh /scratch/scratch/ucahtbe/fmm/logs/ch4-4b low'
#
# Every launch runs `helmholtz_boxes` at one rank with every core under the library's defaults
# (`--m2l defaults`), one ε, `--study levels,m2l` (each level raised and on RCR by an override;
# each level whose M2L is not RCR on RCR; the form's classes by decade below the bound's threshold
# on RCR), the reference at every point to N = 2·10⁵ and above at §9.1's samples, two targets per
# box of the deepest level and every point of the 8 worst boxes of φ and of ∇φ (`--refine 8`) and of
# the boxes at the centres of the domain's faces, W1's tangent points (`--refine-faces`), the sample
# cached for the case's other ε (<out-dir>/cache). The CH4.4a cases whose trees run the
# diagonal form or RCR at k a ≥ 2.8 (its logs, "by level"). Stages:
#   smoke    W1 at kD = 30, 10⁻⁶ (a check of the script)
#   low      W1, W2 at kD = 30 (10⁻⁶), 100, 250; the constant density on W1, W2 at kD = 30
#            (10⁻⁶) and 100; W3 resolved at kD = 30 (10⁻⁶)
#   mid      W1 at kD = 500 (random charges and the plane wave's trace), W2 at kD = 500
#   w3       W3 resolved at kD = 100; W3 at N = 10⁶ for kD = 200 and 500
#   w1-1000  W1 at kD = 1000, random charges
#   pw1000   W1 at kD = 1000, the plane wave's trace
#   remedy   W1 with the plane wave's trace at kD = 500 and 1000, 10⁻⁸: also each class of the
#            threshold's decade alone (`--class-runs`) and the level of the worst ∇φ box raised
#            by 1 and 2 under the defaults' M2L (`--raise-by 1,2`); then every level with V pairs
#            raised at once by 1, 2 and 4 (`--raise-levels v --raise-together`)
# EPS (default "1e-3 1e-6 1e-8"), TIMEOUT per launch (default 5400 s), ONLY (a regular expression
# on the launch names), DRY_RUN (print the launches only). Each launch writes <out-dir>/<name>.md
# and .err and a line in <out-dir>/runs.log.
set -u
out=$1
shift
mkdir -p "$out/cache"
target=${CARGO_TARGET_DIR:-target}/release/examples
B=$target/helmholtz_boxes
if [ "$(uname)" = Darwin ]; then
    export HWLOC_COMPONENTS=${HWLOC_COMPONENTS--darwin}
    cores=$(sysctl -n hw.perflevel0.physicalcpu 2> /dev/null || echo 6)
else
    # GNU nproc honours OMP_NUM_THREADS: count the cores without it.
    cores=$(env -u OMP_NUM_THREADS nproc)
fi
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
TIMEOUT=${TIMEOUT:-5400}
EPS=${EPS:-1e-3 1e-6 1e-8}
COMMON=(--m2l defaults --study levels,m2l --refine 8 --refine-faces --worst 8
    --reference-cache "$out/cache")

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
        echo "$n: $launcher $B ${COMMON[*]} --threads $cores $*"
        return 0
    fi
    local start=$SECONDS
    {
        echo "# $(date '+%Y-%m-%d %H:%M:%S') $(hostname) job ${SLURM_JOB_ID:-none}"
        echo "# $launcher $B ${COMMON[*]} --threads $cores $*"
    } > "$out/$n.err"
    # shellcheck disable=SC2086
    timeout "$TIMEOUT" $launcher "$B" "${COMMON[@]}" --threads "$cores" "$@" < /dev/null \
        > "$out/$n.md" 2>> "$out/$n.err"
    local status=$?
    echo "$n: exit $status, $((SECONDS - start)) s, $(date +%T)" | tee -a "$out/runs.log"
}

# `each <workload> <kd> <density> <eps…>`: one launch per ε.
each() {
    local w=$1 kd=$2 density=$3
    shift 3
    local e
    for e in "$@"; do
        launch "$w-kd$kd-$density-$e" --workload "$w" --kd "$kd" --density "$density" --eps "$e"
    done
}

for stage in "$@"; do
    echo "=== $(date '+%Y-%m-%d %H:%M:%S') stage $stage"
    case "$stage" in
    smoke) each w1 30 random 1e-6 ;;
    low)
        # At kD = 30 the far field is Batched at 10⁻³ and the waves at 10⁻⁸ (CH4.4a's logs).
        for w in w1 w2; do
            each "$w" 30 random 1e-6
            # shellcheck disable=SC2086
            each "$w" 100 random $EPS
            # shellcheck disable=SC2086
            each "$w" 250 random $EPS
            each "$w" 30 constant 1e-6
            # shellcheck disable=SC2086
            each "$w" 100 constant $EPS
        done
        each w3 30 random 1e-6
        ;;
    mid)
        # shellcheck disable=SC2086
        each w1 500 random $EPS
        # shellcheck disable=SC2086
        each w1 500 planewave $EPS
        # shellcheck disable=SC2086
        each w2 500 random $EPS
        ;;
    w3)
        # shellcheck disable=SC2086
        each w3 100 random $EPS
        # shellcheck disable=SC2086
        each w3:1000000 200 random $EPS
        # shellcheck disable=SC2086
        each w3:1000000 500 random $EPS
        ;;
    # shellcheck disable=SC2086
    w1-1000) each w1 1000 random $EPS ;;
    # shellcheck disable=SC2086
    pw1000) each w1 1000 planewave $EPS ;;
    remedy)
        for kd in 500 1000; do
            launch "remedy-w1-kd$kd-planewave-1e-8" --workload w1 --kd "$kd" --density planewave \
                --eps 1e-8 --class-runs --raise-by 1,2
            # Every level with V pairs raised at once (the worst box's error spreads over the
            # levels at a tangent point, CH4.0b), no study.
            launch "remedy-all-w1-kd$kd-planewave-1e-8" --workload w1 --kd "$kd" \
                --density planewave --eps 1e-8 --study single --raise-levels v --raise-together \
                --raise-by 1,2,4
        done
        ;;
    *)
        echo "ch4-4b-runs.sh: unknown stage $stage (smoke, low, mid, w3, w1-1000, pw1000, remedy)" >&2
        exit 2
        ;;
    esac
done
