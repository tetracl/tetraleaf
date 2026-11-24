#!/bin/bash
# tools/helmholtz/ch4-8a-runs.sh: the machine runs of Helmholtz CH4.8a, the constant density on W2
# mapped and a remedy chosen (docs/helmholtz/tasks/CH4.8a.md, "Machine runs"), from the
# repository root after
#
#     cargo build --release -p tetraleaf-validate --example helmholtz_boxes
#     tools/helmholtz/ch4-8a-runs.sh <out-dir> <stage…>
#
# On Kathleen inside a job (tools/kathleen/README.md), one node, the slow nodes excluded, e.g.
#
#     tools/kathleen/submit.sh --time=06:00:00 --exclude=node-c11b-004,node-c11b-013 -- \
#         run.sbatch 'cargo build --release -p tetraleaf-validate --example helmholtz_boxes' \
#         'tools/helmholtz/ch4-8a-runs.sh /scratch/scratch/ucahtbe/fmm/logs/ch4-8a map all survey'
#
# Every launch runs `helmholtz_boxes` at one rank with every core under the library's defaults
# (`--m2l defaults`), one ε, the reference at every point to N = 2·10⁵ and above at §9.1's samples,
# two targets per box of the deepest level and every point of the 8 worst boxes of φ and of ∇φ
# (`--refine 8`) and of the boxes at the centres of the domain's faces (`--refine-faces`), the
# sample cached in CACHE (default <out-dir>/cache; the H4 gate's
# /scratch/scratch/ucahtbe/fmm/logs/h4-gate/cache holds W1's and W2's at kD = 250 and 500).
# Stages:
#   smoke    W2 at kD = 30, the constant density, 10⁻⁶ (a check of the script)
#   map      W2 with the constant density at kD = 250 and 500: `--study levels,m2l --all-levels
#            --raise-by 1,2` (each level's part, each non-RCR level on RCR, the form's classes by
#            decade, the level of the worst ∇φ box's largest part raised by 1 and 2)
#   all      W2 at kD = 250 and 500 and W1 at kD = 500, the constant density: every level with V
#            pairs raised at once by 1 and 2 (`--study single --raise-levels v --raise-together`)
#   survey   W3 with the constant density, resolved at kD = 100 and at N = 10⁶ for kD = 200 and
#            500: `--study levels,m2l --all-levels`
# EPS (default "1e-3 1e-6 1e-8"), TIMEOUT per launch (default 5400 s), ONLY (a regular expression
# on the launch names), DRY_RUN (print the launches only). Each launch writes <out-dir>/<name>.md
# and .err and a line in <out-dir>/runs.log.
set -u
out=$1
shift
CACHE=${CACHE:-$out/cache}
mkdir -p "$out" "$CACHE"
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
COMMON=(--m2l defaults --refine 8 --refine-faces --worst 8 --reference-cache "$CACHE")

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

# `each <prefix> <workload> <kd> <arguments…>`: one launch per ε of EPS, the constant density.
each() {
    local prefix=$1 w=$2 kd=$3
    shift 3
    local e
    for e in $EPS; do
        launch "$prefix-$w-kd$kd-constant-$e" --workload "$w" --kd "$kd" --density constant \
            --eps "$e" "$@"
    done
}

STUDY=(--study levels,m2l --all-levels)
for stage in "$@"; do
    echo "=== $(date '+%Y-%m-%d %H:%M:%S') stage $stage"
    case "$stage" in
    smoke) EPS=1e-6 each smoke w2 30 "${STUDY[@]}" ;;
    map)
        for kd in 250 500; do
            each map w2 "$kd" "${STUDY[@]}" --raise-by 1,2
        done
        ;;
    all)
        for wk in w2:250 w2:500 w1:500; do
            each all "${wk%:*}" "${wk#*:}" --study single --raise-levels v --raise-together \
                --raise-by 1,2
        done
        ;;
    survey)
        each survey w3 100 "${STUDY[@]}"
        each survey w3:1000000 200 "${STUDY[@]}"
        each survey w3:1000000 500 "${STUDY[@]}"
        ;;
    *)
        echo "ch4-8a-runs.sh: unknown stage $stage (smoke, map, all, survey)" >&2
        exit 2
        ;;
    esac
done
