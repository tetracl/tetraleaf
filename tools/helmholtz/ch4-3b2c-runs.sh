#!/bin/bash
# tools/helmholtz/ch4-3b2c-runs.sh: the machine runs of Helmholtz CH4.3b2c, the diagonal form's
# time rule by default (docs/helmholtz/tasks/CH4.3b2c.md, "Machine runs"), from the repository
# root after
#
#     cargo build --release -p tetraleaf-validate --example helmholtz_accuracy
#     tools/helmholtz/ch4-3b2c-runs.sh <out-dir> <stage…>
#
# Every launch runs `helmholtz_accuracy` with `--levels`, the near rules of the defaults, §9.1's
# complex charges and samples, two evaluations (the rows the second's), the interpolated build
# (the library's default), under the library's defaults (`--m2l auto-pw`, CH4.3c2's degrees,
# CH4.3b2d's `PerBox`), with the diagonal form off (the defaults until here) and under the time
# rule (`--diagonal auto`, the proposed default). `tools/helmholtz/ch4-3b2c-tables.py` makes the
# report's tables. Stages:
#   w1        1 × cores: W1 at kD = 100, 200, 500 at ε 10⁻³, 10⁻⁶, 10⁻⁷, 10⁻⁸ and kD = 1000 at
#             10⁻⁶, `--diagonal off,auto --diagonal-costs machine,form,rcr` (the defaults, the
#             time rule, and the time rule forced to the form and to RCR where a class passes:
#             each option timed at the level's degree, CH4.3b2a)
#   w3        1 × cores: W3 at N = 10⁶, kD = 200 and 500, the same ε and settings
#   reverse   the w1 and w3 cases with `--diagonal auto,off` and the machine's costs alone (a
#             run's place in the launch against the setting)
#   repeat    W3 at kD = 200 at the four ε, `--diagonal off,auto,off,auto` (spread: the two `off`
#             runs; placement: `off` after `auto` against `off` before; CH4.3b2a's report §5.2)
#   ranks     cores × 1 (one rank per core): W1 and W3 at kD = 500 at 10⁻⁶ and 10⁻⁸ with the
#             three costs (the constants at many ranks per node), and on a 40-core node W1 at
#             kD = 1000, 10⁻⁶ (the whole-level fallback of levels 2–3; CH4.1b)
#   smoke     W1 at kD = 30, 10⁻⁷, w1's settings (a check of the script)
# TIMEOUT per launch (default 3600 s), ONLY (a regular expression on the launch names),
# DRY_RUN (print the launches only).
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
slurm=0
[ -n "${SLURM_JOB_ID:-}" ] && slurm=1
TIMEOUT=${TIMEOUT:-3600}
COMMON=(--near-v cost --near-wx cost --levels --m2l auto-pw --evaluations 2
    --diagonal-build interpolated)
THREE=(--diagonal off,auto --diagonal-costs machine,form,rcr)

# `mpi <ranks> <threads>`: the launcher of r ranks of t threads each.
mpi() {
    local r=$1 t=$2
    if [ "$(uname)" = Darwin ]; then
        echo "mpirun --mca btl_tcp_if_include lo0 --mca oob_tcp_if_include lo0 --bind-to none --host localhost:$cores --map-by slot -n $r"
    elif [ "$slurm" = 1 ] && [ "$r" -gt 1 ] && [ $((r % 2)) -eq 0 ]; then
        echo "mpirun --report-bindings --map-by ppr:$((r / 2)):socket:PE=$t --bind-to core -n $r"
    else
        echo "mpirun --report-bindings --map-by slot:PE=$t --bind-to core -n $r"
    fi
}

# `launch <name> <ranks> <threads> <arguments…>`: one launch of the example under the timeout.
launch() {
    local n=$1 r=$2 t=$3
    shift 3
    [[ "$n" =~ ${ONLY:-.} ]] || return 0
    local launcher
    launcher=$(mpi "$r" "$t")
    if [ -n "${DRY_RUN:-}" ]; then
        echo "$n: $launcher $A ${COMMON[*]} --threads $t $*"
        return 0
    fi
    local start=$SECONDS
    {
        echo "# $(date '+%Y-%m-%d %H:%M:%S') $(hostname) job ${SLURM_JOB_ID:-none}"
        echo "# $launcher $A ${COMMON[*]} --threads $t $*"
    } > "$out/$n.err"
    # shellcheck disable=SC2086
    timeout "$TIMEOUT" $launcher "$A" "${COMMON[@]}" --threads "$t" "$@" < /dev/null \
        > "$out/$n.md" 2>> "$out/$n.err"
    local status=$?
    echo "$n: exit $status, $((SECONDS - start)) s, $(date +%T)" | tee -a "$out/runs.log"
}

EPS="1e-3 1e-6 1e-7 1e-8"
# `cases <prefix> <settings…>`: the w1 and w3 cases at one rank with every core.
w1() {
    local p=$1 kd eps
    shift
    for kd in 100 200 500; do
        for eps in $EPS; do
            launch "$p-w1-$kd-$eps" 1 "$cores" --workload w1 --kd "$kd" --eps "$eps" "$@"
        done
    done
    launch "$p-w1-1000-1e-6" 1 "$cores" --workload w1 --kd 1000 --eps 1e-6 "$@"
}
w3() {
    local p=$1 kd eps
    shift
    for kd in 200 500; do
        for eps in $EPS; do
            launch "$p-w3-$kd-$eps" 1 "$cores" --workload w3:1000000 --kd "$kd" --eps "$eps" "$@"
        done
    done
}

for stage in "$@"; do
    case $stage in
    smoke) launch smoke 1 "$cores" --workload w1 --kd 30 --eps 1e-7 "${THREE[@]}" ;;
    w1) w1 run "${THREE[@]}" ;;
    w3) w3 run "${THREE[@]}" ;;
    reverse)
        w1 reverse --diagonal auto,off
        w3 reverse --diagonal auto,off
        ;;
    repeat)
        for eps in $EPS; do
            launch "repeat-w3-200-$eps" 1 "$cores" --workload w3:1000000 --kd 200 \
                --eps "$eps" --diagonal off,auto,off,auto
        done
        ;;
    ranks)
        for eps in 1e-6 1e-8; do
            launch "ranks-w1-500-$eps" "$cores" 1 --workload w1 --kd 500 --eps "$eps" \
                "${THREE[@]}"
            launch "ranks-w3-500-$eps" "$cores" 1 --workload w3:1000000 --kd 500 \
                --eps "$eps" "${THREE[@]}"
        done
        if [ "$cores" -eq 40 ]; then
            launch "ranks-w1-1000-1e-6" 40 1 --workload w1 --kd 1000 --eps 1e-6 \
                --diagonal off,auto
        fi
        ;;
    *)
        echo "ch4-3b2c-runs.sh: unknown stage $stage" >&2
        exit 2
        ;;
    esac
done
