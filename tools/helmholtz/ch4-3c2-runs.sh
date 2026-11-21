#!/bin/bash
# tools/helmholtz/ch4-3c2-runs.sh: the machine runs of Helmholtz CH4.3c2
# (docs/helmholtz/tasks/CH4.3c2.md, "Machine runs"), from the repository root after the
# examples were built in release:
#
#     cargo build --release -p tetraleaf-validate --example helmholtz_calibrate \
#         --example helmholtz_accuracy --example helmholtz_boxes
#     tools/helmholtz/ch4-3c2-runs.sh <stage> <out-dir>
#
# Stages, each launch one example on one rank × THREADS threads (default: the machine's
# cores, 6 on the M6, 40 on a Kathleen node):
# - `scans-cube`: level scans (`helmholtz_calibrate --scan levels`, CH3.3's method: one
#   level's degree from ⌊x⌋ − 2 upward, every other level held at CH3.3's hold rule at
#   10⁻¹⁰, Rotation, ε = 10⁻³ … 10⁻⁸ at ε/2) on W4's cube at N = 10⁵: kD = 250 level 2
#   (x = 54.1), kD = 500 levels 2–4 (x = 108.3, 54.1, 27.1), kD = 1000 levels 2–4
#   (x = 216.5, 108.3, 54.1);
# - `scans-w1`: the same on W1 at kD = 250 level 2 (x = 54.1) and kD = 500 levels 2 and 3
#   (x = 108.3, 54.1);
# - `verify-matrix`: CH3.1c's matrix and CH3.3's held-out cases (`helmholtz_accuracy`)
#   under Rotation (`--m2l rotation`, the near rules off, as CH3.3 verified) and the
#   library's defaults (`--m2l auto-pw --near-v cost --near-wx cost`, as CH4.3b2*), against
#   ε, the rule as built (the refit), and with `--rule ch3-3` (CH3.3's rule) for the errors
#   and times before;
# - `verify-w1`: W1 at kD = 500 and 1000 (`helmholtz_boxes`, the samples only) at 10⁻⁸,
#   10⁻⁷, 10⁻⁶ and 10⁻³, the random charges, the plane wave's trace and the constant
#   density, Rotation and `helmholtz_boxes`' default (Auto, the near rules off, as CH4.0c),
#   the rule as built and CH3.3's.
# Each launch's output in <out-dir>/<name>.md (stderr in <name>.err), a line per launch in
# <out-dir>/runs.log. BLAS threads 1. Accuracy runs; times are read on Kathleen only
# (exclusive nodes), never asserted. TIMEOUT (default 18000) bounds a launch; ONLY (a
# regular expression on the launch names) runs a subset.
set -u
stage=$1
out=$2
mkdir -p "$out"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
E=target/release/examples
TIMEOUT=${TIMEOUT:-18000}
if [ "$(uname)" = Darwin ]; then
    export HWLOC_COMPONENTS=${HWLOC_COMPONENTS:--darwin}
    T=${THREADS:-6}
    pre=(mpirun -n 1 --mca btl_tcp_if_include lo0 --mca oob_tcp_if_include lo0
        --host localhost:6 --bind-to none)
else
    T=${THREADS:-40}
    pre=()
fi

# `launch <name> <example> <arguments…>`: one launch under the timeout.
launch() {
    local name=$1 example=$2
    shift 2
    [[ "$name" =~ ${ONLY:-.} ]] || return 0
    local start=$SECONDS
    timeout "$TIMEOUT" "${pre[@]}" "$E/$example" "$@" --threads "$T" < /dev/null \
        > "$out/$name.md" 2> "$out/$name.err"
    echo "$name: exit $?, $((SECONDS - start)) s, $(date +%T)" | tee -a "$out/runs.log"
}

SCAN=(--eps 1e-3,1e-4,1e-5,1e-6,1e-7,1e-8 --scan levels --hold 0.29,0.45,1.65,1
    --hold-eps 1e-10)
case $stage in
scans-cube)
    C=(--workload w4-cube:100000 "${SCAN[@]}")
    launch cube-250 helmholtz_calibrate "${C[@]}" --kd 250 --level 2
    launch cube-500 helmholtz_calibrate "${C[@]}" --kd 500 --min-x 20
    launch cube-1000 helmholtz_calibrate "${C[@]}" --kd 1000 --min-x 50
    ;;
scans-w1)
    launch w1-250 helmholtz_calibrate --workload w1 "${SCAN[@]}" --kd 250 --level 2
    launch w1-500 helmholtz_calibrate --workload w1 "${SCAN[@]}" --kd 500 --min-x 50
    ;;
verify-matrix)
    A=helmholtz_accuracy
    for rule in calibrated ch3-3; do
        for m2l in rotation defaults; do
            if [ $m2l = rotation ]; then
                R=(--rule "$rule" --m2l rotation)
            else
                R=(--rule "$rule" --m2l auto-pw --near-v cost --near-wx cost)
            fi
            n=$rule-$m2l
            # CH3.1c's matrix (its §6) with W3 at kD = 100 also at 10⁻⁸, as CH3.3 verified.
            launch "$n-w1" $A --workload w1 --kd 10,30,100,500 --eps 1e-3,1e-6,1e-8 "${R[@]}"
            launch "$n-w3" $A --workload w3 --kd 10,30,100 --eps 1e-3,1e-6,1e-8 "${R[@]}"
            for d in cube sphere plummer clusters; do
                launch "$n-w4-$d" $A --workload "w4-$d:100000" --kd 1,10,100 \
                    --eps 1e-3,1e-6,1e-8 "${R[@]}"
            done
            launch "$n-w1-real" $A --workload w1 --kd 100 --eps 1e-3,1e-6 --charges real \
                "${R[@]}"
            launch "$n-w1-complex-k" $A --workload w1 --kd 100 --kim 10 --eps 1e-3,1e-6 \
                "${R[@]}"
            # CH3.3's held-out cases.
            launch "$n-w2" $A --workload w2 --kd 30,100 --eps 1e-3,1e-6,1e-8 "${R[@]}"
            launch "$n-w1-held" $A --workload w1 --kd 50,150 --eps 1e-3,1e-6,1e-8 "${R[@]}"
            for d in cube plummer; do
                launch "$n-w4-$d-held" $A --workload "w4-$d:100000" --kd 3,30 \
                    --eps 1e-3,1e-6,1e-8 "${R[@]}"
            done
            launch "$n-w1-mid" $A --workload w1 --kd 100 --eps 1e-4,1e-7 "${R[@]}"
        done
    done
    ;;
verify-w1)
    for kd in 500 1000; do
        for density in random planewave constant; do
            for rule in calibrated ch3-3; do
                launch "w1-$kd-$density-$rule" helmholtz_boxes --workload w1 --kd "$kd" \
                    --eps 1e-8,1e-7,1e-6,1e-3 --m2l rotation,default --density "$density" \
                    --per-box 0 --rule "$rule"
            done
        done
    done
    ;;
*)
    echo "usage: $0 scans-cube|scans-w1|verify-matrix|verify-w1 <out-dir>" >&2
    exit 2
    ;;
esac
