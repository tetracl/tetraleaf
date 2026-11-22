#!/bin/bash
# tools/helmholtz/ch4-3c3-runs.sh: the machine runs of Helmholtz CH4.3c3
# (docs/helmholtz/tasks/CH4.3c3.md, "Machine runs"), from the repository root after the
# examples were built in release:
#
#     cargo build --release -p tetraleaf-validate --example helmholtz_calibrate \
#         --example helmholtz_boxes
#     tools/helmholtz/ch4-3c3-runs.sh <stage> <out-dir>
#
# Stages, each launch one example on one rank × THREADS threads (default: the machine's
# cores, 6 on the M6, 40 on a Kathleen node):
# - `scans`: level scans (`helmholtz_calibrate --scan levels`, CH3.3's and CH4.3c2's method:
#   one level's degree from ⌊x⌋ − 2 upward, every other level held at CH3.3's hold rule at
#   10⁻¹⁰, Rotation, ε = 10⁻³ … 10⁻⁸ at ε/2) on W1 with the constant density at every level of
#   x = 13.5 … 216.5, each a launch of its own: kD = 250 levels 2–4 (`w1c-250`, one launch),
#   kD = 500 levels 2–5 (`w1c-500-l<l>`), kD = 1000 levels 2–6 (`w1c-1000-l<l>`; on a
#   Kathleen node about 20 s per degree at kD = 500 and 110 s at kD = 1000, so kD = 1000's
#   level 2 takes about 70 min); FIRST_<name> (the name's dashes as underscores) names a
#   launch's first degree; `scans-low` the same below the plan's range, at x = 6.8 (kD = 250
#   level 5, kD = 500 level 6), where a refit that lowers the degrees reaches;
# - `verify-w1`: W1 with the constant density at kD = 250, 500 and 1000 (`helmholtz_boxes`, the
#   samples only) at 10⁻⁸, 10⁻⁷, 10⁻⁶ and 10⁻³, Rotation and `helmholtz_boxes`' default, under the
#   rule as built (RULE, default `calibrated`; or a rule by its constants
#   `γ,α,β,C,a,b,d[,cap]`), on today's trees (c_leaf 4);
# - `verify-matrix`: CH4.3c2's verification of a changed rule, CH3.1c's matrix and CH3.3's
#   held-out cases (`helmholtz_accuracy`) under Rotation (the near rules off) and the library's
#   defaults (`--m2l auto-pw --near-v cost --near-wx cost`), against ε, under the rule as built
#   and CH4.3c2's (`--rule ch4-3c2`, the errors and times before);
# - `verify-w1-all`: W1 at kD = 250, 500 and 1000 (`helmholtz_boxes`, the samples only) at 10⁻⁸,
#   10⁻⁷, 10⁻⁶ and 10⁻³, the random charges, the plane wave's trace and the constant density,
#   Rotation and `helmholtz_boxes`' default, under the rule as built and CH4.3c2's;
# - `times`: W1 at kD = 1000 and 500 (`helmholtz_boxes`, the random charges), each case twice
#   under each rule, the rules alternated (the rule, CH4.3c2's, the rule, CH4.3c2's), Rotation
#   and the default, as CH4.3c2's report §7.
# Each launch's output in <out-dir>/<name>.md (stderr in <name>.err), a line per launch in
# <out-dir>/runs.log. BLAS threads 1. Accuracy runs; times are read on Kathleen only
# (exclusive nodes), never asserted. TIMEOUT (default 21000) bounds a launch; ONLY (a
# regular expression on the launch names) runs a subset.
set -u
stage=$1
out=$2
mkdir -p "$out"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
E=target/release/examples
TIMEOUT=${TIMEOUT:-21000}
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

# `scan <name> <kd> <level> <eps> [first]`: one level's scan, its first degree FIRST_<name>
# if set, else the argument, else the example's ⌊x⌋ − 2.
scan() {
    local name=$1 kd=$2 level=$3 eps=$4 first=${5:-}
    local var=FIRST_${name//-/_}
    first=${!var:-$first}
    local F=()
    [ -n "$first" ] && F=(--first "$first")
    launch "$name" helmholtz_calibrate --workload w1 --density constant --kd "$kd" \
        --eps "$eps" --scan levels --hold 0.29,0.45,1.65,1 --hold-eps 1e-10 \
        --level "$level" "${F[@]}"
}

ALL=1e-3,1e-4,1e-5,1e-6,1e-7,1e-8
case $stage in
scans)
    # kD = 250, levels 2–4 in one launch (x = 54.1, 27.1, 13.5; level 5 at x = 6.8 is the
    # plateau's).
    if [[ w1c-250 =~ ${ONLY:-.} ]]; then
        launch w1c-250 helmholtz_calibrate --workload w1 --density constant --kd 250 \
            --eps "$ALL" --scan levels --hold 0.29,0.45,1.65,1 --hold-eps 1e-10 --min-x 13
    fi
    for l in 5 4 3 2; do
        scan "w1c-500-l$l" 500 "$l" "$ALL"
    done
    for l in 6 5 4 3 2; do
        scan "w1c-1000-l$l" 1000 "$l" "$ALL"
    done
    ;;
scans-low)
    scan w1c-250-l5 250 5 "$ALL"
    scan w1c-500-l6 500 6 "$ALL"
    ;;
verify-w1)
    rule=${RULE:-calibrated}
    for kd in 250 500 1000; do
        launch "w1-$kd-constant" helmholtz_boxes --workload w1 --kd "$kd" \
            --eps 1e-8,1e-7,1e-6,1e-3 --m2l rotation,default --density constant \
            --per-box 0 --rule "$rule"
    done
    ;;
verify-matrix)
    A=helmholtz_accuracy
    for rule in calibrated ch4-3c2; do
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
verify-w1-all)
    for kd in 250 500 1000; do
        for density in random planewave constant; do
            for rule in calibrated ch4-3c2; do
                launch "w1-$kd-$density-$rule" helmholtz_boxes --workload w1 --kd "$kd" \
                    --eps 1e-8,1e-7,1e-6,1e-3 --m2l rotation,default --density "$density" \
                    --per-box 0 --rule "$rule"
            done
        done
    done
    ;;
times)
    for kd in 1000 500; do
        for i in 1 2; do
            for rule in calibrated ch4-3c2; do
                launch "w1-$kd-$rule-$i" helmholtz_boxes --workload w1 --kd "$kd" \
                    --eps 1e-3,1e-6,1e-8 --m2l rotation,default --per-box 0 --rule "$rule"
            done
        done
    done
    ;;
*)
    echo "usage: $0 scans|scans-low|verify-w1|verify-matrix|verify-w1-all|times <out-dir>" >&2
    exit 2
    ;;
esac
