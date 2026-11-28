#!/bin/bash
# tools/helmholtz/ch2-1b-runs.sh: the machine runs of Helmholtz CH2.1b, the Yukawa operator and
# FMM (docs/helmholtz/tasks/CH2.1b.md, "Machine runs"), from the repository root after
#
#     cargo build --release -p tetraleaf-validate --example helmholtz_accuracy
#     tools/helmholtz/ch2-1b-runs.sh <out-dir> <stage…>
#
# Every launch runs `helmholtz_accuracy` at one rank with every core (1 x 40 on a Kathleen node)
# unless the stage says otherwise, with gradients and §9.1's samples, under the library's
# defaults named (`--m2l auto-pw --diagonal auto --diagonal-build interpolated`, the V and W/X
# rules of the cost model; on Yukawa Auto runs Rotation and no level takes the waves or the
# diagonal form). Yukawa is `--kd 0 --kim κ` (D = 1, so κD is κ). `tools/helmholtz/ch2-1b-tables.py`
# makes the report's tables. Stages:
#   smoke        Yukawa on W1's sphere at 20,000 points, κD = 1, 10⁻⁶ (a check of the script)
#   accuracy     the acceptance: W1 at 10⁵ points (`w1:100000`), W3 at 10⁵ and W4's four
#                distributions at 10⁵, at κD = KAPPAS (default 0.01 0.1 1 10 100 1000) and
#                ε = 10⁻³, 10⁻⁶, 10⁻⁸ in one launch per (workload, κD), real charges; W3 at 10⁶
#                at 10⁻⁶ (one launch per κD)
#   determinism  W1 and W4's Plummer sphere at 10⁵ points, κD = 1 and 100, 10⁻⁶, at 1 and every
#                thread in one launch, and at 4 ranks x a quarter of the cores (their checksums
#                against `accuracy`'s)
#   dense        W4's cube at 10⁵ points, κD = 1 and 30, 10⁻⁶: Rotation and Dense in one launch
#   complex      W1 at 10⁵ points, κD = 1 and 100, 10⁻⁶, complex charges
#   checksums    real and complex k, `main` against the branch (one launch each): W1 at kD = 30
#                and 100, W3 at 10⁵ at kD = 30, W4's Plummer sphere at 10⁵ at kD = 10, the cube
#                at kD = 20 + 2i (W4 cube at 10⁵), at 10⁻⁶, under the defaults and under Rotation:
#                every checksum must agree (`main`'s binary from MAIN, a directory holding its
#                `helmholtz_accuracy`)
#   time         D4's evidence: the complex path, `main` against the branch, one setting per
#                launch, in the order main, branch, branch, main, twice: W1 at kD = 100 and 250
#                (real k), W3 at 10⁶ at kD = 50 + 5i (complex k), 10⁻⁶, two evaluations (the row
#                the second's); build and evaluation of the slowest rank
# TIMEOUT per launch (default 3600 s), ONLY (a regular expression on the launch names), DRY_RUN
# (print the launches only). Each launch writes <out-dir>/<name>.md and .err and a line in
# <out-dir>/runs.log. Times are read on Kathleen's exclusive nodes only, never asserted; no run
# on the M6 without the user's go-ahead.
set -u
out=$1
shift
mkdir -p "$out"
out=$(cd "$out" && pwd)
target=${CARGO_TARGET_DIR:-target}/release/examples
B=$target/helmholtz_accuracy
A=${MAIN:+$MAIN/helmholtz_accuracy}
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

# `mpi <ranks>`: the launcher of r ranks of cores/r threads each.
mpi() {
    local r=$1
    if [ "$(uname)" = Darwin ]; then
        echo "mpirun --mca btl_tcp_if_include lo0 --mca oob_tcp_if_include lo0 --bind-to none --host localhost:$cores --map-by slot -n $r"
    else
        echo "mpirun --report-bindings --map-by slot:PE=$((cores / r)) --bind-to core -n $r"
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

# `launch <name> <binary> <arguments…>`: one launch at one rank with every core.
launch() {
    local n=$1 bin=$2
    shift 2
    # shellcheck disable=SC2046
    record "$n" $(mpi 1) "$bin" "${COMMON[@]}" "$@"
}

# `yukawa <name> <workload> <κ> <arguments…>`: a Yukawa launch of the branch at every core.
yukawa() {
    local n=$1 w=$2 kappa=$3
    shift 3
    launch "$n" "$B" --workload "$w" --kd 0 --kim "$kappa" --threads "$cores" "$@"
}

# `need_main <stage>`: a refusal of a comparison of `main` and the branch without MAIN.
need_main() {
    if [ -z "$A" ] && [ -z "${DRY_RUN:-}" ]; then
        echo "ch2-1b-runs.sh: the stage $1 needs MAIN (the directory of main's helmholtz_accuracy)" >&2
        exit 2
    fi
}

for stage in "$@"; do
    case $stage in
    smoke)
        yukawa smoke w1:20000 1 --eps 1e-6 --charges real
        ;;
    accuracy)
        for kappa in $KAPPAS; do
            yukawa "acc-w1-$kappa" "w1:$N" "$kappa" --eps 1e-3,1e-6,1e-8 --charges real
            yukawa "acc-w3-$kappa" "w3:$N" "$kappa" --eps 1e-3,1e-6,1e-8 --charges real
            for d in $W4S; do
                yukawa "acc-w4-$d-$kappa" "w4-$d:$N" "$kappa" --eps 1e-3,1e-6,1e-8 --charges real
            done
        done
        for kappa in $KAPPAS; do
            yukawa "acc-w3large-$kappa" w3:1000000 "$kappa" --eps 1e-6 --charges real
        done
        ;;
    determinism)
        for kappa in 1 100; do
            for w in "w1:$N" "w4-plummer:$N"; do
                tag=${w%%:*}
                launch "det-$tag-$kappa-1x" "$B" --workload "$w" --kd 0 --kim "$kappa" \
                    --eps 1e-6 --charges real --threads "1,$cores"
                # shellcheck disable=SC2046
                record "det-$tag-$kappa-4x" $(mpi 4) "$B" "${COMMON[@]}" --workload "$w" \
                    --kd 0 --kim "$kappa" --eps 1e-6 --charges real --threads $((cores / 4))
            done
        done
        ;;
    dense)
        for kappa in 1 30; do
            yukawa "dense-w4-cube-$kappa" "w4-cube:$N" "$kappa" --eps 1e-6 --charges real \
                --m2l rotation,dense
        done
        ;;
    complex)
        for kappa in 1 100; do
            yukawa "complex-w1-$kappa" "w1:$N" "$kappa" --eps 1e-6 --charges complex
        done
        ;;
    checksums)
        need_main checksums
        for c in w1@30 w1@100 "w3:$N@30" "w4-plummer:$N@10" "w4-cube:$N@20+2"; do
            w=${c%@*} k=${c#*@}
            re=${k%+*} im=0
            [ "$re" != "$k" ] && im=${k#*+}
            tag=$(echo "${w%%:*}-$re-$im" | tr -d '.')
            for m2l in auto-pw rotation; do
                for who in main branch; do
                    bin=$B
                    [ "$who" = main ] && bin=$A
                    launch "sum-$tag-$m2l-$who" "$bin" --workload "$w" --kd "$re" --kim "$im" \
                        --eps 1e-6 --threads "$cores" --m2l "$m2l"
                done
            done
        done
        ;;
    time)
        need_main time
        for c in w1@100 w1@250 w3:1000000@50+5; do
            w=${c%@*} k=${c#*@}
            re=${k%+*} im=0
            [ "$re" != "$k" ] && im=${k#*+}
            tag=$(echo "${w%%:*}-$re-$im" | tr -d '.')
            for i in 1 2; do
                for o in "main:$A:a" "branch:$B:a" "branch:$B:b" "main:$A:b"; do
                    IFS=: read -r who bin j <<< "$o"
                    launch "time-$tag-$who-$i$j" "$bin" --workload "$w" --kd "$re" --kim "$im" \
                        --eps 1e-6 --threads "$cores" --evaluations 2
                done
            done
        done
        ;;
    *)
        echo "ch2-1b-runs.sh: unknown stage $stage" >&2
        exit 2
        ;;
    esac
done
