#!/bin/bash
# tools/helmholtz/ch4-8b-runs.sh: the machine runs of Helmholtz CH4.8b, CH4.8a's remedy (one degree
# after the cap of the high-frequency branch) as the default, verified (docs/helmholtz/tasks/CH4.8b.md,
# "Machine runs"), from the repository root after
#
#     cargo build --release -p tetraleaf-validate --example helmholtz_boxes --example helmholtz_fmm3d
#     tools/helmholtz/ch4-8b-runs.sh <out-dir> <stage…>
#
# On Kathleen inside a job (tools/kathleen/README.md), one node, the slow nodes excluded, e.g.
#
#     tools/kathleen/submit.sh '--time=06:00:00 --exclude=node-c11b-004,node-c11b-013' -- \
#         run.sbatch 'cargo build --release -p tetraleaf-validate --example helmholtz_boxes \
#         --example helmholtz_fmm3d' 'tools/helmholtz/ch4-8b-runs.sh <logs>/ch4-8b const'
#
# One rank with every core (1 × 40 on Kathleen, 1 × 72 on locust), BLAS threads 1. Stages:
#   smoke   W2 at kD = 30 with the constant density at 10⁻⁶ (`const`'s launch) and W1 at kD = 30 at
#           10⁻³ and 10⁻⁶ (`random`'s), a check of the script
#   const   `helmholtz_boxes --m2l defaults --study single`, the reference at every point to
#           N = 2·10⁵ and above at §9.1's samples, two targets per box of the deepest level and every
#           point of the 8 worst boxes of φ and of ∇φ (`--refine 8`) and of the boxes at the centres
#           of the domain's faces (`--refine-faces`), the sample cached in CACHE (default
#           <out-dir>/cache; CH4.8a's and the H4 gate's caches may be copied there first), the
#           constant density, one launch per ε of EPS (default the six of 10⁻³ … 10⁻⁸), on CONST
#           (default W1 and W2 at kD = 10, 30, 100, 175, 250, 500, W2 also at 375, W3 resolved at
#           kD = 10, 30, 100, W3 at N = 10⁶ for kD = 200 and 500): the acceptance's constant density
#   before  the same launches under `--rule ch4-3c3` (`DegreeRule::CH4_3C3`, the default before
#           CH4.8b): each row's before (run after `const`, so that the cached samples are shared)
#   random  `helmholtz_fmm3d compare --settings defaults` without FMM3D, every ε of EPS in one
#           launch, on RANDOM (default CH4.4a's random-charge and plane-wave cases): the
#           acceptance's random charges, with the build, the evaluation and the memory
#   cost    `compare --settings defaults --evaluations 2`, one ε of COST_EPS (default 10⁻³, 10⁻⁶,
#           10⁻⁸) per launch, CH4.3c3's rule (A, `--rule ch4-3c3`) and the default (B) in the order
#           A B B A (one setting per launch, CH4.3b2c §6), on COST (default W1 at kD = 100, 250, 500,
#           1000, W2 at 500, W3 at N = 10⁶ for kD = 500)
# The cases are `<workload>@<kD>` with the examples' workload names (`w3:1000000@500`). ONLY (a
# regular expression on the launch names), TIMEOUT per launch (default 5400 s), DRY_RUN (print
# the launches only). Each launch writes <out-dir>/<name>.md and .err and a line in
# <out-dir>/runs.log; off Slurm (locust) the line names the load and the free memory before it.
set -u
out=$1
shift
CACHE=${CACHE:-$out/cache}
mkdir -p "$out" "$CACHE"
target=${CARGO_TARGET_DIR:-target}/release/examples
B=$target/helmholtz_boxes
F=$target/helmholtz_fmm3d
if [ "$(uname)" = Darwin ]; then
    export HWLOC_COMPONENTS=${HWLOC_COMPONENTS--darwin}
    cores=$(sysctl -n hw.perflevel0.physicalcpu 2> /dev/null || echo 6)
else
    # GNU nproc honours OMP_NUM_THREADS: count the cores without it.
    cores=$(env -u OMP_NUM_THREADS nproc)
fi
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
TIMEOUT=${TIMEOUT:-5400}
EPS=${EPS:-1e-3 1e-4 1e-5 1e-6 1e-7 1e-8}
COST_EPS=${COST_EPS:-1e-3 1e-6 1e-8}
const_cases="w1@10 w1@30 w1@100 w1@175 w1@250 w1@500"
const_cases="$const_cases w2@10 w2@30 w2@100 w2@175 w2@250 w2@375 w2@500"
const_cases="$const_cases w3@10 w3@30 w3@100 w3:1000000@200 w3:1000000@500"
CONST=${CONST:-$const_cases}
random_cases="w1@10 w1@30 w1@100 w1@250 w1@500 w1@1000 w2@10 w2@30 w2@100 w2@250 w2@500"
random_cases="$random_cases w3@10 w3@30 w3@100 w3:1000000@200 w3:1000000@500"
random_cases="$random_cases w1@500/planewave w1@1000/planewave"
RANDOM_CASES=${RANDOM_CASES:-$random_cases}
COST=${COST:-w1@100 w1@250 w1@500 w1@1000 w2@500 w3:1000000@500}
BOXES=(--m2l defaults --study single --refine 8 --refine-faces --worst 8 --reference-cache "$CACHE")
if [ -n "${SLURM_JOB_ID:-}" ]; then slurm=1; else slurm=0; fi

# `mpi`: the launcher of one rank with every core.
mpi() {
    if [ "$(uname)" = Darwin ]; then
        echo "mpirun --mca btl_tcp_if_include lo0 --mca oob_tcp_if_include lo0 --bind-to none --host localhost:$cores --map-by slot -n 1"
    else
        echo "mpirun --report-bindings --map-by slot:PE=$cores --bind-to core -n 1"
    fi
}

# `state`: off Slurm, the load and the free memory (locust is shared).
state() {
    [ "$slurm" = 1 ] || [ "$(uname)" = Darwin ] ||
        echo "load $(uptime | sed 's/.*load average: //'), available $(free -g | awk '/Mem:/ { print $7 }') GB"
}

# `parts <case>`: sets workload, kd, density and name (`w3:1000000@500/planewave` ->
# `w3-n1000000-kd500-planewave`).
parts() {
    workload=${1%@*} kd=${1#*@} density=random
    case "$kd" in
    */*) density=${kd#*/} kd=${kd%/*} ;;
    esac
    name=${workload/:/-n}-kd$kd
    [ "$density" = random ] || name=$name-$density
}

# `launch <name> <example> <arguments…>`: one launch of an example under the timeout.
launch() {
    local n=$1 exe=$2
    shift 2
    [[ "$n" =~ ${ONLY:-.} ]] || return 0
    local launcher
    launcher=$(mpi)
    if [ -n "${DRY_RUN:-}" ]; then
        echo "$n: $launcher $exe $* --threads $cores"
        return 0
    fi
    local start=$SECONDS before
    before=$(state)
    {
        echo "# $(date '+%Y-%m-%d %H:%M:%S') $(hostname) job ${SLURM_JOB_ID:-none} $before"
        echo "# $launcher $exe $* --threads $cores"
    } > "$out/$n.err"
    # shellcheck disable=SC2086 # the launcher's words split
    timeout "$TIMEOUT" $launcher "$exe" "$@" --threads "$cores" < /dev/null \
        > "$out/$n.md" 2>> "$out/$n.err"
    local status=$?
    echo "$n: exit $status, $((SECONDS - start)) s, $(date +%T) $before" | tee -a "$out/runs.log"
}

# `boxes <prefix> <case> <arguments…>`: one `helmholtz_boxes` launch per ε of EPS, the constant
# density.
boxes() {
    local prefix=$1 e
    parts "$2"
    shift 2
    for e in $EPS; do
        launch "$prefix-$name-constant-$e" "$B" --workload "$workload" --kd "$kd" \
            --density constant --eps "$e" "${BOXES[@]}" "$@"
    done
}

# `compare <name> <case> <eps,…> <arguments…>`: one `helmholtz_fmm3d compare` launch, ours alone.
compare() {
    local n=$1 e=$3
    parts "$2"
    shift 3
    launch "$n" "$F" compare --workload "$workload" --kd "$kd" --density "$density" \
        --settings defaults --eps "$e" "$@"
}

for stage in "$@"; do
    echo "=== $(date '+%Y-%m-%d %H:%M:%S') stage $stage"
    case "$stage" in
    smoke)
        EPS=1e-6 boxes smoke w2@30
        compare smoke-w1-kd30 w1@30 1e-3,1e-6
        ;;
    const) for c in $CONST; do boxes const "$c"; done ;;
    before) for c in $CONST; do boxes before "$c" --rule ch4-3c3; done ;;
    random)
        for c in $RANDOM_CASES; do
            parts "$c"
            compare "random-$name" "$c" "$(echo "$EPS" | tr ' ' ',')"
        done
        ;;
    cost)
        for c in $COST; do
            parts "$c"
            for e in $COST_EPS; do
                i=0
                for r in A B B A; do
                    i=$((i + 1))
                    if [ "$r" = A ]; then
                        compare "cost-$name-$e-$i-ch4-3c3" "$c" "$e" --evaluations 2 --rule ch4-3c3
                    else
                        compare "cost-$name-$e-$i-default" "$c" "$e" --evaluations 2
                    fi
                done
            done
        done
        ;;
    *)
        echo "ch4-8b-runs.sh: unknown stage $stage (smoke, const, before, random, cost)" >&2
        exit 2
        ;;
    esac
done
