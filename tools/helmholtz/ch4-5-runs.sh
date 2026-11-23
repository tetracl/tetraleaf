#!/bin/bash
# tools/helmholtz/ch4-5-runs.sh: the machine runs of Helmholtz CH4.5, host performance against
# FMM3D at matched accuracy (docs/helmholtz/tasks/CH4.5.md, "Machine runs"), from the repository
# root after
#
#     cargo build --release -p tetraleaf-validate --example helmholtz_fmm3d
#     tools/helmholtz/ch4-5-runs.sh <out-dir> <stage…>
#
# On Kathleen inside a job (tools/kathleen/README.md), one node, slow nodes excluded:
#
#     tools/kathleen/submit.sh '--time=06:00:00 --exclude=node-c11b-004,node-c11b-013' -- \
#         run.sbatch 'tools/helmholtz/ch4-5-runs.sh /scratch/scratch/ucahtbe/fmm/logs/ch4-5 grid'
#
# The cases are CH4.4a's sweep (`<workload>@<kD>[/<density>]`, the examples' workload names);
# CASES (a list) or ONLY (a regular expression on launch names) narrows them. Stages:
#   grid   FMM3D over the eps grid GRID (default 3e-3, 1e-3, ..., 1e-9), each variant of
#          VARIANTS (`default fastker`; `fastker` alone at kD = 1000, CH4.4a's question 3; one
#          variant on the density cases as in CH4.4a, `default`, `fastker` at kD = 1000), REPEAT 1,
#          OpenMP at its defaults, one launch per eps under tools/fmm3d/maxrss.py; then
#          `helmholtz_fmm3d errors` over the grid into grid-<case>.md. Decides the matches.
#   ours   our FMM under the library's defaults (`compare --settings defaults`, no FMM3D), each
#          case at each ε of EPS alone (our memory one ε per launch, CH4.4a's question 4), one
#          build and two evaluations (`--evaluations 2`), on every split of SPLITS in order
#          (launch a) and reversed (b): ours-<case>-<ε>-<split>-<a|b>.md and .err.
#   cleaf  W3 at N = 10⁶, kD = 500 (CH4.3e's Also) at c_leaf 2 and 4, every split of SPLITS
#          and N × 1, a and b: cleaf<c>-w3-n1000000-kd500-<ε>-<split>-<a|b>.md.
#   time   FMM3D at the eps of <out-dir>/time-eps.txt (written by `ch4-5-tables.py eps` from
#          the grid and our errors; lines `<case> <variant> <eps,…>`), OpenMP at its defaults,
#          REPEAT 2 (the best of two calls) below N = 10⁶ points and 1 from there (a call of
#          seconds to minutes, its setup in every call), launch a (the variants in order) and b
#          (reversed), one launch per eps under maxrss.py; then `errors` over them into
#          time-<case>.md (this machine's errors at those eps). `default` at kD = 1000 runs
#          off Slurm only, and only with 250 GB available (Kathleen's node has 180 GB).
#   bind   FMM3D bound (OMP_PLACES=cores OMP_PROC_BIND=close) beside its defaults (review
#          2026-11-19, B8), at the eps of <out-dir>/bind-eps.txt (the same format), each
#          binding in launches a and b in opposite orders, as `time`; into bind-<case>.md.
# SPLITS defaults to "1x40 2x20 4x10 8x5" under Slurm and "1x72 2x36 4x18 8x9" elsewhere; the
# N × 1 split of `cleaf` is the cores'. Under Slurm an even rank count is bound per socket
# (`ppr:<r/2>:socket:PE=<t>`), elsewhere `slot:PE=<t>`. FMM3D_ROOT as for tools/fmm3d/build.sh;
# TIMEOUT per launch (default 7200 s); DRY_RUN prints the launches only. Every launch writes a
# line to <out-dir>/runs.log; off Slurm (locust) also the load and free memory before and
# after it. Times are read on Kathleen's exclusive nodes and a quiet locust only, never
# asserted; no run on the M6 without the user's go-ahead.
set -u
out=$1
shift
mkdir -p "$out/problems" "$out/grid" "$out/time" "$out/bind"
target=${CARGO_TARGET_DIR:-target}/release/examples
F=$target/helmholtz_fmm3d
if [ -z "${FMM3D_ROOT:-}" ]; then
    if [ -n "${TETRALEAF_KATHLEEN_ROOT:-}" ]; then
        FMM3D_ROOT=$TETRALEAF_KATHLEEN_ROOT/fmm3d
    elif [ -d /data/ucahtbe ]; then
        FMM3D_ROOT=/data/ucahtbe/fmm3d
    else
        FMM3D_ROOT=
    fi
fi
# GNU nproc honours OMP_NUM_THREADS: count the cores without it.
cores=$(env -u OMP_NUM_THREADS nproc 2> /dev/null || sysctl -n hw.perflevel0.physicalcpu)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
if [ -n "${SLURM_JOB_ID:-}" ]; then
    slurm=1
    SPLITS=${SPLITS:-"1x40 2x20 4x10 8x5"}
else
    slurm=0
    SPLITS=${SPLITS:-"1x72 2x36 4x18 8x9"}
fi
TIMEOUT=${TIMEOUT:-7200}
EPS=${EPS:-"1e-3 1e-6 1e-8"}
GRID=${GRID:-3e-3,1e-3,3e-4,1e-4,3e-5,1e-5,3e-6,1e-6,3e-7,1e-7,3e-8,1e-8,3e-9,1e-9}
VARIANTS=${VARIANTS:-"default fastker"}
sweep="w1@10 w1@30 w1@100 w1@250 w1@500 w1@1000 w2@10 w2@30 w2@100 w2@250 w2@500"
sweep="$sweep w3@10 w3@30 w3@100 w3:1000000@200 w3:1000000@500 w1@500/planewave w1@1000/planewave"
for c in w1@10 w1@30 w1@100 w2@10 w2@30 w2@100 w3@10 w3@30; do sweep="$sweep $c/constant"; done
CASES=${CASES:-$sweep}
CLEAF_CASE=w3:1000000@500

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

# `variants <case>`: FMM3D's variants for a case.
variants() {
    parts "$1"
    if [ "$kd" = 1000 ]; then
        echo fastker
    elif [ "$density" != random ]; then
        echo default
    else
        echo "$VARIANTS"
    fi
}

# `state`: off Slurm, the load and the free memory (locust is shared).
state() { [ "$slurm" = 1 ] || echo "load $(uptime | sed 's/.*load average: //'), available $(free -g | awk '/Mem:/ { print $7 }') GB"; }

# `note <line>`: a line in runs.log and on stdout.
note() { echo "$*" | tee -a "$out/runs.log"; }

# `selected <launch name>`: whether ONLY admits it.
selected() { [[ "$1" =~ ${ONLY:-.} ]]; }

# `problem <case>`: writes the case's problem file once (sets prob).
problem() {
    parts "$1"
    prob=$out/problems/$name.bin
    if [ ! -f "$prob" ] && [ -z "${DRY_RUN:-}" ]; then
        "$F" problem --workload "$workload" --kd "$kd" --density "$density" --out "$prob" >> \
            "$out/problems.log" || { note "no problem file for $1"; return 1; }
    fi
}

# `repeat`: FMM3D's calls per launch for the problem file prob: 2 below 10⁶ points, else 1. N is
# read with GNU od (Kathleen, locust); macOS's od misreads the 8-byte integer, so DRY_RUN there
# may show the wrong REPEAT.
repeat() {
    local n
    n=$(od -An -t d8 -j 8 -N 8 "$prob" 2> /dev/null | tr -d ' ')
    if [ -n "$n" ] && [ "$n" -lt 1000000 ]; then echo 2; else echo 1; fi
}

# `driver <prefix> <variant> <bind> <eps> <repeat>`: one launch of FMM3D's driver at one eps on
# every core, under maxrss.py, its log <prefix>.log; uses prob and kd.
driver() {
    local pre=$1 variant=$2 bind=$3 e=$4 repeat=$5 start status env available
    local exe=$FMM3D_ROOT/$variant/bin/fmm3d_driver
    env=(OMP_NUM_THREADS="$cores")
    [ "$bind" = bound ] && env+=(OMP_PLACES=cores OMP_PROC_BIND=close)
    selected "$(basename "$pre")-$e" || return 0
    if [ "$variant" = default ] && [ "$kd" = 1000 ] && [ -z "${DRY_RUN:-}" ]; then
        available=0
        [ "$slurm" = 1 ] || available=$(free -g | awk '/Mem:/ { print $7 }')
        if [ "$available" -lt 250 ]; then
            note "$(basename "$pre") $e: skipped, $available GB available (default at kD = 1000)"
            return 0
        fi
    fi
    if [ -n "${DRY_RUN:-}" ]; then
        echo "$(basename "$pre") $e: ${env[*]} $exe helmholtz $prob $pre $kd 0 $e $repeat"
        return 0
    fi
    # A launch that fails must not leave an earlier run's output to be read as its own.
    rm -f "$pre-$e.out" "$pre-$e.full.out"
    {
        echo "# $(date '+%Y-%m-%d %H:%M:%S') $(hostname) job ${SLURM_JOB_ID:-none} $(state)"
        echo "# ${env[*]} $exe helmholtz $prob $pre $kd 0 $e $repeat"
    } >> "$pre.log"
    start=$SECONDS
    env "${env[@]}" timeout "$TIMEOUT" python3 tools/fmm3d/maxrss.py "$pre.log" "token=$e" -- \
        "$exe" helmholtz "$prob" "$pre" "$kd" 0 "$e" "$repeat" >> "$pre.log" 2>&1
    status=$?
    [ "$status" = 0 ] || echo "FAILED: exit $status" >> "$pre.log"
    note "$(basename "$pre") $e: exit $status, $((SECONDS - start)) s, $(date +%T) $(state)"
}

# `errors <out.md> <label=prefix,…> <eps,…>`: FMM3D's outputs against the case's oracle.
errors() {
    [ -z "${DRY_RUN:-}" ] || return 0
    "$F" errors --workload "$workload" --kd "$kd" --density "$density" --eps "$3" --fmm3d "$2" \
        --problem "$prob" --threads "$cores" > "$1" 2>> "${1%.md}.err" ||
        note "errors $(basename "$1"): exit $?"
}

# `mpi <ranks> <threads>`: the launcher of a split.
mpi() {
    local r=$1 t=$2
    if [ "$slurm" = 1 ] && [ $((r % 2)) -eq 0 ]; then
        echo "mpirun --report-bindings --map-by ppr:$((r / 2)):socket:PE=$t --bind-to core -n $r"
    else
        echo "mpirun --report-bindings --map-by slot:PE=$t --bind-to core -n $r"
    fi
}

# `ours <launch name> <split> <ε> [compare options…]`: one launch of ours at the defaults;
# uses prob, workload, kd, density.
ours() {
    local n=$1 split=$2 e=$3 r t start status
    shift 3
    r=${split%x*} t=${split#*x}
    selected "$n" || return 0
    # Word splitting of the launcher is intended.
    # shellcheck disable=SC2046
    set -- $(mpi "$r" "$t") "$F" compare --workload "$workload" --kd "$kd" --density "$density" \
        --settings defaults --eps "$e" --evaluations 2 --problem "$prob" --threads "$t" "$@"
    if [ -n "${DRY_RUN:-}" ]; then
        echo "$n: $*"
        return 0
    fi
    {
        echo "# $(date '+%Y-%m-%d %H:%M:%S') $(hostname) job ${SLURM_JOB_ID:-none} $(state)"
        echo "# $*"
    } > "$out/$n.err"
    start=$SECONDS
    timeout "$TIMEOUT" "$@" < /dev/null > "$out/$n.md" 2>> "$out/$n.err"
    status=$?
    note "$n: exit $status, $((SECONDS - start)) s, $(date +%T) $(state)"
}

# `reversed <words…>`: the words in reverse order.
reversed() { printf '%s\n' "$@" | awk '{ a[NR] = $0 } END { for (i = NR; i > 0; i--) printf "%s%s", a[i], (i > 1 ? " " : "\n") }'; }

stage_grid() {
    local c v labels
    for c in $CASES; do
        problem "$c" || continue
        labels=
        for v in $(variants "$c"); do
            for e in ${GRID//,/ }; do driver "$out/grid/$name-$v" "$v" free "$e" 1; done
            labels="$labels,$v=$out/grid/$name-$v"
        done
        if selected "grid-$name"; then errors "$out/grid-$name.md" "${labels#,}" "$GRID"; fi
    done
}

stage_ours() {
    local c e s
    for c in $CASES; do
        problem "$c" || continue
        for e in $EPS; do
            # shellcheck disable=SC2086 # the splits are words
            for s in $SPLITS; do ours "ours-$name-$e-$s-a" "$s" "$e"; done
            # shellcheck disable=SC2046,SC2086
            for s in $(reversed $SPLITS); do ours "ours-$name-$e-$s-b" "$s" "$e"; done
        done
    done
}

stage_cleaf() {
    local cl e s splits="$SPLITS ${cores}x1"
    problem "$CLEAF_CASE" || return
    for e in $EPS; do
        for cl in 2 4; do
            # shellcheck disable=SC2086
            for s in $splits; do ours "cleaf$cl-$name-$e-$s-a" "$s" "$e" --c-leaf "$cl"; done
        done
        for cl in 4 2; do
            # shellcheck disable=SC2046,SC2086
            for s in $(reversed $splits); do ours "cleaf$cl-$name-$e-$s-b" "$s" "$e" --c-leaf "$cl"; done
        done
    done
}

# `timed <stage> <list> <bindings>`: the `time` and `bind` stages: per case of the list, each
# (variant, binding) at its eps, launch a in the list's order and b reversed; then `errors`
# over every run into <stage>-<case>.md.
timed() {
    local stage=$1 list=$2 bindings=$3 c v eps launch labels tokens runs run bind reps
    [ -f "$list" ] || { echo "ch4-5-runs.sh: $list is missing (ch4-5-tables.py eps)" >&2; return 1; }
    for c in $CASES; do
        parts "$c"
        grep -q "^$name " "$list" || continue
        problem "$c" || continue
        reps=$(repeat)
        labels= tokens= runs=
        while read -r n v eps; do
            [ "$n" = "$name" ] || continue
            for bind in $bindings; do runs="$runs $v:$bind:$eps"; done
        done < "$list"
        for launch in a b; do
            # shellcheck disable=SC2086
            [ "$launch" = b ] && runs=$(reversed $runs)
            for run in $runs; do
                IFS=: read -r v bind eps <<< "$run"
                for e in ${eps//,/ }; do
                    driver "$out/$stage/$name-$v-$bind-$launch" "$v" "$bind" "$e" "$reps"
                done
                labels="$labels,$v-$bind-$launch=$out/$stage/$name-$v-$bind-$launch"
                tokens="$tokens,$eps"
            done
        done
        # The union of the case's eps, each once; `errors` reads every run at every token,
        # so a token a run lacks is a failed row of that run.
        tokens=$(echo "${tokens#,}" | tr , '\n' | awk '!seen[$0]++' | paste -sd, -)
        if selected "$stage-$name"; then errors "$out/$stage-$name.md" "${labels#,}" "$tokens"; fi
    done
}

if [ -z "$FMM3D_ROOT" ] && [ -z "${DRY_RUN:-}" ]; then
    echo "ch4-5-runs.sh: FMM3D_ROOT is not set" >&2
    exit 2
fi
for stage in "$@"; do
    note "=== $(date '+%Y-%m-%d %H:%M:%S') stage $stage on $(hostname), job ${SLURM_JOB_ID:-none}, $cores cores"
    case "$stage" in
    grid) stage_grid ;;
    ours) stage_ours ;;
    cleaf) stage_cleaf ;;
    time) timed time "$out/time-eps.txt" free ;;
    bind) timed bind "$out/bind-eps.txt" "free bound" ;;
    *)
        echo "ch4-5-runs.sh: unknown stage $stage (grid, ours, cleaf, time, bind)" >&2
        exit 2
        ;;
    esac
done
