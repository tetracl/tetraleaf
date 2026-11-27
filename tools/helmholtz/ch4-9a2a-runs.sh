#!/bin/bash
# tools/helmholtz/ch4-9a2a-runs.sh: the machine runs of Helmholtz CH4.9a2a, the leaf harmonics and
# the real radial functions per point against a leaf point's reference parts
# (docs/helmholtz/tasks/CH4.9a2a.md, "Machine runs"), from the repository root after
#
#     cargo build --release -p tetraleaf-validate --example helmholtz_leaf
#     tools/helmholtz/ch4-9a2a-runs.sh <out-dir> <stage…>
#
# Each launch runs `helmholtz_leaf` on one thread for one precision and one degree, at every k̂ of
# KHAT (one setting of what a time decides per launch: the degree; review 2026-11-24, B12), pinned
# to one core where `taskset` exists (CORE, default 2), with LANES lanes (default 8).
# `tools/helmholtz/ch4-9a2a-tables.py` makes the report's tables. Stages:
#   smoke   p = 8 and 17 at k̂ = 2, f64, 64 points (a check of the script)
#   parts   PS (default the leaf levels' degrees of CH4.6a2's cases, 8 12 15 17 19 26 32 41 45)
#           × KHAT (default 0.5,2,4,8,16), f64 then f32, real k
#   lanes   p = 17 and 45 at k̂ = 4 with 4, 8 and 16 lanes, f64 and f32 (the lane count)
# TIMEOUT per launch (default 1800 s), ONLY (a regular expression on the launch names), DRY_RUN
# (print the launches only). Each launch writes <out-dir>/<name>.md and .err and a line in
# <out-dir>/runs.log.
set -u
out=$1
shift
mkdir -p "$out"
out=$(cd "$out" && pwd)
L=${CARGO_TARGET_DIR:-target}/release/examples/helmholtz_leaf
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 \
    BLIS_NUM_THREADS=1
TIMEOUT=${TIMEOUT:-1800}
PS=${PS:-8 12 15 17 19 26 32 41 45}
KHAT=${KHAT:-0.5,2,4,8,16}
LANES=${LANES:-8}
pin=()
command -v taskset > /dev/null && pin=(taskset -c "${CORE:-2}")

# `launch <name> <arguments…>`: one launch under the timeout.
launch() {
    local n=$1
    shift
    [[ "$n" =~ ${ONLY:-.} ]] || return 0
    if [ -n "${DRY_RUN:-}" ]; then
        echo "$n: ${pin[*]+${pin[*]}} $L $*"
        return 0
    fi
    local start=$SECONDS
    {
        echo "# $(date '+%Y-%m-%d %H:%M:%S') $(hostname) job ${SLURM_JOB_ID:-none}"
        echo "# load: $(uptime)"
        echo "# ${pin[*]+${pin[*]}} $L $*"
    } > "$out/$n.err"
    timeout "$TIMEOUT" ${pin[@]+"${pin[@]}"} "$L" "$@" < /dev/null > "$out/$n.md" 2>> "$out/$n.err"
    local status=$?
    echo "$n: exit $status, $((SECONDS - start)) s, $(date +%T)" | tee -a "$out/runs.log"
}

for stage in "$@"; do
    case $stage in
    smoke) launch smoke --ps 8,17 --khat 2 --points 64 --lanes "$LANES" ;;
    parts)
        for precision in f64 f32; do
            for p in $PS; do
                launch "parts-$precision-p$p" --ps "$p" --khat "$KHAT" --precision "$precision" \
                    --lanes "$LANES"
            done
        done
        ;;
    lanes)
        for precision in f64 f32; do
            for b in 4 8 16; do
                launch "lanes-$precision-b$b" --ps 17,45 --khat 4 --precision "$precision" --lanes "$b"
            done
        done
        ;;
    *)
        echo "ch4-9a2a-runs.sh: unknown stage $stage" >&2
        exit 2
        ;;
    esac
done
