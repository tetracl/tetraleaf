#!/bin/bash
# tools/helmholtz/ch4-6a2-runs.sh: the machine runs of Helmholtz CH4.6a2, the φ step's time per
# level, the product against rustfft (docs/helmholtz/tasks/CH4.6a2.md, "Machine runs"), from the
# repository root after
#
#     cargo build --release -p tetraleaf-validate --features rustfft --example helmholtz_accuracy
#     tools/helmholtz/ch4-6a2-runs.sh <out-dir> <stage…>
#
# Every `helmholtz_accuracy` launch runs at one rank with every core, `--levels`, §9.1's complex
# charges and samples, gradients, two evaluations (the rows the second's), under the library's
# defaults named (`--m2l auto-pw --diagonal auto --diagonal-build interpolated`, the V and W/X
# rules of the cost model), and one φ step per launch (review 2026-11-24, B12: on Kathleen a
# second setting in one launch runs 1.09–1.20× slower). `tools/helmholtz/ch4-6a2-tables.py`
# makes the report's tables. Stages:
#   smoke   W1 at kD = 30, 10⁻⁶, product then rustfft (a check of the script and the tables)
#   levels  W1 at kD = 250, 500, 1000 at 10⁻⁶, W1 at kD = 500 at 10⁻³ and 10⁻⁸, W3 at N = 10⁶,
#           kD = 500 at 10⁻⁶: per case product, rustfft, rustfft, product (ABBA), then
#           `--diagonal-costs form` with the product and with rustfft (every level the bound
#           reaches on the form, for the time rule's question)
#   waves   the cases where CH4.3b2c found the waves kept against the form (W1 at kD = 100, 200
#           at 10⁻⁷ and 10⁻⁸; W3 at N = 10⁶, kD = 200 at 10⁻⁶ and 10⁻⁸ and kD = 500 at 10⁻⁸): per
#           case the defaults with the product and with rustfft, then `--diagonal-costs form`
#           with each (the time rule's question on levels the waves keep)
#   ffts    spike S9's `times` (spikes/fft-plans, its own workspace and target directory), one
#           thread: `rustfft product`, `rustfft smooth` and `product product` at the degrees PS
#           (comma-separated, required: the diagonal levels' degrees of `levels`)
# TIMEOUT per launch (default 3600 s), ONLY (a regular expression on the launch names),
# DRY_RUN (print the launches only). Each launch writes <out-dir>/<name>.md and .err and a line
# in <out-dir>/runs.log.
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
COMMON=(--near-v cost --near-wx cost --levels --m2l auto-pw --evaluations 2
    --diagonal auto --diagonal-build interpolated)

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

# `launch <name> <arguments…>`: one launch of the example at one rank with every core.
launch() {
    local n=$1
    shift
    # shellcheck disable=SC2046
    record "$n" $(mpi) "$A" "${COMMON[@]}" --threads "$cores" "$@"
}

# `case_runs <name> <arguments…>`: one case's six launches (the stage `levels`).
case_runs() {
    local n=$1
    shift
    launch "$n-product-1" --phi-step product "$@"
    launch "$n-rustfft-1" --phi-step rustfft "$@"
    launch "$n-rustfft-2" --phi-step rustfft "$@"
    launch "$n-product-2" --phi-step product "$@"
    launch "$n-form-product" --phi-step product --diagonal-costs form "$@"
    launch "$n-form-rustfft" --phi-step rustfft --diagonal-costs form "$@"
}

for stage in "$@"; do
    case $stage in
    smoke)
        launch smoke-product --phi-step product --workload w1 --kd 30 --eps 1e-6
        launch smoke-rustfft --phi-step rustfft --workload w1 --kd 30 --eps 1e-6
        ;;
    levels)
        for kd in 250 500 1000; do
            case_runs "w1-$kd-1e-6" --workload w1 --kd "$kd" --eps 1e-6
        done
        for eps in 1e-3 1e-8; do
            case_runs "w1-500-$eps" --workload w1 --kd 500 --eps "$eps"
        done
        case_runs "w3-500-1e-6" --workload w3:1000000 --kd 500 --eps 1e-6
        ;;
    waves)
        for c in w1:100:1e-7 w1:100:1e-8 w1:200:1e-7 w1:200:1e-8 w3:1000000:200:1e-6 \
            w3:1000000:200:1e-8 w3:1000000:500:1e-8; do
            eps=${c##*:}
            rest=${c%:*}
            kd=${rest##*:}
            w=${rest%:*}
            n="waves-${w%%:*}-$kd-$eps"
            for s in "product-1 product" "rustfft-1 rustfft" "form-product product" \
                "form-rustfft rustfft"; do
                # shellcheck disable=SC2086
                set -- $s
                costs=()
                [[ $1 == form-* ]] && costs=(--diagonal-costs form)
                launch "$n-$1" --phi-step "$2" ${costs[@]+"${costs[@]}"} --workload "$w" --kd "$kd" \
                    --eps "$eps"
            done
        done
        ;;
    ffts)
        if [ -z "${PS:-}" ]; then
            echo "ch4-6a2-runs.sh: ffts needs PS, the degrees" >&2
            exit 2
        fi
        spike=$(cd spikes/fft-plans && pwd)
        bin=${SPIKE_TARGET:-$(pwd)/target/spikes-fft-plans}
        if [ -z "${DRY_RUN:-}" ]; then
            (cd "$spike" && CARGO_TARGET_DIR=$bin cargo build --release --example times) ||
                exit 1
        fi
        for setting in "rustfft product" "rustfft smooth" "product product"; do
            # shellcheck disable=SC2086
            set -- $setting
            record "ffts-$1-$2" "$bin/release/examples/times" "$1" "$2" "$PS"
        done
        ;;
    *)
        echo "ch4-6a2-runs.sh: unknown stage $stage" >&2
        exit 2
        ;;
    esac
done
