#!/bin/sh
# spikes/fft-plans/run.sh: S9's runs on one machine (Helmholtz CH0.9), into results/<machine>/.
#
#     run.sh <machine> [all|tests|record|times|fingerprints] [mpirun options...]
#
# On locust: tools/gh200/remote.sh 'spikes/fft-plans/run.sh locust'. On Kathleen, in jobs:
#     tools/kathleen/submit.sh run.sbatch 'spikes/fft-plans/run.sh kathleen local'
#     tools/kathleen/submit.sh --qos=test --nodes=2 -- run.sbatch \
#         'spikes/fft-plans/run.sh kathleen fingerprints --map-by ppr:1:node'
# `local` is everything but the fingerprints. The times run one back end and grid per launch
# (a second setting in one launch runs slower on Kathleen; review 2026-11-24, B12), each on one
# thread; the BLAS thread variables at 1.
set -eu
cd "$(dirname "$0")"
machine=$1
what=${2:-all}
[ $# -ge 2 ] && shift 2 || shift $#
out=results/$machine
mkdir -p "$out"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
bin=${CARGO_TARGET_DIR:-../../target/spikes-fft-plans}/release/examples

cargo build --release --examples
cargo test --release --no-run

note() {
    { date -u '+%Y-%m-%dT%H:%M:%SZ'; hostname; uptime; } >>"$out/load.txt"
}

case $what in all | local | tests)
    note
    cargo test --release -- --nocapture --test-threads=1 >"$out/tests.txt" 2>&1
    ;;
esac
case $what in all | local | record)
    note
    "$bin/record" >"$out/record.md"
    ;;
esac
case $what in all | local | times)
    note
    for setting in "product product" "fftw product" "fftw smooth" "fftw pow2" \
        "fftw-even smooth" "fftw-split smooth" "rustfft product" "rustfft smooth" \
        "rustfft pow2"; do
        # shellcheck disable=SC2086
        set -- $setting
        "$bin/times" "$1" "$2" >"$out/times-$1-$2.md"
    done
    note
    ;;
esac
case $what in all | fingerprints)
    note
    mpirun -n 2 "$@" "$bin/fingerprints" >"$out/fingerprints.md"
    ;;
esac
echo "S9 runs ($what) on $machine: $out"
