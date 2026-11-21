#!/usr/bin/env bash
set -euo pipefail

if (( $# > 1 )); then
    printf 'Usage: %s [OUTPUT_DIR]\n' "$0" >&2
    exit 2
fi

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
eqclib=$(cd -- "$here/../../.." && pwd)
out=$(realpath -m -- "${1:-$eqclib/test/out/chimes/shielded/h2}")
gsl=$(realpath -- "${GSL_DIR:-$HOME/lib/gsl}")
chimes_dir=$(realpath -- "${CHIMES_DIR:-$eqclib/../Chimes}")
chimes_local=$(realpath -- "${CHIMES_LOCAL:-$chimes_dir/.local}")
chimes_data=$(realpath -- "${CHIMES_DATA:-$chimes_dir/chimes-data}")
pegase=$(realpath -- "${PEGASE_DIR:-$eqclib/data/PEGASE2}")
build=$(mktemp -d -t eqclib-h2-comparison.XXXXXX)
trap 'rm -rf -- "$build"' EXIT

# Reproduce the original diagnosis, including the earlier H2 prescriptions.
mkdir -- "$build/src" "$build/results"
cp -- "$eqclib/src/Makefile" "$eqclib"/src/*.[ch] "$build/src/"
printf '\nCoolingCIE.o: CPPFLAGS += -DEQC_CIE_TRACKED_ONLY\n' >> "$build/src/Makefile"
printf 'CPPFLAGS += -DEQC_H2_LEGACY_FORMATION -DEQC_H2_LEGACY_SHIELDING\n' >> "$build/src/Makefile"
make -s -j4 -C "$here" OUT="$build" EQCLIB_DIR="$build/src" GSL_DIR="$gsl" \
    CHIMES_DIR="$chimes_dir" CHIMES_LOCAL="$chimes_local" "$build/h2_eqclib" "$build/h2_chimes"

export LD_LIBRARY_PATH="$chimes_local/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
for mode in equilibrium cooling; do
    prefix=$mode
    if [[ "$mode" == equilibrium ]]; then prefix=eq; fi
    "$build/h2_eqclib" "$mode" default "$pegase" "$build/results/${prefix}_eqclib.csv"
    "$build/h2_eqclib" "$mode" form025 "$pegase" "$build/results/${prefix}_eqclib_form025.csv"
    "$build/h2_chimes" "$mode" matched "$chimes_data" "$build/results/${prefix}_chimes_matched.csv"
done

mkdir -p -- "$out"
mv -- "$build/results/"*.csv "$out/"
printf 'Comparison CSVs saved in %s\n' "$out"
