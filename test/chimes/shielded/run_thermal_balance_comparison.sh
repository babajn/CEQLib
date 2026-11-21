#!/usr/bin/env bash
set -euo pipefail

current_only=false
if [[ ${1:-} == --current-only ]]; then
    current_only=true
    shift
fi
if (( $# > 1 )); then
    printf 'Usage: %s [--current-only] [OUTPUT_DIR]\n' "$0" >&2
    exit 2
fi

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
eqclib=$(cd -- "$here/../../.." && pwd)
default_dir=thermal_balance
if "$current_only"; then
    default_dir=thermal_balance_current
fi
out=$(realpath -m -- "${1:-$eqclib/test/out/chimes/shielded/$default_dir}")
gsl=$(realpath -- "${GSL_DIR:-$HOME/lib/gsl}")
pegase=$(realpath -- "${PEGASE_DIR:-$eqclib/data/PEGASE2}")
chimes_dir=$(realpath -- "${CHIMES_DIR:-$eqclib/../Chimes}")
chimes_local=$(realpath -- "${CHIMES_LOCAL:-$chimes_dir/.local}")
chimes_data=$(realpath -- "${CHIMES_DATA:-$chimes_dir/chimes-data}")
build=$(mktemp -d -t eqclib-thermal-balance.XXXXXX)
trap 'rm -rf -- "$build"' EXIT

mkdir -- "$build/current-src" "$build/results"
cp -- "$eqclib/src/Makefile" "$eqclib"/src/*.[ch] "$build/current-src/"
printf '\nCoolingCIE.o: CPPFLAGS += -DEQC_CIE_TRACKED_ONLY\n' >> "$build/current-src/Makefile"
if ! "$current_only"; then
    cp -R -- "$build/current-src" "$build/legacy-src"
    printf '\nCPPFLAGS += -DEQC_H2_LEGACY_FORMATION -DEQC_H2_LEGACY_SHIELDING\n' >> "$build/legacy-src/Makefile"
fi
make -s -j4 -C "$here" OUT="$build/current" EQCLIB_DIR="$build/current-src" GSL_DIR="$gsl" \
    "$build/current/thermal_eqclib"
if ! "$current_only"; then
    make -s -j4 -C "$here" OUT="$build/legacy" EQCLIB_DIR="$build/legacy-src" GSL_DIR="$gsl" \
        "$build/legacy/thermal_eqclib"
    "$build/legacy/thermal_eqclib" default "$pegase" "$build/results/native_legacy.csv"
fi
make -s -C "$here" OUT="$build/chimes" CHIMES_DIR="$chimes_dir" CHIMES_LOCAL="$chimes_local" \
    "$build/chimes/thermal_chimes"
export LD_LIBRARY_PATH="$chimes_local/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

"$build/current/thermal_eqclib" default "$pegase" "$build/results/native_eqclib.csv"
"$build/chimes/thermal_chimes" matched "$chimes_data" "$build/results/native_chimes.csv"

"${PYTHON:-python3}" - "$build" <<'PY'
import csv
import sys
from pathlib import Path

build = Path(sys.argv[1])
results = build / "results"

def load(name):
    with (results / f"native_{name}.csv").open(newline="") as data:
        return {(row["nH"], row["T"]): row for row in csv.DictReader(data) if row["status"] == "0"}

def write_state(output, name, key, row, molecular_fraction):
    ionized = min(1.0, max(0.0, float(row["xHII"])))
    molecular = min(1.0 - ionized, max(0.0, molecular_fraction))
    metals = [min(1.0, max(0.0, float(row[ion]))) for ion in ("xCII", "xOI", "xSiII")]
    print(name, *key, molecular, float(row["xe"]), ionized, *metals, file=output)

eq, chimes = (load(name) for name in ("eqclib", "chimes"))
legacy = load("legacy") if (results / "native_legacy.csv").is_file() else {}
with (build / "compositions.txt").open("w") as output:
    for key in sorted(eq, key=lambda point: tuple(map(float, point))):
        row = eq[key]
        write_state(output, "eq_state", key, row, float(row["fH2"]))
        if key in legacy:
            write_state(output, "eq_legacy_h2", key, row, float(legacy[key]["fH2"]))
        if key in chimes:
            row = chimes[key]
            write_state(output, "chimes_state", key, row, float(row["fH2"]))
PY

"$build/current/thermal_eqclib" default "$pegase" "$build/results/shared_eqclib.csv" "$build/compositions.txt"
"$build/chimes/thermal_chimes" matched "$chimes_data" "$build/results/shared_chimes.csv" "$build/compositions.txt"
mkdir -p -- "$out"
mv -- "$build/results/"*.csv "$out/"
printf 'Thermal-budget CSVs saved in %s\n' "$out"
