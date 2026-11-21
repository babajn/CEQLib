#!/usr/bin/env bash
set -euo pipefail

if (( $# > 1 )); then
    printf 'Usage: %s [OUTPUT_DIR]\n' "$0" >&2
    exit 2
fi

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
eqclib=$(cd -- "$here/../../.." && pwd)
out=$(realpath -m -- "${1:-$eqclib/test/out/chimes/shielded/dust}")
baseline=${DUST_BASE_REF:-bc72f906e7d2d0a32565694ab34684231f82d2a3}
gsl=$(realpath -- "${GSL_DIR:-$HOME/lib/gsl}")
chimes_dir=$(realpath -- "${CHIMES_DIR:-$eqclib/../Chimes}")
chimes_local=$(realpath -- "${CHIMES_LOCAL:-$chimes_dir/.local}")
chimes_data=$(realpath -- "${CHIMES_DATA:-$chimes_dir/chimes-data}")
pegase=$(realpath -- "${PEGASE_DIR:-$eqclib/data/PEGASE2}")
build=$(mktemp -d -t eqclib-dust-comparison.XXXXXX)
trap 'rm -rf -- "$build"' EXIT

# Vary the grain implementation with the same build flags and historical CIE table.
mkdir -- "$build/after-src" "$build/results"
cp -- "$eqclib/src/Makefile" "$eqclib"/src/*.[ch] "$build/after-src/"
printf '\nCPPFLAGS += -DEQC_H2_LEGACY_FORMATION -DEQC_H2_LEGACY_SHIELDING\n' >> "$build/after-src/Makefile"
git -C "$eqclib" show "$baseline:src/CoolingCIE.c" > "$build/after-src/CoolingCIE.c"
cp -R -- "$build/after-src" "$build/before-src"
git -C "$eqclib" show "$baseline:src/GrainSurfaceChemistry.c" > "$build/before-src/GrainSurfaceChemistry.c"

for version in before after; do
    make -s -j4 -C "$here" OUT="$build/$version" EQCLIB_DIR="$build/$version-src" GSL_DIR="$gsl" EQCLIB_GAS_Z=0.019 \
        "$build/$version/compare_eqclib"
    "$build/$version/compare_eqclib" equilibrium default "$pegase" "$build/results/eq_$version.csv"
    "$build/$version/compare_eqclib" cooling default "$pegase" "$build/results/cooling_$version.csv"
done

"${PYTHON:-python3}" - "$build/results" <<'PY'
import csv
from pathlib import Path
import sys

directory = Path(sys.argv[1])

def read_low_temperature_grid(version):
    with (directory / f"eq_{version}.csv").open(newline="") as file:
        return {(row["nH"], row["T0"]): row for row in csv.DictReader(file) if float(row["T0"]) <= 1e4}

def read_low_start_histories(version):
    with (directory / f"cooling_{version}.csv").open(newline="") as file:
        return [{key: value for key, value in row.items() if key != "cpu_us"}
                for row in csv.DictReader(file) if float(row["T0"]) <= 1e4]

before = read_low_temperature_grid("before")
after = read_low_temperature_grid("after")
if not before or before.keys() != after.keys():
    raise SystemExit("Low-temperature grids differ.")
for cell, old in before.items():
    new = after[cell]
    if old["status"] != "0" or new["status"] != "0" or any(old[key] != new[key] for key in ("cool", "heat")):
        raise SystemExit(f"Low-temperature cooling or heating differs at nH={cell[0]}, T={cell[1]}.")
print(f"Identical cooling and heating at all {len(before)} fixed-temperature cells up to 10,000 K.")

cold_before = read_low_start_histories("before")
cold_after = read_low_start_histories("after")
if not cold_before or cold_before != cold_after:
    raise SystemExit("Histories starting at or below 10,000 K differ.")
print(f"Identical low-start histories ({len(cold_before)} rows excluding CPU times), "
      "including prefixes through any failed update.")
PY

make -s -C "$here" OUT="$build/chimes" CHIMES_DIR="$chimes_dir" CHIMES_LOCAL="$chimes_local" \
    "$build/chimes/compare_chimes"
export LD_LIBRARY_PATH="$chimes_local/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
"$build/chimes/compare_chimes" equilibrium full "$chimes_data" "$build/results/eq_chimes_full.csv"
"$build/chimes/compare_chimes" cooling full "$chimes_data" "$build/results/cooling_chimes_full.csv"

mkdir -p -- "$out"
mv -- "$build/results/"*.csv "$out/"
printf 'Comparison CSVs saved in %s\n' "$out"
