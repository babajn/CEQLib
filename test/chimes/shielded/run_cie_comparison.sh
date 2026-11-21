#!/usr/bin/env bash
set -euo pipefail

if (( $# > 1 )); then
    printf 'Usage: %s [OUTPUT_DIR]\n' "$0" >&2
    exit 2
fi

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
eqclib=$(cd -- "$here/../../.." && pwd)
out=$(realpath -m -- "${1:-$eqclib/test/out/chimes/shielded/cie_supplement}")
baseline=${CIE_BASE_REF:-bc72f906e7d2d0a32565694ab34684231f82d2a3}
gsl=$(realpath -- "${GSL_DIR:-$HOME/lib/gsl}")
chimes_dir=$(realpath -- "${CHIMES_DIR:-$eqclib/../Chimes}")
chimes_local=$(realpath -- "${CHIMES_LOCAL:-$chimes_dir/.local}")
chimes_data=$(realpath -- "${CHIMES_DATA:-$chimes_dir/chimes-data}")
pegase=$(realpath -- "${PEGASE_DIR:-$eqclib/data/PEGASE2}")
cpu=${BENCH_CPU:-$(awk '/^Cpus_allowed_list:/ {split($2, cpus, /[-,]/); print cpus[1]}' /proc/self/status)}
build=$(mktemp -d -t eqclib-cie-comparison.XXXXXX)
trap 'rm -rf -- "$build"' EXIT

# Keep the updated grain physics and thermal solver across all three builds.
mkdir -- "$build/supplemented-src" "$build/results" "$build/repeats"
cp -- "$eqclib/src/Makefile" "$eqclib"/src/*.[ch] "$build/supplemented-src/"
printf '\nCPPFLAGS += -DEQC_H2_LEGACY_FORMATION -DEQC_H2_LEGACY_SHIELDING\n' >> "$build/supplemented-src/Makefile"
cp -R -- "$build/supplemented-src" "$build/before-src"
cp -R -- "$build/supplemented-src" "$build/after-src"
git -C "$eqclib" show "$baseline:src/CoolingCIE.c" > "$build/before-src/CoolingCIE.c"
printf '\nCoolingCIE.o: CPPFLAGS += -DEQC_CIE_TRACKED_ONLY\n' >> "$build/after-src/Makefile"

versions=(before after supplemented)
for version in "${versions[@]}"; do
    gas_z=0.019
    if [[ "$version" == supplemented ]]; then
        gas_z=ASPLUND_2009_SOLAR_METALLICITY
    fi
    make -s -j4 -C "$here" OUT="$build/$version" EQCLIB_DIR="$build/$version-src" GSL_DIR="$gsl" \
        EQCLIB_GAS_Z="$gas_z" \
        "$build/$version/compare_eqclib"
    "$build/$version/compare_eqclib" equilibrium default "$pegase" "$build/results/eq_$version.csv"
    "$build/$version/compare_eqclib" cooling default "$pegase" "$build/results/cooling_$version.csv"
done

make -s -C "$here" OUT="$build/chimes" CHIMES_DIR="$chimes_dir" CHIMES_LOCAL="$chimes_local" \
    "$build/chimes/compare_chimes"
export LD_LIBRARY_PATH="$chimes_local/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
for network in matched full; do
    "$build/chimes/compare_chimes" equilibrium "$network" "$chimes_data" "$build/results/eq_chimes_$network.csv"
    "$build/chimes/compare_chimes" cooling "$network" "$chimes_data" "$build/results/cooling_chimes_$network.csv"
done

# The driver reports a per-state median from five trials. Repeat full passes on
# one core and rotate their order to reduce the effect of CPU frequency drift.
repeats=7
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for (( repeat=1; repeat<=repeats; repeat++ )); do
    offset=$(( (repeat-1) % ${#versions[@]} ))
    order=("${versions[@]:offset}" "${versions[@]:0:offset}")
    for mode in timing courant; do
        for version in "${order[@]}"; do
            taskset -c "$cpu" "$build/$version/compare_eqclib" "$mode" default "$pegase" \
                "$build/repeats/${mode}_${version}_${repeat}.csv"
        done
    done
done

for mode in timing courant; do
    for version in "${versions[@]}"; do
        cp -- "$build/repeats/${mode}_${version}_1.csv" "$build/results/${mode}_${version}.csv"
    done
    for network in matched full; do
        taskset -c "$cpu" "$build/chimes/compare_chimes" "$mode" "$network" "$chimes_data" \
            "$build/results/${mode}_chimes_$network.csv"
    done
done

"${PYTHON:-python3}" - "$build" "$repeats" <<'PY'
import csv
from collections import defaultdict
from pathlib import Path
from statistics import median
import sys

build = Path(sys.argv[1])
repeat_count = int(sys.argv[2])
results = build / "results"
versions = ("before", "after", "supplemented")


def timing_groups(stage, group_field, version, repeat):
    groups = defaultdict(dict)
    path = build / "repeats" / f"{stage}_{version}_{repeat}.csv"
    with path.open(newline="") as file:
        for row in csv.DictReader(file):
            if row["method"] == "eqclib_eq":
                groups[float(row[group_field])][(row["nH"], row["T0"], row["dt_yr"])] = row
    return groups


with (results / "performance.csv").open("w", newline="") as file:
    writer = csv.writer(file, lineterminator="\n")
    writer.writerow(("stage", "version", "repeat", "group", "cells", "comparable_cells", "cpu_us", "median_cpu_us"))
    for stage, group_field in (("timing", "dt_yr"), ("courant", "mass_msun")):
        for repeat in range(1, repeat_count + 1):
            runs = {version: timing_groups(stage, group_field, version, repeat) for version in versions}
            groups = runs["before"].keys()
            if not groups or any(run.keys() != groups for run in runs.values()):
                raise SystemExit(f"Different benchmark groups for {stage}, repeat {repeat}.")
            for group in sorted(groups):
                tables = {version: run[group] for version, run in runs.items()}
                cells = tables["before"].keys()
                if any(table.keys() != cells for table in tables.values()):
                    raise SystemExit(f"Different benchmark states for {stage}, group {group}, repeat {repeat}.")
                common = sorted(key for key in cells if all(table[key]["status"] == "0" for table in tables.values()))
                if not common:
                    raise SystemExit(f"No comparable successes for {stage}, group {group}, repeat {repeat}.")
                for version, table in tables.items():
                    samples = [float(table[key]["cpu_us_median"]) for key in common]
                    writer.writerow((stage, version, repeat, f"{group:g}", len(cells), len(common),
                                     f"{sum(samples):.12g}", f"{median(samples):.12g}"))
PY

mkdir -p -- "$out"
mv -- "$build/results/"*.csv "$out/"
printf 'Comparison CSVs saved in %s (CPU %s, %d paired timing passes)\n' "$out" "$cpu" "$repeats"
