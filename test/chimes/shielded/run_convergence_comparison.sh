#!/usr/bin/env bash
set -euo pipefail

if (( $# > 1 )); then
    printf 'Usage: %s [OUTPUT_DIR]\n' "$0" >&2
    exit 2
fi

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
eqclib=$(cd -- "$here/../../.." && pwd)
out=$(realpath -m -- "${1:-$eqclib/test/out/chimes/shielded/convergence}")
baseline=${CONVERGENCE_BASE_REF:-bc72f906e7d2d0a32565694ab34684231f82d2a3}
gsl=$(realpath -- "${GSL_DIR:-$HOME/lib/gsl}")
pegase=$(realpath -- "${PEGASE_DIR:-$eqclib/data/PEGASE2}")
cpu=${BENCH_CPU:-$(awk '/^Cpus_allowed_list:/ {split($2, cpus, /[-,]/); print cpus[1]}' /proc/self/status)}
build=$(mktemp -d -t eqclib-convergence.XXXXXX)
trap 'rm -rf -- "$build"' EXIT

# Hold the dust model and all build flags fixed across the two executables.
mkdir -- "$build/after-src" "$build/results" "$build/repeats"
cp -- "$eqclib/src/Makefile" "$eqclib"/src/*.[ch] "$build/after-src/"
printf '\nCPPFLAGS += -DEQC_H2_LEGACY_FORMATION -DEQC_H2_LEGACY_SHIELDING\n' >> "$build/after-src/Makefile"
# Reproduce the original two solver stops with the historical CIE table.
git -C "$eqclib" show "$baseline:src/CoolingCIE.c" > "$build/after-src/CoolingCIE.c"
cp -R -- "$build/after-src" "$build/before-src"
git -C "$eqclib" show "$baseline:src/ThermalEquilibrium.c" > "$build/before-src/ThermalEquilibrium.c"

for version in before after; do
    make -s -j4 -C "$here" OUT="$build/$version" EQCLIB_DIR="$build/$version-src" GSL_DIR="$gsl" EQCLIB_GAS_Z=0.019 \
        "$build/$version/compare_eqclib"
    "$build/$version/compare_eqclib" equilibrium default "$pegase" "$build/results/eq_$version.csv"
    "$build/$version/compare_eqclib" cooling default "$pegase" "$build/results/cooling_$version.csv"
done

# Each driver measures five copies of each starting state and reports their median.
# Repeat whole passes on one core and alternate the order to reduce frequency drift.
repeats=7
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for (( repeat=1; repeat<=repeats; repeat++ )); do
    order=(before after)
    if (( repeat % 2 == 0 )); then
        order=(after before)
    fi
    for mode in timing courant; do
        for version in "${order[@]}"; do
            taskset -c "$cpu" "$build/$version/compare_eqclib" "$mode" default "$pegase" \
                "$build/repeats/${mode}_${version}_${repeat}.csv"
        done
    done
done

for mode in timing courant; do
    for version in before after; do
        cp -- "$build/repeats/${mode}_${version}_1.csv" "$build/results/${mode}_${version}.csv"
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

def read_csv(path):
    with path.open(newline="") as file:
        return list(csv.DictReader(file))

with (results / "performance.csv").open("w", newline="") as file:
    writer = csv.writer(file, lineterminator="\n")
    writer.writerow(("stage", "version", "repeat", "group", "cells", "comparable_cells", "cpu_us", "median_cpu_us"))
    for stage, group_field in (("timing", "dt_yr"), ("courant", "mass_msun")):
        for repeat in range(1, repeat_count + 1):
            by_version = {}
            for version in ("before", "after"):
                groups = defaultdict(dict)
                rows = read_csv(build / "repeats" / f"{stage}_{version}_{repeat}.csv")
                for row in rows:
                    if row["method"] == "eqclib_eq":
                        groups[float(row[group_field])][(row["nH"], row["T0"])] = row
                by_version[version] = groups
            before, after = by_version["before"], by_version["after"]
            if not before or before.keys() != after.keys():
                raise SystemExit(f"Different benchmark groups for {stage}, repeat {repeat}.")
            for group, cells in sorted(before.items()):
                other = after[group]
                if cells.keys() != other.keys():
                    raise SystemExit(f"Different benchmark states for {stage}, group {group}, repeat {repeat}.")
                common = sorted(key for key in cells if cells[key]["status"] == other[key]["status"] == "0")
                if not common:
                    raise SystemExit(f"No comparable successes for {stage}, group {group}, repeat {repeat}.")
                for version, table in (("before", cells), ("after", other)):
                    samples = [float(table[key]["cpu_us_median"]) for key in common]
                    writer.writerow((stage, version, repeat, f"{group:g}", len(cells), len(samples),
                                     f"{sum(samples):.12g}", f"{median(samples):.12g}"))

for version in ("before", "after"):
    histories = defaultdict(list)
    for row in read_csv(results / f"cooling_{version}.csv"):
        if row["method"] == "eqclib_eq":
            histories[(row["nH"], row["T0"])].append(row)
    stopped = {key for key, rows in histories.items() if any(row["status"] != "0" for row in rows)}
    if version == "before" and stopped != {("1", "10000"), ("1", "1e+07")}:
        raise SystemExit(f"Unexpected baseline stops: {sorted(stopped)}.")
    if version == "after" and (len(histories) != 6 or stopped or
                               any(float(rows[-1]["t_yr"]) != 1e9 for rows in histories.values())):
        raise SystemExit(f"The updated equilibrium histories did not all reach 1 Gyr: {sorted(stopped)}.")
    print(f"{version}: {len(histories) - len(stopped)}/{len(histories)} completed equilibrium histories.")
PY

mkdir -p -- "$out"
mv -- "$build/results/"*.csv "$out/"
printf 'Comparison CSVs saved in %s (CPU %s, %d paired timing passes)\n' "$out" "$cpu" "$repeats"
