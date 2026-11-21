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
versions=(legacy formation photodissociation shielding revised revised_b87)
default_dir=h2_prescription
repeats=7
if "$current_only"; then
    versions=(eqclib)
    default_dir=h2_current
    repeats=1
fi
out=$(realpath -m -- "${1:-$eqclib/test/out/chimes/shielded/$default_dir}")
gsl=$(realpath -- "${GSL_DIR:-$HOME/lib/gsl}")
chimes_dir=$(realpath -- "${CHIMES_DIR:-$eqclib/../Chimes}")
chimes_local=$(realpath -- "${CHIMES_LOCAL:-$chimes_dir/.local}")
chimes_data=$(realpath -- "${CHIMES_DATA:-$chimes_dir/chimes-data}")
pegase=$(realpath -- "${PEGASE_DIR:-$eqclib/data/PEGASE2}")
cpu=${BENCH_CPU:-$(awk '/^Cpus_allowed_list:/ {split($2, cpus, /[-,]/); print cpus[1]}' /proc/self/status)}
build=$(mktemp -d -t eqclib-h2-prescription.XXXXXX)
trap 'rm -rf -- "$build"' EXIT

mkdir -- "$build/src" "$build/results"
cp -- "$eqclib/src/Makefile" "$eqclib"/src/*.[ch] "$build/src/"
for version in "${versions[@]}"; do
    source="$build/$version-src"
    cp -R -- "$build/src" "$source"
    flags=()
    case "$version" in
        legacy) flags=(-DEQC_H2_LEGACY_FORMATION -DEQC_H2_LEGACY_SHIELDING) ;;
        formation) flags=(-DEQC_H2_LEGACY_SHIELDING) ;;
        photodissociation) flags=(-DEQC_H2_LEGACY_FORMATION -DEQC_H2_LEGACY_SHIELDING -DEQC_H2_B87_PHOTO) ;;
        shielding) flags=(-DEQC_H2_LEGACY_FORMATION) ;;
        revised_b87) flags=(-DEQC_H2_B87_PHOTO) ;;
    esac
    # Limit CIE to C, O and Si to match the CHIMES network in this test.
    printf '\nCoolingCIE.o: CPPFLAGS += -DEQC_CIE_TRACKED_ONLY\nCPPFLAGS += %s\n' "${flags[*]}" >> "$source/Makefile"
    make -s -j4 -C "$here" OUT="$build/$version" EQCLIB_DIR="$source" GSL_DIR="$gsl" "$build/$version/h2_eqclib"
    "$build/$version/h2_eqclib" native default "$pegase" "$build/results/eq_$version.csv"
    "$build/$version/h2_eqclib" cooling default "$pegase" "$build/results/cooling_$version.csv"
done

# Use EQCLib's 3 km s^-1 turbulent width in the matched CHIMES reference.
make -s -C "$here" OUT="$build/chimes" CHIMES_DIR="$chimes_dir" CHIMES_LOCAL="$chimes_local" \
    H2_CHIMES_CPPFLAGS=-DCMP_CHIMES_DOPPLER_B=3.0 "$build/chimes/h2_chimes"
export LD_LIBRARY_PATH="$chimes_local/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
"$build/chimes/h2_chimes" native matched "$chimes_data" "$build/results/eq_chimes_matched.csv"
"$build/chimes/h2_chimes" cooling matched "$chimes_data" "$build/results/cooling_chimes_matched.csv"

# Rotate the run order so that each EQCLib variant is timed at different points
# in the sequence. All timings are per-thread CPU time, pinned to one core.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for (( repeat=1; repeat<=repeats; repeat++ )); do
    for (( i=0; i<${#versions[@]}; i++ )); do
        version=${versions[(i+repeat-1)%${#versions[@]}]}
        taskset -c "$cpu" "$build/$version/h2_eqclib" courant default "$pegase" "$build/timing_$version.csv"
        if (( repeat == 1 && i == 0 )); then
            header=$(head -n 1 -- "$build/timing_$version.csv")
            printf 'variant,pass,%s\n' "$header" > "$build/results/timing.csv"
        fi
        sed "1d;s/^/$version,$repeat,/" "$build/timing_$version.csv" >> "$build/results/timing.csv"
    done
done
taskset -c "$cpu" "$build/chimes/h2_chimes" courant matched "$chimes_data" "$build/timing_chimes_matched.csv"
sed '1d;s/^/chimes_matched,1,/' "$build/timing_chimes_matched.csv" >> "$build/results/timing.csv"

mkdir -p -- "$out"
mv -- "$build/results/"*.csv "$out/"
printf 'Comparison CSVs saved in %s (CPU %s, %d EQCLib passes)\n' "$out" "$cpu" "$repeats"
