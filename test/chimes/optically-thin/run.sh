#!/usr/bin/env bash
set -euo pipefail

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
eqclib=$(cd -- "$here/../../.." && pwd)
chimes_dir=$(realpath -- "${CHIMES_DIR:-$eqclib/../Chimes}")
out=$(realpath -m -- "${1:-$eqclib/test/out/chimes/optically-thin}")
gsl=$(realpath -- "${GSL_DIR:-$HOME/lib/gsl}")
chimes_local=$(realpath -- "${CHIMES_LOCAL:-$chimes_dir/.local}")
build=$(mktemp -d)
trap 'rm -rf -- "$build"' EXIT

mkdir -p -- "$out" "$build/eqclib" "$build/chimes/src"
cp -- "$eqclib/src/Makefile" "$eqclib"/src/*.[ch] "$build/eqclib/"
chimes_makefile="$chimes_dir/chimes/Makefile"
if [[ ! -f "$chimes_makefile" ]]; then
    chimes_makefile="$chimes_dir/chimes/Makefile_template"
fi
cp -- "$chimes_makefile" "$build/chimes/Makefile"
cp -- "$chimes_dir"/chimes/src/*.[ch] "$build/chimes/src/"
# Match the gas composition used by CHIMES in the mixed cases.
make -s -C "$build/eqclib" -j4 GSL_DIR="$gsl" CPPFLAGS="-I$gsl/include -DEQC_CIE_TRACKED_ONLY"
make -s -C "$build/chimes" -j4 EXEC=libchimes.so CFLAGS='-O3 -DCHIMES_USE_DOUBLE_PRECISION' \
    SUNDIALS_PATH="$chimes_local" HDF5_INCL_PATH="$chimes_local/include" HDF5_LIB_PATH="$chimes_local/lib"

gcc -std=c99 -O3 -Wall -Wextra -DCHIMES_USE_DOUBLE_PRECISION \
    -I"$build/eqclib" -I"$gsl/include" -I"$build/chimes/src" -I"$chimes_local/include" \
    "$here/compare.c" "$build/eqclib/libEQCLib.a" "$build/chimes/libchimes.so" \
    -L"$gsl/lib" -Wl,-rpath,"$gsl/lib" -Wl,-rpath-link,"$chimes_local/lib" \
    -lgsl -lgslcblas -lm -o "$build/compare"

export LD_LIBRARY_PATH="$chimes_local/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cpu=${BENCH_CPU:-$(awk '/^Cpus_allowed_list:/ {split($2, cpus, /[-,]/); print cpus[1]}' /proc/self/status)}
for group in hydrogen mixed; do
    if ! taskset -c "$cpu" "$build/compare" "$eqclib/data/PEGASE2" \
        "$chimes_dir/chimes-data/chimes_main_data.hdf5" "$build" "$group" 2>"$build/$group.log" \
        || [[ -s "$build/$group.log" ]]; then
        head -30 "$build/$group.log" >&2
        exit 1
    fi
done

for table in history steps timings; do
    {
        cat -- "$build/hydrogen-$table.csv"
        tail -n +2 -- "$build/mixed-$table.csv"
    } > "$out/$table.csv"
done
"${PYTHON:-python3}" "$here/summarize.py" "$out"
