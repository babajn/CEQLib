#CC = mpicc
CC = gcc
EQCLIB_DIR ?= ../src
GSL_DIR ?= /home/babajn/lib/gsl
PYTHON ?= python3

OPTIMIZE += -O3 -std=c99
OPTIMIZE += -finline-functions -ffast-math -funroll-loops -frerun-cse-after-loop -fexpensive-optimizations -fomit-frame-pointer -msse3 -mfpmath=sse
OPTIMIZE += -Winline -Wundef -Wmissing-noreturn -Wstrict-prototypes 
OPTIMIZE += -g 
OPTIMIZE += -DDEBUG_LEVEL=0

INCLUDEPATH += -I$(GSL_DIR)/include -I$(EQCLIB_DIR)

LDFLAGS += -L$(EQCLIB_DIR) -lEQCLib
LDFLAGS += -L$(GSL_DIR)/lib -Wl,-rpath,$(GSL_DIR)/lib -lgsl -lgslcblas

#include makeflags

ifeq ($(COMPILER),"craycc")
OPTIMIZE += -msse3 -m3dnow -m64 -mcpu=opteron -march=opteron 
endif
ifeq ($(COMPILER),"gcc")
OPTIMIZE += -msse3 -m3dnow -m64 -mcpu=opteron -march=opteron 
endif

CFLAGS = $(OPTIMIZE) $(OPTIONS) $(INCLUDEPATH) 

LDFLAGS += -lm

OBJS	=	main.o \
        RandomNumberGenerator.o \

SRCS = main.c \
        RandomNumberGenerator.c \
        
HEADS = Astro.h \
        Constants.h \
        RandomNumberGenerator.h \

ALL:   start

.PHONY: start eqclib check cell-check benchmark hii hii-benchmark hii-plot co co-plot molecular-settings-benchmark molecular-settings-plot wolfire-plot noneq noneq-benchmark noneq-plot
eqclib:
	$(MAKE) -C $(EQCLIB_DIR) GSL_DIR=$(GSL_DIR) libEQCLib.a

start: $(OBJS) eqclib
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) -o test.out

check: EquilibriumTest.o HydrogenChemistryTest.o ThermochemistryModeTest.o eqclib
	mkdir -p out
	$(CC) $(CFLAGS) -fno-fast-math EquilibriumTest.o $(LDFLAGS) -o out/equilibrium_test.out
	$(CC) $(CFLAGS) -fno-fast-math HydrogenChemistryTest.o $(LDFLAGS) -o out/hydrogen_chemistry_test.out
	$(CC) $(CFLAGS) -fno-fast-math ThermochemistryModeTest.o $(LDFLAGS) -o out/thermochemistry_mode_test.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/check.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/equilibrium_test.out" "$(abspath ../data/PEGASE2)"; \
	"$(CURDIR)/out/hydrogen_chemistry_test.out" "$(abspath ../data/PEGASE2)"; \
	"$(CURDIR)/out/thermochemistry_mode_test.out" "$(abspath ../data/PEGASE2)"

EquilibriumTest.o HydrogenChemistryTest.o ThermochemistryModeTest.o: override CFLAGS += -fno-fast-math
HydrogenChemistryTest.o ThermochemistryModeTest.o: $(EQCLIB_DIR)/ChemicalEquilibriumCooling.h

cell-check: eqclib
	$(MAKE) -C $(EQCLIB_DIR) GSL_DIR=$(GSL_DIR) eqclib-cell
	$(PYTHON) CellRunnerTest.py "$(EQCLIB_DIR)/eqclib-cell"

benchmark: ThermalBenchmark.c eqclib
	mkdir -p out
	$(CC) $(CFLAGS) -fno-fast-math -fopenmp ThermalBenchmark.c $(LDFLAGS) -o out/thermal_benchmark.out
	./out/thermal_benchmark.out atomic
	./out/thermal_benchmark.out molecular

hii: HIIRegionTest.c HIIRegionFixture.h eqclib
	mkdir -p out/hii
	$(CC) $(CFLAGS) -fno-fast-math HIIRegionTest.c $(LDFLAGS) -o out/hii_region_test.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/hii.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/hii_region_test.out" "$(abspath ../data/PEGASE2)"; \
	mv equilibrium.csv thermal_curve.csv reactions.csv "$(CURDIR)/out/hii/"

hii-benchmark: HIIRegionBenchmark.c HIIRegionFixture.h eqclib
	mkdir -p out/hii
	$(CC) $(CFLAGS) -fno-fast-math HIIRegionBenchmark.c $(LDFLAGS) -o out/hii_region_benchmark.out
	@set -e; csv=$$(mktemp "$(CURDIR)/out/hii/benchmark.XXXXXX"); \
	trap 'rm -f "$$csv"' EXIT; \
	./out/hii_region_benchmark.out "$(abspath ../data/PEGASE2)" > "$$csv"; \
	mv "$$csv" out/hii/hii_benchmark.csv
	@echo "H II region timings saved in out/hii/hii_benchmark.csv"

hii-plot: hii
	$(PYTHON) ../plot/plot_hii_region.py

co: MolecularTest.c MolecularFixture.h eqclib
	mkdir -p out/co
	$(CC) $(CFLAGS) -fno-fast-math MolecularTest.c $(LDFLAGS) -o out/molecular_test.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/co.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/molecular_test.out" "$(abspath ../data/PEGASE2)"; \
	mv equilibrium.csv column_scan.csv "$(CURDIR)/out/co/"

co-plot: co
	$(PYTHON) ../plot/plot_co_chemistry.py

molecular-settings-benchmark: MolecularSettingsBenchmark.c MolecularFixture.h eqclib
	mkdir -p out/co
	$(CC) $(CFLAGS) -fno-fast-math MolecularSettingsBenchmark.c $(LDFLAGS) -o out/molecular_settings_benchmark.out
	@set -e; csv=$$(mktemp "$(CURDIR)/out/co/settings.XXXXXX"); \
	trap 'rm -f "$$csv"' EXIT; \
	./out/molecular_settings_benchmark.out "$(abspath ../data/PEGASE2)" > "$$csv"; \
	mv "$$csv" out/co/molecular_settings.csv

molecular-settings-plot: molecular-settings-benchmark
	$(PYTHON) ../plot/plot_molecular_settings.py

wolfire-plot:
	$(PYTHON) ../plot/plot_wolfire1995.py

noneq: NonEquilibriumTest.c eqclib
	mkdir -p out/noneq
	$(CC) $(CFLAGS) -fno-fast-math NonEquilibriumTest.c $(LDFLAGS) -o out/non_equilibrium_test.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/noneq.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/non_equilibrium_test.out" "$(abspath ../data/PEGASE2)"; \
	mv evolution.csv evolution_direct.csv recombination.csv recombination_direct.csv \
	   convergence.csv convergence_direct.csv "$(CURDIR)/out/noneq/"

noneq-benchmark: NonEquilibriumBenchmark.c eqclib
	mkdir -p out/noneq
	$(CC) $(CFLAGS) -fno-fast-math NonEquilibriumBenchmark.c $(LDFLAGS) -o out/non_equilibrium_benchmark.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/noneq-bench.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/non_equilibrium_benchmark.out" "$(abspath ../data/PEGASE2)"; \
	mv benchmark.csv "$(CURDIR)/out/noneq/"

noneq-plot: noneq noneq-benchmark
	$(PYTHON) ../plot/plot_noneq.py

clean:
	rm -f $(OBJS) test.out EquilibriumTest.o HydrogenChemistryTest.o ThermochemistryModeTest.o out/equilibrium_test.out out/hydrogen_chemistry_test.out out/thermochemistry_mode_test.out out/thermal_benchmark.out out/hii_region_test.out out/hii_region_benchmark.out out/molecular_test.out out/molecular_settings_benchmark.out out/non_equilibrium_test.out out/non_equilibrium_benchmark.out

tag:
	ctags *.[c,h]

pkg: 
	tar cvf plot_anime.tar ./*.[ch] ./*.rb Makefile tags
	gzip plot_anime.tar


main.o: main.c Makefile $(HEADS) 
