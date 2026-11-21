#CC = mpicc
CC = gcc
EQCLIB_DIR ?= ../src
GSL_DIR ?= $(or $(shell gsl-config --prefix 2>/dev/null),$(wildcard $(HOME)/lib/gsl),/usr)
export PEGASE_DIR ?= $(abspath ../data/PEGASE2)
PYTHON ?= python3

OPTIMIZE += -O3 -std=c99
OPTIMIZE += -finline-functions -ffast-math -funroll-loops -frerun-cse-after-loop -fexpensive-optimizations -fomit-frame-pointer
OPTIMIZE += -Winline -Wundef -Wmissing-noreturn -Wstrict-prototypes 
OPTIMIZE += -g 
OPTIMIZE += -DDEBUG_LEVEL=0

INCLUDEPATH += -Icommon -I$(GSL_DIR)/include -I$(EQCLIB_DIR)

LDFLAGS += -L$(EQCLIB_DIR) -lEQCLib
LDFLAGS += -L$(GSL_DIR)/lib -Wl,-rpath,$(GSL_DIR)/lib -lgsl -lgslcblas

#include makeflags

ifneq ($(filter x86_64% i386% i486% i586% i686%,$(shell $(CC) -dumpmachine)),)
OPTIMIZE += -msse3 -mfpmath=sse
ifeq ($(COMPILER),"craycc")
OPTIMIZE += -msse3 -m3dnow -m64 -mcpu=opteron -march=opteron 
endif
ifeq ($(COMPILER),"gcc")
OPTIMIZE += -msse3 -m3dnow -m64 -mcpu=opteron -march=opteron 
endif
endif

CFLAGS = $(OPTIMIZE) $(OPTIONS) $(INCLUDEPATH) 

LDFLAGS += -lm

WOLFIRE_OBJS = out/wolfire1995_main.o out/RandomNumberGenerator.o
UNIT_OBJS = out/EquilibriumTest.o out/HydrogenChemistryTest.o out/ThermochemistryModeTest.o
COMMON_HEADERS = common/Astro.h common/Constants.h
TEST_SOURCES = $(wildcard benchmark/*.[ch] co/*.[ch] common/*.[ch] hii/*.[ch] noneq/*.[ch] unit/*.[ch] wolfire1995/*.[ch])

ALL:   start

.PHONY: start eqclib check cell-check notebook-data benchmark hii hii-benchmark hii-plot co co-plot molecular-settings-benchmark molecular-settings-plot wolfire-plot noneq noneq-benchmark noneq-hot noneq-plot clean tag pkg
eqclib:
	$(MAKE) -C $(EQCLIB_DIR) GSL_DIR=$(GSL_DIR) libEQCLib.a

out:
	mkdir -p $@

start: $(WOLFIRE_OBJS) eqclib
	$(CC) $(CFLAGS) $(WOLFIRE_OBJS) $(LDFLAGS) -o test.out

out/wolfire1995_main.o: wolfire1995/main.c $(COMMON_HEADERS) Makefile | out
	$(CC) $(CFLAGS) -c $< -o $@

out/RandomNumberGenerator.o: wolfire1995/RandomNumberGenerator.c wolfire1995/RandomNumberGenerator.h Makefile | out
	$(CC) $(CFLAGS) -c $< -o $@

check: $(UNIT_OBJS) eqclib
	$(CC) $(CFLAGS) -fno-fast-math out/EquilibriumTest.o $(LDFLAGS) -o out/equilibrium_test.out
	$(CC) $(CFLAGS) -fno-fast-math out/HydrogenChemistryTest.o $(LDFLAGS) -o out/hydrogen_chemistry_test.out
	$(CC) $(CFLAGS) -fno-fast-math out/ThermochemistryModeTest.o $(LDFLAGS) -o out/thermochemistry_mode_test.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/check.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/equilibrium_test.out" "$$PEGASE_DIR"; \
	"$(CURDIR)/out/hydrogen_chemistry_test.out" "$$PEGASE_DIR"; \
	"$(CURDIR)/out/thermochemistry_mode_test.out" "$$PEGASE_DIR"

$(UNIT_OBJS): out/%.o: unit/%.c $(COMMON_HEADERS) $(EQCLIB_DIR)/ChemicalEquilibriumCooling.h Makefile | out
	$(CC) $(CFLAGS) -fno-fast-math -c $< -o $@

out/HydrogenChemistryTest.o: $(EQCLIB_DIR)/CoolingCIE.h $(EQCLIB_DIR)/CoolingHeatingRates.h

cell-check: eqclib
	$(MAKE) -C $(EQCLIB_DIR) GSL_DIR=$(GSL_DIR) eqclib-cell
	$(PYTHON) unit/CellRunnerTest.py "$(EQCLIB_DIR)/eqclib-cell"

# Keep timing runs sequential even when the outer make uses -j.
notebook-data:
	$(MAKE) -j1 check cell-check hii co noneq hii-benchmark molecular-settings-benchmark noneq-benchmark noneq-hot
	$(MAKE) -j1 -C wolfire2003 EQCLIB_DIR="$(abspath $(EQCLIB_DIR))" GSL_DIR="$(GSL_DIR)" run

benchmark: benchmark/ThermalBenchmark.c eqclib
	mkdir -p out
	$(CC) $(CFLAGS) -fno-fast-math -fopenmp $< $(LDFLAGS) -o out/thermal_benchmark.out
	./out/thermal_benchmark.out atomic
	./out/thermal_benchmark.out molecular

hii: hii/HIIRegionTest.c hii/HIIRegionFixture.h eqclib
	mkdir -p out/hii
	$(CC) $(CFLAGS) -fno-fast-math $< $(LDFLAGS) -o out/hii_region_test.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/hii.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/hii_region_test.out" "$$PEGASE_DIR"; \
	mv equilibrium.csv thermal_curve.csv reactions.csv "$(CURDIR)/out/hii/"

hii-benchmark: hii/HIIRegionBenchmark.c hii/HIIRegionFixture.h eqclib
	mkdir -p out/hii
	$(CC) $(CFLAGS) -fno-fast-math $< $(LDFLAGS) -o out/hii_region_benchmark.out
	@set -e; csv=$$(mktemp "$(CURDIR)/out/hii/benchmark.XXXXXX"); \
	trap 'rm -f "$$csv"' EXIT; \
	./out/hii_region_benchmark.out "$$PEGASE_DIR" > "$$csv"; \
	mv "$$csv" out/hii/hii_benchmark.csv
	@echo "H II region timings saved in out/hii/hii_benchmark.csv"

hii-plot: hii
	$(PYTHON) ../plot/plot_hii_region.py

co: co/MolecularTest.c co/MolecularFixture.h eqclib
	mkdir -p out/co
	$(CC) $(CFLAGS) -fno-fast-math $< $(LDFLAGS) -o out/molecular_test.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/co.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/molecular_test.out" "$$PEGASE_DIR"; \
	mv equilibrium.csv column_scan.csv "$(CURDIR)/out/co/"

co-plot: co
	$(PYTHON) ../plot/plot_co_chemistry.py

molecular-settings-benchmark: co/MolecularSettingsBenchmark.c co/MolecularFixture.h eqclib
	mkdir -p out/co
	$(CC) $(CFLAGS) -fno-fast-math $< $(LDFLAGS) -o out/molecular_settings_benchmark.out
	@set -e; csv=$$(mktemp "$(CURDIR)/out/co/settings.XXXXXX"); \
	trap 'rm -f "$$csv"' EXIT; \
	./out/molecular_settings_benchmark.out "$$PEGASE_DIR" > "$$csv"; \
	mv "$$csv" out/co/molecular_settings.csv

molecular-settings-plot: molecular-settings-benchmark
	$(PYTHON) ../plot/plot_molecular_settings.py

wolfire-plot:
	$(PYTHON) ../plot/plot_wolfire1995.py

noneq: noneq/NonEquilibriumTest.c eqclib
	mkdir -p out/noneq
	$(CC) $(CFLAGS) -fno-fast-math $< $(LDFLAGS) -o out/non_equilibrium_test.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/noneq.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/non_equilibrium_test.out" "$$PEGASE_DIR"; \
	mv evolution.csv evolution_direct.csv recombination.csv recombination_direct.csv \
	   convergence.csv convergence_direct.csv "$(CURDIR)/out/noneq/"

noneq-benchmark: noneq/NonEquilibriumBenchmark.c eqclib
	mkdir -p out/noneq
	$(CC) $(CFLAGS) -fno-fast-math $< $(LDFLAGS) -o out/non_equilibrium_benchmark.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/noneq-bench.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/non_equilibrium_benchmark.out" "$$PEGASE_DIR"; \
	mv benchmark.csv "$(CURDIR)/out/noneq/"

# Use the production archive for all variants. Override just its cooling
# module in the reference executables.
noneq-hot: noneq/NonEquilibriumHotCoolingTest.c eqclib
	mkdir -p out/noneq
	$(CC) $(CFLAGS) -fno-fast-math -c $< -o out/noneq_hot_driver.o
	$(MAKE) -C $(EQCLIB_DIR) GSL_DIR=$(GSL_DIR) ChemicalEquilibriumCoolingLegacy.o ChemicalEquilibriumCoolingLinearCIE.o
	$(CC) out/noneq_hot_driver.o $(EQCLIB_DIR)/libEQCLib.a $(filter-out -lEQCLib,$(LDFLAGS)) -o out/noneq_hot_current.out
	$(CC) out/noneq_hot_driver.o $(EQCLIB_DIR)/ChemicalEquilibriumCoolingLegacy.o $(EQCLIB_DIR)/libEQCLib.a $(filter-out -lEQCLib,$(LDFLAGS)) -o out/noneq_hot_legacy.out
	$(CC) out/noneq_hot_driver.o $(EQCLIB_DIR)/ChemicalEquilibriumCoolingLinearCIE.o $(EQCLIB_DIR)/libEQCLib.a $(filter-out -lEQCLib,$(LDFLAGS)) -o out/noneq_hot_linear.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/noneq-hot.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	mkdir legacy cie_linear cie_bridge; \
	(cd legacy && "$(CURDIR)/out/noneq_hot_legacy.out" "$$PEGASE_DIR" legacy); \
	(cd cie_linear && "$(CURDIR)/out/noneq_hot_linear.out" "$$PEGASE_DIR" cie_linear); \
	(cd cie_bridge && "$(CURDIR)/out/noneq_hot_current.out" "$$PEGASE_DIR" cie_bridge); \
	for name in hot_cooling hot_benchmark oxygen_after_switch; do \
	    cat "legacy/$$name.csv" > "$$name.csv"; \
	    tail -n +2 "cie_linear/$$name.csv" >> "$$name.csv"; \
	    tail -n +2 "cie_bridge/$$name.csv" >> "$$name.csv"; \
	done; \
	mv hot_cooling.csv hot_benchmark.csv oxygen_after_switch.csv "$(CURDIR)/out/noneq/"

noneq-plot: noneq noneq-benchmark noneq-hot
	$(PYTHON) ../plot/plot_noneq.py

clean:
	rm -f $(WOLFIRE_OBJS) $(UNIT_OBJS) test.out out/equilibrium_test.out out/hydrogen_chemistry_test.out out/thermochemistry_mode_test.out out/thermal_benchmark.out out/hii_region_test.out out/hii_region_benchmark.out out/molecular_test.out out/molecular_settings_benchmark.out out/non_equilibrium_test.out out/non_equilibrium_benchmark.out out/noneq_hot_driver.o out/noneq_hot_current.out out/noneq_hot_legacy.out out/noneq_hot_linear.out

tag:
	ctags $(TEST_SOURCES)

pkg:
	tar cvf plot_anime.tar $(TEST_SOURCES) Makefile tags
	gzip plot_anime.tar
