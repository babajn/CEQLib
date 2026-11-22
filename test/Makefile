#CC = mpicc
CC = gcc
EQCLIB_DIR ?= ../src
GSL_DIR ?= $(or $(shell gsl-config --prefix 2>/dev/null),$(wildcard $(HOME)/lib/gsl),/usr)
export PEGASE_DIR ?= $(abspath ../data/PEGASE2)
UVB_DIR ?= $(abspath ../data/FG20)
CHIMES_DATA_DIR ?= $(abspath ../../Chimes/chimes-data)
CLOUDY_DIR ?= $(abspath ../../c25.00)
BENCHMARK_CPU ?= $(shell awk '/^Cpus_allowed_list:/ {split($$2, cpus, /[,-]/); print cpus[1]}' /proc/self/status)
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
UNIT_OBJS = out/EquilibriumTest.o out/HydrogenChemistryTest.o out/ThermochemistryModeTest.o out/UVBackgroundTest.o
COMMON_HEADERS = common/Astro.h common/Constants.h
TEST_SOURCES = $(wildcard benchmark/*.[ch] co/*.[ch] common/*.[ch] hii/*.[ch] noneq/*.[ch] unit/*.[ch] uvb/*.[ch] wolfire1995/*.[ch])

ALL:   start

.PHONY: start eqclib check cell-check xray-legacy-check notebook-data benchmark hii hii-benchmark hii-plot co co-plot molecular-settings-benchmark molecular-settings-plot wolfire-plot noneq noneq-benchmark noneq-hot noneq-plot uvb uvb-metal uvb-metal-sensitivity uvb-metal-photoionization uvb-metal-photoionization-disk uvb-metal-cloudy uvb-charge-transfer-disk uvb-charge-transfer-bdf uvb-charge-transfer-benchmark uvb-benchmark uvb-asura-benchmark clean tag pkg
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
	$(CC) $(CFLAGS) -fno-fast-math out/UVBackgroundTest.o $(LDFLAGS) -o out/uv_background_test.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/check.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/equilibrium_test.out" "$$PEGASE_DIR"; \
	"$(CURDIR)/out/hydrogen_chemistry_test.out" "$$PEGASE_DIR"; \
	"$(CURDIR)/out/thermochemistry_mode_test.out" "$$PEGASE_DIR"; \
	"$(CURDIR)/out/uv_background_test.out" "$$PEGASE_DIR"

$(UNIT_OBJS): out/%.o: unit/%.c $(COMMON_HEADERS) $(EQCLIB_DIR)/ChemicalEquilibriumCooling.h Makefile | out
	$(CC) $(CFLAGS) -fno-fast-math -c $< -o $@

out/EquilibriumTest.o: $(EQCLIB_DIR)/Xray.h
out/HydrogenChemistryTest.o: $(EQCLIB_DIR)/CoolingCIE.h $(EQCLIB_DIR)/CoolingHeatingRates.h
out/UVBackgroundTest.o: $(EQCLIB_DIR)/CollisionalGasPhaseChemistry.h $(EQCLIB_DIR)/PhotochemicalGasPhaseChemistry.h \
                       $(EQCLIB_DIR)/EquilibriumAbundance.h $(EQCLIB_DIR)/MoleculeChemistry.h \
                       $(EQCLIB_DIR)/CoolingHeatingRates.h $(EQCLIB_DIR)/PhotoelectronDeposition.h $(EQCLIB_DIR)/UVBackground.h

xray-legacy-check: unit/EquilibriumTest.c eqclib | out
	$(MAKE) -C $(EQCLIB_DIR) GSL_DIR=$(GSL_DIR) libEQCLibLegacyXray.a
	$(CC) $(CFLAGS) -fno-fast-math -DEQC_LEGACY_XRAY_CUTOFF $< $(EQCLIB_DIR)/libEQCLibLegacyXray.a $(filter-out -lEQCLib,$(LDFLAGS)) -o out/equilibrium_legacy_xray_test.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/check.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/equilibrium_legacy_xray_test.out" "$$PEGASE_DIR"

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

out/uvb:
	mkdir -p $@

uvb: uvb/UVBackgroundReference.c uvb/compare.py eqclib | out/uvb
	$(CC) $(CFLAGS) -fno-fast-math $< $(LDFLAGS) -o out/uvb_reference.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/uvb.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/uvb_reference.out" "$(abspath $(UVB_DIR))"; \
	mv rates.csv equilibrium.csv soft-xray.csv secondary.csv disk-updates.csv disk-grid.csv EQCooling/UVBackgroundRates.dat "$(CURDIR)/out/uvb/"
	$(PYTHON) uvb/compare.py "$(abspath $(UVB_DIR))/fg20_cloudy.ascii" "$(CHIMES_DATA_DIR)" out/uvb

uvb-metal: uvb/MetalXrayReference.c uvb/metal_compare.py eqclib | out/uvb
	$(CC) $(CFLAGS) -fno-fast-math $< $(LDFLAGS) -o out/uvb_metal_reference.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/uvb-metal.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/uvb_metal_reference.out" "$(abspath $(UVB_DIR))"; \
	mv metal-*.csv "$(CURDIR)/out/uvb/"
	$(PYTHON) uvb/metal_compare.py "$(abspath $(UVB_DIR))/fg20_cloudy.ascii" out/uvb

uvb-metal-sensitivity: uvb/UVBackgroundReference.c eqclib | out/uvb
	$(CC) $(CFLAGS) -fno-fast-math $< $(LDFLAGS) -o out/uvb_metal_sensitivity.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/uvb-metal-sensitivity.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/uvb_metal_sensitivity.out" "$(abspath $(UVB_DIR))" --metal-sensitivity; \
	mv disk-500ev.csv "$(CURDIR)/out/uvb/"

out/uvb_metal_photoionization_reference.out: uvb/MetalPhotoionizationReference.c eqclib | out
	$(CC) $(CFLAGS) -fno-fast-math $< $(LDFLAGS) -o $@

uvb-metal-photoionization: out/uvb_metal_photoionization_reference.out uvb/metal_photoionization_compare.py | out/uvb
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/uvb-metal-pi.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/uvb_metal_photoionization_reference.out" "$(abspath $(UVB_DIR))"; \
	mv metal-photoionization-*.csv "$(CURDIR)/out/uvb/"
	$(PYTHON) uvb/metal_photoionization_compare.py "$(abspath $(UVB_DIR))/fg20_cloudy.ascii" "$(abspath $(CHIMES_DATA_DIR))" out/uvb --atomic-data "$(abspath $(CLOUDY_DIR))/data"

uvb-metal-cloudy: out/uvb_metal_photoionization_reference.out uvb/metal_cloudy.py | out/uvb
	$(PYTHON) uvb/metal_cloudy.py "$(abspath $(CLOUDY_DIR))" "$(abspath $(UVB_DIR))/fg20_cloudy.ascii" out/uvb "$(CURDIR)/out/uvb_metal_photoionization_reference.out"

uvb-metal-photoionization-disk: uvb/UVBackgroundReference.c eqclib | out/uvb
	$(CC) $(CFLAGS) -fno-fast-math $< $(LDFLAGS) -o out/uvb_metal_photoionization_disk.out
	$(CC) $(CFLAGS) -fno-fast-math -DEQC_EXPERIMENTAL_UVB_OXYGEN_COOLING $< $(EQCLIB_DIR)/ChemicalEquilibriumCooling.c $(LDFLAGS) -o out/uvb_metal_photoionization_oxygen.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/uvb-metal-disk.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/uvb_metal_photoionization_disk.out" "$(abspath $(UVB_DIR))" --metal-photoionization-disk; \
	mv disk-metal-photoionization.csv "$(CURDIR)/out/uvb/"; \
	"$(CURDIR)/out/uvb_metal_photoionization_oxygen.out" "$(abspath $(UVB_DIR))" --metal-photoionization-disk; \
	mv disk-metal-photoionization.csv "$(CURDIR)/out/uvb/disk-metal-photoionization-oxygen-weight.csv"

uvb-charge-transfer-disk: uvb/UVBackgroundReference.c eqclib | out/uvb
	$(CC) $(CFLAGS) -fno-fast-math $< $(LDFLAGS) -o out/uvb_charge_transfer_disk.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/uvb-charge-transfer-disk.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/uvb_charge_transfer_disk.out" "$(abspath $(UVB_DIR))" --charge-transfer-disk; \
	mv disk-charge-transfer.csv disk-charge-transfer-budget.csv "$(CURDIR)/out/uvb/"

uvb-charge-transfer-bdf: uvb/UVBChargeTransferReference.c eqclib | out/uvb
	$(CC) $(CFLAGS) -fno-fast-math $< $(LDFLAGS) -o out/uvb_charge_transfer_reference.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/uvb-charge-transfer-bdf.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/uvb_charge_transfer_reference.out" "$(abspath $(UVB_DIR))"; \
	"$(CURDIR)/out/uvb_charge_transfer_reference.out" "$(abspath $(UVB_DIR))" --regression "$$PEGASE_DIR"; \
	mv uvb-charge-transfer-bdf.csv uvb-charge-transfer-regression.csv "$(CURDIR)/out/uvb/"

uvb-charge-transfer-benchmark: uvb/ASURABenchmark.c eqclib | out/uvb
	$(CC) $(CFLAGS) -fno-fast-math $< $(LDFLAGS) -o out/uvb_charge_transfer_benchmark.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/uvb-charge-transfer-bench.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	taskset -c "$(BENCHMARK_CPU)" "$(CURDIR)/out/uvb_charge_transfer_benchmark.out" "$$PEGASE_DIR" "$(abspath $(UVB_DIR))" --charge-transfer; \
	taskset -c "$(BENCHMARK_CPU)" "$(CURDIR)/out/uvb_charge_transfer_benchmark.out" "$$PEGASE_DIR" "$(abspath $(UVB_DIR))" --charge-transfer-bidirectional; \
	mv asura-charge-transfer.csv asura-charge-transfer-bidirectional.csv "$(CURDIR)/out/uvb/"

uvb-benchmark: noneq/NonEquilibriumBenchmark.c uvb/benchmark.py eqclib | out/uvb
	$(CC) $(CFLAGS) -fno-fast-math $< $(LDFLAGS) -o out/uvb_benchmark.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/uvb-bench.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/uvb_benchmark.out" "$$PEGASE_DIR" --uvb-off; \
	mv benchmark.csv "$(CURDIR)/out/uvb/benchmark-off.csv"; \
	"$(CURDIR)/out/uvb_benchmark.out" "$$PEGASE_DIR" "$(abspath $(UVB_DIR))"; \
	mv benchmark.csv "$(CURDIR)/out/uvb/benchmark-on.csv"
	$(PYTHON) uvb/benchmark.py out/uvb

uvb-asura-benchmark: uvb/ASURABenchmark.c eqclib | out/uvb
	$(CC) $(CFLAGS) -fno-fast-math $< $(LDFLAGS) -o out/uvb_asura_benchmark.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/uvb-asura.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	taskset -c "$(BENCHMARK_CPU)" "$(CURDIR)/out/uvb_asura_benchmark.out" "$$PEGASE_DIR" "$(abspath $(UVB_DIR))"; \
	mv asura-benchmark.csv "$(CURDIR)/out/uvb/"
	@echo "ASURA benchmark saved in out/uvb/asura-benchmark.csv (CPU $(BENCHMARK_CPU))"

clean:
	rm -f $(WOLFIRE_OBJS) $(UNIT_OBJS) test.out out/equilibrium_test.out out/equilibrium_legacy_xray_test.out out/hydrogen_chemistry_test.out out/thermochemistry_mode_test.out out/uv_background_test.out out/uvb_reference.out out/uvb_benchmark.out out/uvb_asura_benchmark.out out/thermal_benchmark.out out/hii_region_test.out out/hii_region_benchmark.out out/molecular_test.out out/molecular_settings_benchmark.out out/non_equilibrium_test.out out/non_equilibrium_benchmark.out out/noneq_hot_driver.o out/noneq_hot_current.out out/noneq_hot_legacy.out out/noneq_hot_linear.out
	rm -f out/uvb_metal_reference.out out/uvb_metal_sensitivity.out

tag:
	ctags $(TEST_SOURCES)

pkg:
	tar cvf plot_anime.tar $(TEST_SOURCES) Makefile tags
	gzip plot_anime.tar
