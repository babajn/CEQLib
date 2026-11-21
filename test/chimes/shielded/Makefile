# EQCLib vs CHIMES comparison.
#   make -C test/chimes/shielded CHIMES_DIR=/path/to/Chimes run
# CHIMES runs without > 13.6 eV photons (EQCLib uses NLyn = 0); the *-euv runs keep
# the full B87 spectrum as a control.
# CHIMES_DIR holds chimes/ (built libchimes.so), chimes-data/ and .local/ (HDF5, SUNDIALS).
CC = gcc
PYTHON ?= python3
GSL_DIR ?= $(HOME)/lib/gsl
EQCLIB_DIR ?= $(abspath ../../../src)
EQCLIB_GAS_Z ?= ASPLUND_2009_SOLAR_METALLICITY
CHIMES_DIR ?= $(abspath ../../../../Chimes)
CHIMES_LIB_DIR = $(CHIMES_DIR)/chimes
CHIMES_DATA ?= $(CHIMES_DIR)/chimes-data
CHIMES_LOCAL ?= $(CHIMES_DIR)/.local
PEGASE_DIR ?= $(abspath ../../../data/PEGASE2)
OUT ?= ../../out/chimes/shielded
BENCH_CPU ?= 2

CFLAGS = -std=c99 -O2 -Wall -Wextra
EQC_CFLAGS = $(CFLAGS) -DCMP_GAS_Z=$(EQCLIB_GAS_Z) -I$(EQCLIB_DIR) -I$(GSL_DIR)/include
EQC_LIBS = -L$(EQCLIB_DIR) -lEQCLib -L$(GSL_DIR)/lib -Wl,-rpath,$(GSL_DIR)/lib -lgsl -lgslcblas -lm
CHIMES_CFLAGS = $(CFLAGS) -DCHIMES_USE_DOUBLE_PRECISION -I$(CHIMES_LIB_DIR)/src -I$(CHIMES_LOCAL)/include
CHIMES_LIBS = -L$(CHIMES_LIB_DIR) -Wl,-rpath,$(CHIMES_LIB_DIR) -lchimes \
              -L$(CHIMES_LOCAL)/lib -Wl,-rpath,$(CHIMES_LOCAL)/lib -lhdf5 -lm

.PHONY: all eqclib run equilibrium timing courant cooling plot clean
all: $(OUT)/compare_eqclib $(OUT)/compare_chimes

eqclib:
	$(MAKE) -C $(EQCLIB_DIR) GSL_DIR=$(GSL_DIR) libEQCLib.a

$(OUT)/compare_eqclib: CompareDriver.c EQCLibBackend.c CompareBackend.h Makefile eqclib
	mkdir -p $(OUT)
	$(CC) $(EQC_CFLAGS) CompareDriver.c EQCLibBackend.c $(EQC_LIBS) -o $@

$(OUT)/compare_chimes: CompareDriver.c ChimesBackend.c CompareBackend.h Makefile
	mkdir -p $(OUT)
	$(CC) $(CHIMES_CFLAGS) CompareDriver.c ChimesBackend.c $(CHIMES_LIBS) -o $@

$(OUT)/h2_eqclib: H2CompareDriver.c H2Diagnostics.h EQCLibBackend.c CompareBackend.h Makefile eqclib
	mkdir -p $(OUT)
	$(CC) $(EQC_CFLAGS) -DEQC_H2_DIAGNOSTICS H2CompareDriver.c EQCLibBackend.c $(EQC_LIBS) -o $@

$(OUT)/h2_chimes: H2CompareDriver.c H2Diagnostics.h ChimesBackend.c CompareBackend.h Makefile
	mkdir -p $(OUT)
	$(CC) $(H2_CHIMES_CPPFLAGS) $(CHIMES_CFLAGS) -DEQC_H2_DIAGNOSTICS H2CompareDriver.c ChimesBackend.c $(CHIMES_LIBS) -o $@

$(OUT)/thermal_eqclib: ThermalBudgetDriver.c ThermalDiagnostics.h H2Diagnostics.h EQCLibBackend.c CompareBackend.h Makefile eqclib
	mkdir -p $(OUT)
	$(CC) $(EQC_CFLAGS) -DEQC_H2_DIAGNOSTICS -DEQC_THERMAL_DIAGNOSTICS ThermalBudgetDriver.c EQCLibBackend.c $(EQC_LIBS) -o $@

$(OUT)/thermal_chimes: ThermalBudgetDriver.c ThermalDiagnostics.h H2Diagnostics.h ChimesBackend.c CompareBackend.h Makefile
	mkdir -p $(OUT)
	$(CC) $(CHIMES_CFLAGS) -DCMP_CHIMES_DOPPLER_B=3.0 -DEQC_H2_DIAGNOSTICS -DEQC_THERMAL_DIAGNOSTICS ThermalBudgetDriver.c ChimesBackend.c $(CHIMES_LIBS) -o $@

# Accuracy runs may run in parallel; timing runs one after another on one core.
equilibrium: all
	cd $(OUT) && ./compare_eqclib equilibrium default $(PEGASE_DIR) eq_eqclib.csv
	cd $(OUT) && ./compare_eqclib equilibrium co $(PEGASE_DIR) eq_eqclib_co.csv
	cd $(OUT) && ./compare_chimes equilibrium matched $(CHIMES_DATA) eq_chimes_matched.csv > /dev/null
	cd $(OUT) && ./compare_chimes equilibrium full $(CHIMES_DATA) eq_chimes_full.csv > /dev/null
	cd $(OUT) && ./compare_chimes equilibrium matched-euv $(CHIMES_DATA) eq_chimes_matched_euv.csv > /dev/null
	cd $(OUT) && ./compare_chimes equilibrium full-euv $(CHIMES_DATA) eq_chimes_full_euv.csv > /dev/null

timing: all
	cd $(OUT) && taskset -c $(BENCH_CPU) ./compare_eqclib timing default $(PEGASE_DIR) timing_eqclib.csv
	cd $(OUT) && taskset -c $(BENCH_CPU) ./compare_eqclib timing co $(PEGASE_DIR) timing_eqclib_co.csv
	cd $(OUT) && taskset -c $(BENCH_CPU) ./compare_chimes timing matched $(CHIMES_DATA) timing_chimes_matched.csv > /dev/null
	cd $(OUT) && taskset -c $(BENCH_CPU) ./compare_chimes timing full $(CHIMES_DATA) timing_chimes_full.csv > /dev/null

courant: all
	cd $(OUT) && taskset -c $(BENCH_CPU) ./compare_eqclib courant default $(PEGASE_DIR) courant_eqclib.csv
	cd $(OUT) && taskset -c $(BENCH_CPU) ./compare_chimes courant matched $(CHIMES_DATA) courant_chimes_matched.csv > /dev/null
	cd $(OUT) && taskset -c $(BENCH_CPU) ./compare_chimes courant full $(CHIMES_DATA) courant_chimes_full.csv > /dev/null

cooling: all
	cd $(OUT) && ./compare_eqclib cooling default $(PEGASE_DIR) cooling_eqclib.csv
	cd $(OUT) && ./compare_chimes cooling matched $(CHIMES_DATA) cooling_chimes_matched.csv > /dev/null
	cd $(OUT) && ./compare_chimes cooling full $(CHIMES_DATA) cooling_chimes_full.csv > /dev/null

run: equilibrium timing courant cooling

plot:
	$(PYTHON) ../../../plot/plot_chimes_comparison.py --data-dir "$(OUT)" --out-dir "$(OUT)"

clean:
	rm -f $(OUT)/compare_eqclib $(OUT)/compare_chimes
