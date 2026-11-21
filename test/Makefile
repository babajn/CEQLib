#CC = mpicc
CC = gcc
CEQLIB_DIR ?= ../src
GSL_DIR ?= /home/babajn/lib/gsl
PYTHON ?= python3

OPTIMIZE += -O3 -std=c99
OPTIMIZE += -finline-functions -ffast-math -funroll-loops -frerun-cse-after-loop -fexpensive-optimizations -fomit-frame-pointer -msse3 -mfpmath=sse
OPTIMIZE += -Winline -Wundef -Wmissing-noreturn -Wstrict-prototypes 
OPTIMIZE += -g 
OPTIMIZE += -DDEBUG_LEVEL=0

INCLUDEPATH += -I$(GSL_DIR)/include -I$(CEQLIB_DIR)

LDFLAGS += -L$(CEQLIB_DIR) -lCEQLib
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

.PHONY: start ceqlib check benchmark hii hii-benchmark hii-plot co co-plot wolfire-plot
ceqlib:
	$(MAKE) -C $(CEQLIB_DIR) GSL_DIR=$(GSL_DIR) libCEQLib.a

start: $(OBJS) ceqlib
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) -o test.out

check: EquilibriumTest.o ceqlib
	mkdir -p out
	$(CC) $(CFLAGS) -fno-fast-math EquilibriumTest.o $(LDFLAGS) -o out/equilibrium_test.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/check.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/equilibrium_test.out" "$(abspath ../data/PEGASE2)"

EquilibriumTest.o: override CFLAGS += -fno-fast-math

benchmark: ThermalBenchmark.c ceqlib
	mkdir -p out
	$(CC) $(CFLAGS) -fno-fast-math -fopenmp ThermalBenchmark.c $(LDFLAGS) -o out/thermal_benchmark.out
	./out/thermal_benchmark.out atomic
	./out/thermal_benchmark.out molecular

hii: HIIRegionTest.c HIIRegionFixture.h ceqlib
	mkdir -p out/hii
	$(CC) $(CFLAGS) -fno-fast-math HIIRegionTest.c $(LDFLAGS) -o out/hii_region_test.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/hii.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/hii_region_test.out" "$(abspath ../data/PEGASE2)"; \
	mv equilibrium.csv thermal_curve.csv reactions.csv "$(CURDIR)/out/hii/"

hii-benchmark: HIIRegionBenchmark.c HIIRegionFixture.h ceqlib
	mkdir -p out/hii
	$(CC) $(CFLAGS) -fno-fast-math HIIRegionBenchmark.c $(LDFLAGS) -o out/hii_region_benchmark.out
	@set -e; csv=$$(mktemp "$(CURDIR)/out/hii/benchmark.XXXXXX"); \
	trap 'rm -f "$$csv"' EXIT; \
	./out/hii_region_benchmark.out "$(abspath ../data/PEGASE2)" > "$$csv"; \
	mv "$$csv" out/hii/hii_benchmark.csv
	@echo "H II region timings saved in out/hii/hii_benchmark.csv"

hii-plot: hii
	$(PYTHON) ../plot/plot_hii_region.py

co: MolecularTest.c ceqlib
	mkdir -p out/co
	$(CC) $(CFLAGS) -fno-fast-math MolecularTest.c $(LDFLAGS) -o out/molecular_test.out
	@set -e; run_dir=$$(mktemp -d "$(CURDIR)/out/co.XXXXXX"); \
	trap 'rm -rf "$$run_dir"' EXIT; cd "$$run_dir"; \
	"$(CURDIR)/out/molecular_test.out" "$(abspath ../data/PEGASE2)"; \
	mv equilibrium.csv column_scan.csv "$(CURDIR)/out/co/"

co-plot: co
	$(PYTHON) ../plot/plot_co_chemistry.py

wolfire-plot:
	$(PYTHON) ../plot/plot_wolfire1995.py

clean:
	rm -f $(OBJS) test.out EquilibriumTest.o out/equilibrium_test.out out/thermal_benchmark.out out/hii_region_test.out out/hii_region_benchmark.out out/molecular_test.out

tag:
	ctags *.[c,h]

pkg: 
	tar cvf plot_anime.tar ./*.[ch] ./*.rb Makefile tags
	gzip plot_anime.tar


main.o: main.c Makefile $(HEADS) 
