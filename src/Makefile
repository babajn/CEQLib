# Makefile for EQCLib
.SUFFIX: .c .o

CC = gcc
#CC = cc
GSL_DIR ?= $(or $(shell gsl-config --prefix 2>/dev/null),$(wildcard $(HOME)/lib/gsl),/usr)

CFLAGS = -std=c99 -O3 -g -Wall -ffast-math  #-fPIC
#CFLAGS += -I/home/babajn/lib/gsl-1.9/include 
#CFLAGS += -I/home/babajn/lib/dbg_malloc 
CPPFLAGS += -I$(GSL_DIR)/include

#GIT_VERSION = $(shell git log -1 --format="%H")
#CFLAGS += -DGITVERSION=\"$(GIT_VERSION)\"

LDLIBS += -L$(GSL_DIR)/lib -Wl,-rpath,$(GSL_DIR)/lib -lgsl -lgslcblas -lm

OBJS	= RandomNumberGenerator.o \
          Utilities.o \
		ChemicalEquilibriumCooling.o \
		ChemicalEquilibriumSolver.o \
		ChemistryMode.o \
		HydrogenChemistry.o \
		TemperatureTable.o \
		ThermalEquilibrium.o \
		ChemicalElements.o \
		EquilibriumAbundance.o \
		CollisionalGasPhaseChemistry.o \
		PhotochemicalGasPhaseChemistry.o \
		UVBackground.o \
		UVBackgroundSpectrum.o \
		MetalPhotoionization.o \
		MetalXray.o \
		PhotoelectronDeposition.o \
		MoleculeChemistry.o \
		GrainSurfaceChemistry.o \
		CosmicRay.o \
		Xray.o \
		Bremsstrahlung.o \
		CoolingCIE.o \
		HIIregionLineCooling.o \
		FineStructureLines.o \
		HyperFineStructureLines.o \
		StatisticalEquilibrium.o \

XRAY_OBJS = Xray.o ChemicalEquilibriumSolver.o ThermalEquilibrium.o
LEGACY_XRAY_OBJS = $(XRAY_OBJS:.o=LegacyXray.o)

SRCS 	= RandomNumberGenerator.c \
          Utilities.c \
		ChemicalEquilibriumCooling.c \
		ChemicalEquilibriumSolver.c \
		ChemistryMode.c \
		HydrogenChemistry.c \
		TemperatureTable.c \
		ThermalEquilibrium.c \
		ChemicalElements.c \
		EquilibriumAbundance.c \
		CollisionalGasPhaseChemistry.c \
		MoleculeChemistry.c \
		PhotochemicalGasPhaseChemistry.c \
		UVBackground.c \
		UVBackgroundSpectrum.c \
		MetalPhotoionization.c \
		MetalXray.c \
		PhotoelectronDeposition.c \
		GrainSurfaceChemistry.c \
		CosmicRay.c \
		Xray.c \
		Bremsstrahlung.c \
		CoolingCIE.c \
		HIIregionLineCooling.c \
		FineStructureLines.c \
		HyperFineStructureLines.c \
		StatisticalEquilibrium.c \
        
HEADS = Astro.h \
        Constants.h \
        RandomNumberGenerator.h \
        Utilities.h \
		ChemicalEquilibriumCooling.h \
		ChemicalElements.h \
		EquilibriumAbundance.h \
		CollisionalGasPhaseChemistry.h \
		PhotochemicalGasPhaseChemistry.h \
		UVBackground.h \
		MetalPhotoionization.h \
		MetalXray.h \
		PhotoelectronDeposition.h \
		MoleculeChemistry.h \
		GrainSurfaceChemistry.h \
		CosmicRay.h \
		Xray.h \
		Bremsstrahlung.h \
		CoolingCIE.h \
		CoolingHeatingRates.h \
		HIIregionLineCooling.h \
		HIIregionOxygenData.h \
		FineStructureLines.h \
		HyperFineStructureLines.h \
		StatisticalEquilibrium.h \


.PHONY: all
all:	depend libEQCLib.a

.PHONY: clean
clean:
	rm *.o *.s *.a *.swp depend.inc libEQCLib.a eqclib-cell -f

libEQCLib.a: $(OBJS) Makefile
	$(RM) libEQCLib.tmp.a
	$(AR) rcs libEQCLib.tmp.a $(OBJS)
	mv -f libEQCLib.tmp.a $@

libEQCLibLegacyXray.a: $(filter-out $(XRAY_OBJS),$(OBJS)) $(LEGACY_XRAY_OBJS) Makefile
	$(RM) $@.tmp
	$(AR) rcs $@.tmp $(filter %.o,$^)
	mv -f $@.tmp $@

eqclib-cell: CellRunner.o libEQCLib.a
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

# Use the production object's dependency chain and compile flags for both references.
ChemicalEquilibriumCooling.o: Makefile
ChemicalEquilibriumCoolingLegacy.o: ChemicalEquilibriumCooling.c ChemicalEquilibriumCooling.o
	$(COMPILE.c) -DEQC_LEGACY_NONEQUILIBRIUM_OXYGEN_COOLING $(OUTPUT_OPTION) $<
ChemicalEquilibriumCoolingLinearCIE.o: ChemicalEquilibriumCooling.c ChemicalEquilibriumCooling.o
	$(COMPILE.c) -DEQC_LINEAR_NONEQUILIBRIUM_CIE $(OUTPUT_OPTION) $<

$(LEGACY_XRAY_OBJS): %LegacyXray.o: %.c %.o Makefile
	$(COMPILE.c) -DEQC_LEGACY_XRAY_CUTOFF $(OUTPUT_OPTION) $<

# These modules check for non-finite values or return NaN for invalid trials.
CellRunner.o ChemicalEquilibriumCooling.o ChemicalEquilibriumCoolingLegacy.o ChemicalEquilibriumCoolingLinearCIE.o \
ChemicalEquilibriumSolver.o ChemistryMode.o HydrogenChemistry.o TemperatureTable.o ThermalEquilibrium.o \
EquilibriumAbundance.o MoleculeChemistry.o GrainSurfaceChemistry.o Xray.o PhotochemicalGasPhaseChemistry.o \
CoolingCIE.o CosmicRay.o HIIregionLineCooling.o UVBackground.o UVBackgroundSpectrum.o PhotoelectronDeposition.o MetalPhotoionization.o MetalXray.o \
$(LEGACY_XRAY_OBJS): override CFLAGS += -fno-fast-math

HIIregionLineCooling.o: HIIregionLineCooling.h HIIregionOxygenData.h

tag:
	ctags *.[c,h]

.PHONY: depend
depend: depend.inc

depend.inc: $(OBJS:.o=.c) CellRunner.c $(wildcard *.h) Makefile
	$(CC) $(CPPFLAGS) $(CFLAGS) -MM $(OBJS:.o=.c) CellRunner.c > $@.tmp && mv $@.tmp $@

ifneq ($(MAKECMDGOALS),clean)
-include depend.inc
endif
