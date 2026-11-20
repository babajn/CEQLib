# Makefile for CEQLib
.SUFFIX: .c .o

CC = gcc
#CC = cc

CFLAGS = -std=c99 -O3 -g -Wall -ffast-math  #-fPIC
#CFLAGS += -I/home/babajn/lib/gsl-1.9/include 
#CFLAGS += -I/home/babajn/lib/dbg_malloc 
CFLAGS += -I/home/babajn/lib/gsl/include 

#GIT_VERSION = $(shell git log -1 --format="%H")
#CFLAGS += -DGITVERSION=\"$(GIT_VERSION)\"

LDFLAGS += -lm
LDFLAGS += -L/usr/lib64 -lX11 -lgfortran -lm
#LDFLAGS += -L/home/babajn/lib/gsl-1.9/lib -lgsl -lgslcblas 
#LDFLAGS += -L/home/babajn/lib/dbg_malloc -ldbg_malloc
LDFLAGS += -L/home/babajn/lib/gsl/lib -lgsl -lgslcblas 

OBJS	= RandomNumberGenerator.o \
          Utilities.o \
		  nrutil.o \
		  NewtonRaphson.o \
		ChemicalEquilibriumCooling.o \
		ChemicalElements.o \
		EquilibriumAbundance.o \
		CollisionalGasPhaseChemistry.o \
		PhotochemicalGasPhaseChemistry.o \
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

SRCS 	= RandomNumberGenerator.c \
          Utilities.c \
		  nrutil.c \
		  NewtonRaphson.c \
		ChemicalEquilibriumCooling.c \
		ChemicalElements.c \
		EquilibriumAbundance.c \
		CollisionalGasPhaseChemistry.c \
		MoleculeChemistry.c \
		PhotochemicalGasPhaseChemistry.c \
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
		nrutil.h \
		NewtonRaphson.h \
		ChemicalEquilibriumCooling.h \
		ChemicalElements.h \
		EquilibriumAbundance.h \
		CollisionalGasPhaseChemistry.h \
		PhotochemicalGasPhaseChemistry.h \
		MoleculeChemistry.h \
		GrainSurfaceChemistry.h \
		CosmicRay.h \
		Xray.h \
		Bremsstrahlung.h \
		CoolingCIE.h \
		HIIregionLineCooling.h \
		FineStructureLines.h \
		HyperFineStructureLines.h \
		StatisticalEquilibrium.h \


.PHONY: all
all:	depend libCEQLib.a

.PHONY: clean
clean:
	rm *.o *.s *.a *.swp depend.inc libCEQLib.a -f

libCEQLib.a: $(OBJS)
	$(AR) -r $@ $^
	ranlib $@

tag:
	ctags *.[c,h]

.PHONY: depend
depend: $(OBJS:.o=.c)
	-@ $(RM) depend.inc
	-@ for i in $^; do cpp -MM $$i | sed "s/\ [_a-zA-Z0-9][_a-zA-Z0-9]*\.c//g" >> depend.inc; done

-include depend.inc
