#CC = mpicc
CC = gcc

OPTIMIZE += -O3 -std=c99
OPTIMIZE += -finline-functions -ffast-math -funroll-loops -frerun-cse-after-loop -fexpensive-optimizations -fomit-frame-pointer -msse3 -mfpmath=sse
OPTIMIZE += -Winline -Wundef -Wmissing-noreturn -Wstrict-prototypes 
OPTIMIZE += -g 
OPTIMIZE += -DDEBUG_LEVEL=0

INCLUDEPATH += -I/home/babajn/lib/gsl/include 
INCLUDEPATH += -I/home/babajn/lib/libCEQ/src 

LDFLAGS += -L/usr/lib64 -lX11 -lgfortran -lm
LDFLAGS += -L/home/babajn/lib/gsl/lib -lgsl -lgslcblas 
LDFLAGS += -L/home/babajn/lib/libCEQ/src -lCEQLib

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

start: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) -o test.out

clean:
	rm $(OBJS) test.out

tag:
	ctags *.[c,h]

pkg: 
	tar cvf plot_anime.tar ./*.[ch] ./*.rb Makefile tags
	gzip plot_anime.tar


main.o: main.c Makefile $(HEADS) 
