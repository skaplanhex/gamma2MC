
#CPP=cxx
#CFLAGS=-O -lcomplex -lm
CPP=g++
F77=g77
AM_CPPFLAGS = $(shell lhapdf-config --cppflags) 
CFLAGS=-O3 -lm
#LIBS=-lm -lstdc++ 
LIBS=-lm -L/afs/cern.ch/sw/lcg/external/gcc/4.7.0/x86_64-slc6-gcc47-opt/lib64 -lstdc++ -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6 -lg2c $(AM_LDFLAGS)
CCFLAGS= -O3 -c $(AM_CPPFLAGS)

.SUFFIXES: .o .C .cpp

all: random.o particle.o vegasGrid.o integration.o pdfqcd.o gauss.o

.C.o:
	$(CPP) $(CCFLAGS) $< 
.f.o:
	$(F77) -O -c $< 

clean:
	-rm *.o
