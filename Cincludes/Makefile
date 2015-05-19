
#CPP=cxx
#CFLAGS=-O -lcomplex -lm
CPP=g++
F77=g77
AM_CPPFLAGS = $(shell lhapdf-config --cppflags) 
CFLAGS=-O3 -lm
LIBS=-lm -lstdc++ 
CCFLAGS= -O3 -c $(AM_CPPFLAGS)

.SUFFIXES: .o .C .cpp

all: random.o particle.o vegasGrid.o integration.o pdfqcd.o gauss.o

.C.o:
	$(CPP) $(CCFLAGS) $< 
.f.o:
	$(F77) -O -c $< 

clean:
	-rm *.o
