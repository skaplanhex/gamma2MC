
#CPP=xlC
#CFLAGS=-O -lcomplex -lm
CPP=g++
F77=f77
AM_CPPFLAGS = $(shell lhapdf-config --cppflags) 
AM_LDFLAGS = $(shell lhapdf-config --ldflags) 
CFLAGS= -O3 
CFLAGS=
#LIBS=-lm -lstdc++ -lg2c $(AM_LDFLAGS)
LIBS=-lm -L/afs/cern.ch/sw/lcg/external/gcc/4.7.0/x86_64-slc6-gcc47-opt/lib64 -lstdc++ -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6 -lg2c $(AM_LDFLAGS)
CCFLAGS= -O3 -c  $(AM_CPPFLAGS)
#CCFLAGS= -c

.SUFFIXES: .o .C .f

INCL=../Cincludes
CINCLUDES= $(INCL)/random.o $(INCL)/particle.o $(INCL)/vegasGrid.o $(INCL)/integration.o $(INCL)/pdfqcd.o $(INCL)/gauss.o
BASE= gamma2NLO.o HiggsNLO.o ggNLO.o qqNLO.o Cuts.o Distribution.o

gamma2MC: gamma2MC.o $(BASE) Makefile
	$(CPP) $(CFLAGS) -o gamma2MC gamma2MC.o $(BASE) $(CINCLUDES) $(LIBS)

.C.o:
	$(CPP) $(CCFLAGS) $< 
.f.o:
	$(F77) $(CCFLAGS) $< 

clean:
	-rm *.o


