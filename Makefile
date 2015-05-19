
#CPP=xlC
#CFLAGS=-O -lcomplex -lm
CPP=g++
F77=f77
AM_CPPFLAGS = $(shell lhapdf-config --cppflags) 
AM_LDFLAGS = $(shell lhapdf-config --ldflags) 
CFLAGS= -O 
CFLAGS=
LIBS=-lm -lstdc++ -lg2c $(AM_LDFLAGS)
CCFLAGS= -O -c $(AM_CPPFLAGS)
#CCFLAGS= -c

.SUFFIXES: .o .C .f

INCL=./Cincludes
SRC=./src
BIN=./bin

install:  
	cd $(INCL); make all
	cd $(SRC); make gamma2MC
	-mv $(SRC)/gamma2MC $(BIN)/gamma2MC

.C.o:
	$(CPP) $(CCFLAGS) $< 
.f.o:
	$(F77) $(CCFLAGS) $< 

clean:
	cd $(INCL); make clean
	cd $(SRC); make clean
	cd $(BIN); rm gamma2MC


