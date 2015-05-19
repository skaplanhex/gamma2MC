#ifndef DISTHDR
#define DISTHDR

#include <math.h>
#include <iostream>
#include <complex>
#include "../Cincludes/tdio.h"
#include "../Cincludes/pdfqcd.h"
#include "../Cincludes/particle.h"
#include "../Cincludes/integration.h"
#include "../Cincludes/random.h"
#include "../Cincludes/vegasGrid.h"
#include "../Cincludes/average.h"
#include "../Cincludes/histogram.h"

using namespace std;

#define PI M_PI

class Distribution {
 protected:

  histogram *h;
  double hsize;

 public:

  Distribution(int nbins=1, double hmin=0.0, double hmax=1.0) {
    hsize=(hmax-hmin)/((double) nbins);
    h = new histogram(hmin,hsize,nbins);
  }
  virtual ~Distribution(void) {}

  virtual double Param(particle p1, particle p2) {return 0.5;}

  void bin(int n, double *x, double *wgt) {h->bin(n,x,wgt);}

  void printDistribution(double norm) {h->printVegas(norm*hsize);}
  
};

// CosThetaStar distribution:

class CosThetaStar : public Distribution {

 public:

  CosThetaStar(int nBins=25, double hmin=0.0, double hmax=1.0) :
                      Distribution(nBins,hmin,hmax) {}
  virtual ~CosThetaStar(void) {}

  double Param(particle p1, particle p2);
  
};

// YStar distribution:

class YStar : public Distribution {

 public:

  YStar(int nBins=25, double hmin=0.0, double hmax=1.5) :
                      Distribution(nBins,hmin,hmax) {}
  virtual ~YStar(void) {}

  double Param(particle p1, particle p2);
  
};

// TanhYStar distribution:

class TanhYStar : public Distribution {

 public:

  TanhYStar(int nBins=25, double hmin=0.0, double hmax=1.0) :
                      Distribution(nBins,hmin,hmax) {}
  virtual ~TanhYStar(void) {}

  double Param(particle p1, particle p2);
  
};

// Ygg distribution:

class Ygg : public Distribution {

 public:

  Ygg(int nBins=20, double hmin=0.0, double hmax=2.8) :
                      Distribution(nBins,hmin,hmax) {}
  virtual ~Ygg(void) {}

  double Param(particle p1, particle p2);
  
};

// Phigg distribution:

class Phigg : public Distribution {

 public:

  Phigg(int nBins=5, double hmin=0.0, double hmax=PI+1.0E-8) :
                      Distribution(nBins,hmin,hmax) {}
  virtual ~Phigg(void) {}

  double Param(particle p1, particle p2);
  
};

// qT distribution:

class qT : public Distribution {

 public:

  qT(int nBins=20, double hmin=-1.0E-8, double hmax=40.0) :
                      Distribution(nBins,hmin,hmax) {}
  virtual ~qT(void) {}

  double Param(particle p1, particle p2);
  
};

// Diphoton Invariant Mass Distribution:

class Mgamgam : public Distribution {

 public:

  Mgamgam(int nBins=20, double hmin=80.0, double hmax=140.0) :
                      Distribution(nBins,hmin,hmax) {}
  virtual ~Mgamgam(void) {}

  double Param(particle p1, particle p2);
  
};

// User-Defined Distribution:

class distUser : public Distribution {

 public:

  distUser(int nBins=20, double hmin=80.0, double hmax=140.0) :
                      Distribution(nBins,hmin,hmax) {}
  virtual ~distUser(void) {}

  double Param(particle p1, particle p2);
  
};


#endif
