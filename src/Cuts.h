#ifndef CUTSHDR
#define CUTSHDR

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

class Cuts {
 public:
  virtual int cut(particle k1, particle k2)=0;
  virtual int cut(particle k1, particle k2, particle p)=0;
};

// No Cuts:

class cutNone : public Cuts {

 public:

  cutNone(void) {}
  virtual ~cutNone(void) {}

  int cut(particle k1, particle k2);
  int cut(particle k1, particle k2, particle p);
  
};

// Rapidity cuts on gamma-gamma system:

class cutHiggs : public Cuts {
 protected:
  double yHcut;

 public:

  cutHiggs(void) : Cuts() {yHcut=2.0;}
  virtual ~cutHiggs(void) {}

  void setYcut(double y) {yHcut=y;}

  int cut(particle k1, particle k2);
  int cut(particle k1, particle k2, particle p);
  
};

// Minimal pt, rapidity cuts on photons, no Isolation cut:

class cutPhoton : public Cuts {
 protected:
  double ycut, ptcut2, ptcut1;  

  virtual int cutIsolation(particle p1, particle p2, particle p3);

 public:

  cutPhoton(void) {ycut=2.5; ptcut2=25.0; ptcut1=40.0;}

  virtual ~cutPhoton(void) {}

  void setPhotonCuts(double y, double p2, double p1) {
                         ycut=y; ptcut2=p2; ptcut1=p1;}

  int cut(particle k1, particle k2);
  int cut(particle k1, particle k2, particle p);
  
};

// Photon cuts plus standard Isolation:

class cutStandard : public cutPhoton {
 protected:
  double Rcut, Etcut;  

  virtual int cutIsolation(particle p1, particle p2, particle p3);

 public:

  cutStandard(void) : cutPhoton() 
                         {Rcut=0.4; Etcut=5.0;}

  virtual ~cutStandard(void) {}

  void setIsolation(double R, double E) {Rcut=R; Etcut=E;}

};

// Standard cuts plus an observed parton pt cut:
// (If there is no parton, then it is automatically cut.)

class cutPt : public cutStandard {
 protected:
  double ptgcut;

 public:

  cutPt(void) : cutStandard() 
                         {ptgcut=10.0;}

  virtual ~cutPt(void) {}

  void setPtCut(double p) {ptgcut=p;}

  int cut(particle k1, particle k2);
  int cut(particle k1, particle k2, particle p);
  
};

// Photon cuts plus Frixione (smooth) Isolation:

class cutFrixione : public cutPhoton {
 protected:
  double Rcut, Ccut, epsilon;  

  int cutIsolation(particle p1, particle p2, particle p3);

 public:

  cutFrixione(void) : cutPhoton() 
                         {Rcut=0.4; Ccut=1.0-cos(Rcut);epsilon=1.0;}

  virtual ~cutFrixione(void) {}

  void setIsolation(double R, double E) {Rcut=R; Ccut=1.0-cos(Rcut);
                                 epsilon=E;}

};

// Photon cuts plus D0 Isolation cuts (similar to Frixione):

class cutD0 : public cutPhoton {
 protected:
  double Rcut, epsilon;  

  int cutIsolation(particle p1, particle p2, particle p3);

 public:

  cutD0(void) : cutPhoton() 
                         {Rcut=0.4; epsilon=.07;}

  virtual ~cutD0(void) {}

  void setIsolation(double R, double E) {Rcut=R; epsilon=E;}

};

// Standard cuts plus a cut 
// on jets in a larger annulus around the photon.

class cutJet : public cutStandard {
 protected:
  double Rjet, Etjet;
  // Note: Rjet is not the size of the cone around the jet,
  //       but the distance between the jet and the photon
  //       in which the jet may be cut.
  //       At NLO, there can be at most one colored parton in
  //       the final-state, so there can be no dependence on
  //       the jet cone size.  At NNLO and beyond, this must
  //       be corrected.
  //
  // It is assumed that Rjet>Rcut and Etjet>Etcut.  

  int cutIsolation(particle p1, particle p2, particle p3);

 public:

  cutJet(void) : cutStandard() 
                         {Rjet=2.0; Etjet=30.0;}

  virtual ~cutJet(void) {}

  void setJetVeto(double R, double E) {Rjet=R; Etjet=E;}

};

// Cuts to only allow events with a parton (E>Etjet) in an annulus
// between Rcut and Rjet.
// (This is what must be subtracted from Standard Cuts to give
// the Jet Veto.)

class cutAnnulus : public cutJet {

 public:

  cutAnnulus(void) : cutJet() {} 

  virtual ~cutAnnulus(void) {}

  int cut(particle k1, particle k2);
  int cut(particle k1, particle k2, particle p);
  
};

// Photon cuts plus User-defined Isolation Cuts.
// This class is a place-holder, for a user of the code to write a new set 
// of isolation cuts, which may depend on up to 4 parameters.
// Currently it is exactly the same as the Photon cuts plus standard Isolation
// (i.e., "cutStandard" with P1=Rcut, P2=Etcut, P3 and P4 ignored)

// The user needs to modify the class subroutine "cutUser::cutIsolation"
// in the file Cuts.C as desired.

class cutUser : public cutPhoton {
 protected:
  double P1, P2, P3, P4;  

  virtual int cutIsolation(particle p1, particle p2, particle p3);

 public:

  cutUser(void) : cutPhoton() 
                         {P1=0.4; P2=5.0; P3=0.0; P4=0.0;}

  virtual ~cutUser(void) {}

  void setIsolation(double Pa1, double Pa2, double Pa3, double Pa4) {
                        P1=Pa1; P2=Pa2; P3=Pa3; P4=Pa4;}

};

#endif
