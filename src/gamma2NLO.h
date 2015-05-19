#ifndef GAMMA2NLOHDR
#define GAMMA2NLOHDR 

/*

PROGRAM:  gamma2NLO,  version 1.1

DATE:     September 25, 2001  
          Version 1.0 -- Updated October 25, 2007
	  Version 1.1 -- Updated January 24, 2011 
	              -- (Now uses LHAPDF)

AUTHOR:   Carl Schmidt    

LANGUAGE: C++

   A.  DESCRIPTION OF THE PROGRAM

Gamma Gamma production Monte Carlo in pp (or ppbar) collider
This file contains the generic program for
p p(bar) --> gamma gamma X at NLO.
It can be used to derive the production via different processes,
including

1) Through Higgs resonance: p p(bar) --> H X, followed by H--> gamma gamma.
   (Gives cross sections in fb.  Branching ratio to gamma gamma is
    not included--yet.)

2) Through q q(bar) and its associated processes through NLO:
    q (qbar) --> gamma gamma X.
   (Gives cross sections in fb.)

   ******* Subtraction of the final state photon singularities has
   been done, but no fragmentation contributions have been included.

3) Through the (formally NNLO) process g g -->  gamma gamma X,
   which occurs at one loop.  (This box contribution is then treated
   as the LO part of a NLO calculation.)
   (Gives cross sections in fb.)


   B.  HEADER FILES USED 

This header file should be included by the driver program with the
following command:

       ==>    #include "gamma2NLO.h"

It in turn includes other standard and additional 
header files.
    
To run the program, the files "gamma2NLO.C", and the
driver program and other necessary ".C" files should be compiled and 
linked together with the math library.

   C.  FURTHER INFORMATION

See the class definitions below for more information.

*/

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
#include "Cuts.h"
#include "Distribution.h"

using namespace std;

#define MXDIM 50
#define PI M_PI

const double Nc=3.0;
const double CF=4.0/3.0;
const double Zeta3=1.202056903;
const complex<double> I(0.0,1.0);
const double nf=5.0;
const double beta0=(11.0*Nc-2.0*nf)/3.0;

enum event_type { LO = 0, NLO = 1, EASYNLO = 2, HARDNLO = 3 };
enum qq_type { QQTOTAL = 0, QQONLY = 1, QGONLY = 2 };

class gamma2NLO : public vegasGrid {
 protected:
 
     event_type evt; // LO = 0, NLO = 1, EASYNLO = 2, HARDNLO = 3

     collider col;  //(col==pp or col==ppbar)

     double rs, ss;       // rs=sqrt(s)

     double Qsq, Q;    // Qsq is the gamma gamma invariant mass-squared
                          // Q =sqrt(Qsq)

     double Q0, Qbin;  // Q0 is the center of a bin in Q of size Qbin.

     double chiR, chiF, chifr;
 
     double muR, muF, mufr;   // muR=chiR*Q, muF=chiF*Q, mufr=chifr*Q

     double tau;        // tau=Qsq/ss;

     //  The following are the variables for the photons:

     double w1, phi1;   //   w1=(1-cos(theta))/2 in gamma-gamma
                        //          COM frame.
                        //   phi1 = azimuthal angle

     //  The following are the variables for the radiated gluon:

     double z, w, phi;   //   z=tau/xA/xB 
                         //   w=(1-cos(theta))/2 in gluon-Q
                         //      COM frame.
                         //   phi= azimuthal angle

     double zm;           //   =1.0-z;

     double wm;           //   =1.0-w;

     double w1m;          //   =1.0-w1;

     double logzmin;      //   =log(1-zmin);
                          //   To be set in gggamMC.

     double logHF;        //   =log(Qsq/muF/muF);

     double logHfr;        //   =log(Qsq/mufr/mufr);

     double logRF;        //   =log(muR*muR/muF/muF);

     double logz;         //   =log(z); 

     double logzm;        //   =log(zm);

     double xA, xB;  // PDF momentum fractions

     double xA0, xB0;     

     double as2pi;   // alphas(Q)/2/Pi

     double sigma0;   // A universal prefactor.

     double sigLO;    // The LO cross section.

     virtual void setSigma0(void)=0;

     particle k10, k20;  // photons with 2->2 kinematics

     particle k11, k21, p1;  // photons and gluon/quark with 2->3 kinematics

     particle k1A, k2A;  // photons with 2->2 kinematics (and boost A)

     particle k1B, k2B;  // photons with 2->2 kinematics (and boost B)

     virtual void initializeSubprocess(void)=0;

     complex<double> Spa[6][6];

     complex<double> Spb[6][6];

     double S[6][6];

     double jacob, jacob0;

     void setParameters(void);
     void buildProducts(void);
     void buildParticles(void);

     double xwgtLO, xwgtV, xwgtVplus, xwgtVa, xwgtVb,
            xwgtSa, xwgtSb, xwgtR;

     virtual double wgtLO(void)=0;
     virtual double wgtV(void)=0;
     virtual double wgtVplus(void)=0;
     virtual double wgtVa(void)=0;
     virtual double wgtVb(void)=0;
     virtual double wgtR(void)=0;
     virtual double wgtSa(void)=0;
     virtual double wgtSb(void)=0;

     Cuts *cc;  // Default is no cuts.

     Distribution *dist;

     qq_type qqflag;

     virtual void fillBins(void);

     void initialize(void) { setParameters();
                             buildProducts();
                             buildParticles();
                             initializeSubprocess();
                           }

 public:

     gamma2NLO(double energy=14000.0, collider coll=pp, 
	       double m=115.0, event_type et=NLO,  
	       int np=0) : vegasGrid(np) { evt=et; 
			     rs=energy; ss=rs*rs; col=coll; Q0=m; Qbin=0.0;
                           chiR=chiF=chifr=1.0; 
			   cc = new cutNone;
                           dist = new Distribution;
                           }

     virtual ~gamma2NLO(void) {}

     void ppbarCollider(void) {col=ppbar;}

     void ppCollider(void) {col=pp;}

     void newChiR(double cc) {chiR=cc; muR=chiR*Q; }

     void newChiF(double cc) {chiF=cc; muF=chiF*Q; }

     void newChifr(double cc) {chifr=cc; mufr=chifr*Q; }

     void newQ(double m) {Q0=m; }

     void newQbin(double dQ) {Qbin=dQ;}

     void setQQType(qq_type et) { qqflag=et;}

     void printDistribution(double norm) 
       {dist->printDistribution(norm);}

     void setEventType(event_type et) { evt=et;}
                       
     virtual double wgtnorm(void) {return 1000.0*Qbin;}  
// = 1000.0*Qbin so that output is fb, 
//  For Higgs, = 1000.0/2.0/Q, so output is fb.




     //  wgtLO, wgtV are integrated over xA, w1, phi1 only.
     //  wgtVplus, wgtVa, wgtVb are integrated over xA, w1, phi1, and z.
     //  wgtR, wgtSa, wgtSb are integrated over xA, w1, phi1, z, w, phi.

     //  wgtLO, wgtV, wgtVplus have 2->1 kinematics.
     //  wgtVa, wgtSa have 2->1 kinematics but with a z-boost along
     //              parton A.
     //  wgtVb, wgtSb have 2->1 kinematics but with a z-boost along
     //              parton B.
     //  wgtR has full 2->2 kinematics.

     virtual double fxn(void);

     // On integration, fxn() returns Integral(dsigma/dQ,dQ)
     // in fb.  Except for production through a Higgs boson.

     // For Higgs Production:  YOU MUST SET QBIN=0. Q0=HIGGS MASS
     // (The Higgs width is assumed neglibible.)
     // However, cuts on the gamma-gamma in the decay
     // can be applied.
     // Then on integration, fxn() returns sigma in fb. 
     
     virtual void processEvent(void) {
                            double xwgt=fxn();
                            xsec.add(vwgt*xwgt);
                            fillBins();
                            }

     average xsec;

     void reset(void) {xsec.reset();} 
    
     void setCuts(Cuts *c) {cc=c;}

     void setDistribution(Distribution *d) {dist=d;}
};

#endif


