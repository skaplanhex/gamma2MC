#ifndef QQNLOHDR
#define QQNLOHDR 

#include "gamma2NLO.h"

class qqNLO : public gamma2NLO {
 protected:
 
     virtual void setSigma0(void) { sigma0=hbarCsqPB*PI/Qsq/Qsq/3.0*
				      alpha*alpha;}

     //  Note that an extra Bose factor of 1/2 is included in sigma0.
     //  Thus, we will integrate over the full phase space of k1 and k2
     //  as if they were distinguishable.

     particle k1f, k2f, pf;  // fragmentation particles.
     // In this definition, Bose symmetry was used to allow
     // photon 1 to always be the collinear photon.


     double fqq, fqq0, fgq, fqg;

     double pqg, pqq, pgammaq;

     double sigLOqg, sigLOgq;

     double x1,y;  //     x1=t/s,  y=u/s
    
     double X1, Y;  //    X1=log(-t/s),  Y=log(-u/s) 

     virtual void initializeSubprocess(void);

     double xwgtSf, xwgtVf;

 public:

     qqNLO(double energy=14000.0, collider coll=pp,
	   double m=115.0, event_type et=NLO, int np=0) :
                     gamma2NLO(energy,coll,m,et,np) 
                      {qqflag=QQTOTAL; setSigma0();}

     virtual ~qqNLO(void) {}

     virtual double wgtLO(void);
     virtual double wgtV(void);
     virtual double wgtVplus(void);
     virtual double wgtVa(void);
     virtual double wgtVb(void);
     virtual double wgtR(void);
     virtual double wgtSa(void);
     virtual double wgtSb(void);
     virtual double wgtSf(void);
     virtual double wgtVf(void);

     virtual double fxn(void);
     virtual void fillBins(void);
};

#endif


