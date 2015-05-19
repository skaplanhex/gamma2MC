#ifndef HIGGSNLOHDR
#define HIGGSNLOHDR 

#include "gamma2NLO.h"

class HiggsNLO : public gamma2NLO {
 protected:
 
     virtual void setSigma0(void) 
                        {sigma0=hbarCsqPB*PI/64.0*pow(2.0*as2pi/3.0/vev,2);
                           sigLO=sigma0*A1sq()*9.0/4.0;}

     complex<double> W2(double s);

     double A1sq(void);

     double fgg, fgq, fqg, fqiqj, fqiqi, fqiai, fgg0;

     double pgg, pgq;

     virtual void initializeSubprocess(void);
     virtual double wgtLO(void);
     virtual double wgtV(void);
     virtual double wgtVplus(void);
     virtual double wgtVa(void);
     virtual double wgtVb(void);
     virtual double wgtR(void);
     virtual double wgtSa(void);
     virtual double wgtSb(void);

 public:

     HiggsNLO(double energy=14000.0, collider coll=pp,
              double m=115.0, event_type et=NLO, int np=0) :
                     gamma2NLO(energy,coll,m,et,np) {}

     virtual ~HiggsNLO(void) {}

     virtual double wgtnorm(void) {return 1000.0/2.0/Q;}
  
     double XsecLO(void);
     double XsecNLO(void);

   //Note that the NLO cross sections above only contain mt-dependence in
   //the overall sigma0 term, but not in the NLO coefficients.

     double XsecLOlargeMt(void);
     double XsecNLOlargeMt(void);

};

#endif


