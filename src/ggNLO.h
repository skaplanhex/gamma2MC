#ifndef GGNLOHDR
#define GGNLOHDR 

#include "gamma2NLO.h"

class ggNLO : public gamma2NLO {
 protected:
 
     virtual void setSigma0(void) { sigma0=hbarCsqPB*PI/Qsq/Qsq/16.0*
			   pow(alpha*as2pi*(3.0*(1.0/9.0)+2.0*(4.0/9.0)),2);}

     //  Note that an extra Bose factor of 1/2 is included in sigma0.
     //  Thus, we will integrate over the full phase space of k1 and k2
     //  as if they were distinguishable.

     double fgg, fgg0;

     double pgg;

     double x1,y;  //     x1=t/s,  y=u/s
    
     double X1, Y;  //    X1=log(-t/s),  Y=log(-u/s) 

     complex<double> Mpppp,Mmppp,Mmmpp,Mmpmp,Mmppm;

     complex<double> Mppppp(int i1, int i2, int i3, int i4, int i5);
     complex<double> Mmpppp(int i1, int i2, int i3, int i4, int i5);
     complex<double> Mmmppp(int i1, int i2, int i3, int i4, int i5);
     complex<double> Mmpmpp(int i1, int i2, int i3, int i4, int i5);
     complex<double> Mmmmmm(int i1, int i2, int i3, int i4, int i5);
     complex<double> Mpmmmm(int i1, int i2, int i3, int i4, int i5);
     complex<double> Mppmmm(int i1, int i2, int i3, int i4, int i5);
     complex<double> Mpmpmm(int i1, int i2, int i3, int i4, int i5);

     double Cmpppp(int i1, int i2, int i3, int i4, int i5);
     double Cppppm(int i1, int i2, int i3, int i4, int i5);
     double CFmmppp(int i1, int i2, int i3, int i4, int i5);
     double CBmmppp(int i1, int i2, int i3, int i4, int i5);
     double CF1pppmm(int i1, int i2, int i3, int i4, int i5);
     double CF2pppmm(int i1, int i2, int i3, int i4, int i5);
     double CB1pppmm(int i1, int i2, int i3, int i4, int i5);
     double CB2apppmm(int i1, int i2, int i3, int i4, int i5);
     double CB2bpppmm(int i1, int i2, int i3, int i4, int i5);
     double C1pppmm(int i1, int i2, int i3, int i4, int i5);

     double CF1mpppm(int i1, int i2, int i3, int i4, int i5);
     double CF2mpppm(int i1, int i2, int i3, int i4, int i5);
     double CF3mpppm(int i1, int i2, int i3, int i4, int i5);
     double CF4ampppm(int i1, int i2, int i3, int i4, int i5);
     double CF4bmpppm(int i1, int i2, int i3, int i4, int i5);
     double CB1mpppm(int i1, int i2, int i3, int i4, int i5);
     double CB2mpppm(int i1, int i2, int i3, int i4, int i5);
     double CB3mpppm(int i1, int i2, int i3, int i4, int i5);
     double CB4ampppm(int i1, int i2, int i3, int i4, int i5);
     double CB4bmpppm(int i1, int i2, int i3, int i4, int i5);
     double C1mpppm(int i1, int i2, int i3, int i4, int i5);


     virtual void initializeSubprocess(void);

 public:

     ggNLO(double energy=14000.0, collider coll=pp, 
	   double m=115.0, event_type et=NLO, int np=0) :
                     gamma2NLO(energy,coll,m,et,np) {}

     virtual ~ggNLO(void) {}

     virtual double wgtLO(void);
     virtual double wgtV(void);
     virtual double wgtVplus(void);
     virtual double wgtVa(void);
     virtual double wgtVb(void);
     virtual double wgtR(void);
     virtual double wgtSa(void);
     virtual double wgtSb(void);


};

#endif


