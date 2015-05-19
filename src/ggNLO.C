#include "ggNLO.h"

complex<double> L0(double r1, double r2);

complex<double> L1(double r1, double r2);

complex<double> L2(double r1, double r2);

complex<double> Ls1(double r1, double r2, double r3, double r4);

void ggNLO::initializeSubprocess(void)
{ 

  long int iparton=0;
  fgg=pdf(iparton,xA,muF)*pdf(iparton,xB,muF);
  fgg0=pdf(iparton,xA0,muF)*pdf(iparton,xB0,muF);

  pgg=Nc*(1.0+pow(zm,4)+pow(z,4))/z/zm;

  x1=-w1m;
  y=-w1;
  X1=log(-x1);
  Y=log(-y);

  Mpppp=1.0;

  Mmppp=-1.0;  // The minus sign here comes from spinor phases.

  Mmmpp=-1.0-(x1-y)*(X1-Y)-0.5*(x1*x1+y*y)*((X1-Y)*(X1-Y)+PI*PI);
  Mmpmp=-1.0-(x1-1.0)/y*(X1+I*PI)-0.5*(x1*x1+1.0)/y/y*(X1*X1+2.0*I*PI*X1);
  Mmppm=-1.0-(y-1.0)/x1*(Y+I*PI)-0.5*(y*y+1.0)/x1/x1*(Y*Y+2.0*I*PI*Y);

  sigLO=sigma0*2.0*(
           norm(Mpppp)+4.0*norm(Mmppp)+norm(Mmmpp)+norm(Mmpmp)+norm(Mmppm) );

}


double ggNLO::wgtLO(void)
{ 
  return sigLO*fgg0;
}

double ggNLO::wgtV(void)
{ 
  double nloV=2.0/3.0*PI*PI*Nc+beta0*logRF
              +4.0*Nc*logzmin*(logzmin+logHF);

complex<double> FLpppp,FLmppp,FLppmp,FLmmpp,FLmpmp,FLmppm;
complex<double> FSLpppp,FSLmppp,FSLppmp,FSLmmpp,FSLmpmp,FSLmppm;

  FLpppp=0.5;

// The minus sign in mppp and ppmp comes from spinor phases.

  FLmppp=-1.0/8.0*(  (2.0+4.0*x1/y/y-5.0*x1*x1/y/y)*(X1*X1+2.0*I*PI*X1)
                   -(1.0-x1*y)*((X1-Y)*(X1-Y)+PI*PI)
                   +2.0*(9.0/y-10.0*x1)*(X1+I*PI)
                   +(2.0+4.0*y/x1/x1-5.0*y*y/x1/x1)*(Y*Y+2.0*I*PI*Y)
                   -(1.0-y*x1)*((X1-Y)*(X1-Y)+PI*PI)
                   +2.0*(9.0/x1-10.0*y)*(Y+I*PI));

  FLppmp=-1.0/8.0*(  (2.0+6.0*x1/y/y-3.0*x1*x1/y/y)*(X1*X1+2.0*I*PI*X1)
                   -(x1-y)*(x1-y)*((X1-Y)*(X1-Y)+PI*PI)
                   +2.0*(9.0/y-8.0*x1)*(X1+I*PI)
                   +(2.0+6.0*y/x1/x1-3.0*y*y/x1/x1)*(Y*Y+2.0*I*PI*Y)
                   -(x1-y)*(x1-y)*((X1-Y)*(X1-Y)+PI*PI)
                   +2.0*(9.0/x1-8.0*y)*(Y+I*PI));

  FLmmpp=-(x1*x1+y*y)*(4.0*Li4(-x1)+(Y-3.0*X1)*Li3(-x1)
         +X1*X1*(Li2(-x1)-PI*PI/12.0)+1.0/48.0*pow(X1+Y,4)
                -109.0/720.0*pow(PI,4))
         -2.0*x1*y*(Li3(-x1/y)-(X1-Y)*Li2(-x1/y)-Zeta3
                +0.5*Y*((X1-Y)*(X1-Y)+PI*PI))
         -0.5*x1*x1*(Li3(-x1/y)-(X1-Y)*Li2(-x1/y)-Zeta3
                -0.5*pow(X1-Y,3)-Y*((X1-Y)*(X1-Y)+PI*PI))
         +1.0/8.0*(y*y+12.0*x1*y-27.0*x1*x1-8.0/y+9.0/y/y)*X1*X1
         -1.0/8.0*(38.0*x1*y-13.0)*X1*Y +PI*PI/48.0*(114.0*x1*y-43.0)
         -9.0/4.0*(1.0/y+2.0*x1)*X1+0.25
         +I*PI*(-(x1*x1+y*y)*(Li3(-x1/y)-(X1-Y)*Li2(-x1/y)-Zeta3
                +0.5*(Y-0.75)*((X1-Y)*(X1-Y)+PI*PI))
                +0.25*(14.0*(x1-y)-8.0/y+9.0/y/y)*X1-9.0/4.0*(1.0/y-1.0))
         -(y*y+x1*x1)*(4.0*Li4(-y)+(X1-3.0*Y)*Li3(-y)
         +Y*Y*(Li2(-y)-PI*PI/12.0)+1.0/48.0*pow(Y+X1,4)
                -109.0/720.0*pow(PI,4))
         -2.0*y*x1*(Li3(-y/x1)-(Y-X1)*Li2(-y/x1)-Zeta3
                +0.5*X1*((Y-X1)*(Y-X1)+PI*PI))
         -0.5*y*y*(Li3(-y/x1)-(Y-X1)*Li2(-y/x1)-Zeta3
                -0.5*pow(Y-X1,3)-X1*((Y-X1)*(Y-X1)+PI*PI))
         +1.0/8.0*(x1*x1+12.0*y*x1-27.0*y*y-8.0/x1+9.0/x1/x1)*Y*Y
         -1.0/8.0*(38.0*y*x1-13.0)*Y*X1 +PI*PI/48.0*(114.0*y*x1-43.0)
         -9.0/4.0*(1.0/x1+2.0*y)*Y+0.25
         +I*PI*(-(y*y+x1*x1)*(Li3(-y/x1)-(Y-X1)*Li2(-y/x1)-Zeta3
                +0.5*(X1-0.75)*((Y-X1)*(Y-X1)+PI*PI))
                +0.25*(14.0*(y-x1)-8.0/x1+9.0/x1/x1)*Y-9.0/4.0*(1.0/x1-1.0));

  FLmpmp=2.0*(x1*x1+1.0)/y/y*(Li4(-x1/y)-Li4(-y)+0.5*(X1-2.0*Y)*(Li3(-x1)-Zeta3)
         -1.0/48.0*pow(X1,4)-1.0/6.0*X1*Y*Y*Y+1.0/24.0*pow(Y,4)
         +PI*PI/24.0*(7.0*X1*X1+2.0*Y*Y)+7.0/360.0*pow(PI,4))
         +4.0*x1*(x1-3.0)/y/y*(Li4(-x1)+Li4(-x1/y)-Li4(-y)-Y*(Li3(-x1)-Zeta3)
         +PI*PI/6.0*(Li2(-x1)+0.5*Y*Y)-1.0/6.0*X1*Y*Y*Y+1.0/24.0*pow(Y,4)
         -7.0/360.0*pow(PI,4))
         +2.0/3.0*(8.0-x1+30.0*x1/y)*(Li3(-x1)-Li3(-y)-(X1+Y)*Li2(-x1)-0.5*X1*Y*Y)
         -1.0/6.0*(47.0+154.0*x1/y-4.0*x1*x1/y/y)*(Li3(-x1)-X1*Li2(-x1)-Zeta3)
         +1.0/12.0*(3.0-2.0/y-12.0*x1/y/y)*X1*(X1*X1+3.0*PI*PI)-1.0/3.0*y*X1*Y*Y
         +PI*PI/3.0*(y-x1+2*x1/y)*Y
         +5.0/9.0*PI*PI*x1*(1.0-2.0/y-10.0*x1/y/y)*X1+2.0*(1.0+2.0/y)*Zeta3
         -1.0/24.0*(15.0-14.0*x1/y-48.0*x1/y/y)*X1*X1
         +1.0/24.0*(y*y-24.0*y+44.0-8.0*x1*x1*x1/y)*((X1-Y)*(X1-Y)+PI*PI)
         +4.0/9.0*PI*PI*x1/y
         +1.0/24.0*(8.0*x1/y+60.0-24.0*y/x1+27*y*y/x1/x1)*Y*Y
         +1.0/12.0*(2*x1*x1-54*x1-27.0*y*y)*(X1/y+Y/x1)
         +I*PI*((18.0*x1/y/y+4.0*x1/y-1.0)*(Li3(-x1)-Zeta3)
         -(x1*x1+1.0)/6.0/y/y*X1*(2.0*X1*X1-PI*PI)
         -1.0/6.0*(13.0-8.0*x1+78.0*x1/y+4.0/y/y)*Li2(-x1)
         -1.0/3.0*x1*(3.0+2.0*x1/y-10.0*x1/y/y)*X1*X1
         -1.0/3.0*y*Y*(Y+2.0*X1)-PI*PI/36.0*(7.0+4.0*x1+18.0*x1/y-4.0/y/y)
         -1.0/12.0*(15.0-14.0*x1/y-48.0*x1/y/y)*X1
         +1.0/12.0*(8.0*x1/y+60.0-24.0*y/x1+27*y*y/x1/x1)*Y
         -1.0/12.0*(2.0*x1/y-54.0/y-27.0*y/x1));

  FLmppm=2.0*(y*y+1.0)/x1/x1*(Li4(-y/x1)-Li4(-x1)+0.5*(Y-2.0*X1)*(Li3(-y)-Zeta3)
         -1.0/48.0*pow(Y,4)-1.0/6.0*Y*X1*X1*X1+1.0/24.0*pow(X1,4)
         +PI*PI/24.0*(7.0*Y*Y+2.0*X1*X1)+7.0/360.0*pow(PI,4))
         +4.0*y*(y-3.0)/x1/x1*(Li4(-y)+Li4(-y/x1)-Li4(-x1)-X1*(Li3(-y)-Zeta3)
         +PI*PI/6.0*(Li2(-y)+0.5*X1*X1)-1.0/6.0*Y*X1*X1*X1+1.0/24.0*pow(X1,4)
         -7.0/360.0*pow(PI,4))
         +2.0/3.0*(8.0-y+30.0*y/x1)*(Li3(-y)-Li3(-x1)-(Y+X1)*Li2(-y)-0.5*Y*X1*X1)
         -1.0/6.0*(47.0+154.0*y/x1-4.0*y*y/x1/x1)*(Li3(-y)-Y*Li2(-y)-Zeta3)
         +1.0/12.0*(3.0-2.0/x1-12.0*y/x1/x1)*Y*(Y*Y+3.0*PI*PI)-1.0/3.0*x1*Y*X1*X1
         +PI*PI/3.0*(x1-y+2*y/x1)*X1
         +5.0/9.0*PI*PI*y*(1.0-2.0/x1-10.0*y/x1/x1)*Y+2.0*(1.0+2.0/x1)*Zeta3
         -1.0/24.0*(15.0-14.0*y/x1-48.0*y/x1/x1)*Y*Y
         +1.0/24.0*(x1*x1-24.0*x1+44.0-8.0*y*y*y/x1)*((Y-X1)*(Y-X1)+PI*PI)
         +4.0/9.0*PI*PI*y/x1
         +1.0/24.0*(8.0*y/x1+60.0-24.0*x1/y+27*x1*x1/y/y)*X1*X1
         +1.0/12.0*(2*y*y-54*y-27.0*x1*x1)*(Y/x1+X1/y)
         +I*PI*((18.0*y/x1/x1+4.0*y/x1-1.0)*(Li3(-y)-Zeta3)
         -(y*y+1.0)/6.0/x1/x1*Y*(2.0*Y*Y-PI*PI)
         -1.0/6.0*(13.0-8.0*y+78.0*y/x1+4.0/x1/x1)*Li2(-y)
         -1.0/3.0*y*(3.0+2.0*y/x1-10.0*y/x1/x1)*Y*Y
         -1.0/3.0*x1*X1*(X1+2.0*Y)-PI*PI/36.0*(7.0+4.0*y+18.0*y/x1-4.0/x1/x1)
         -1.0/12.0*(15.0-14.0*y/x1-48.0*y/x1/x1)*Y
         +1.0/12.0*(8.0*y/x1+60.0-24.0*x1/y+27*x1*x1/y/y)*X1
         -1.0/12.0*(2.0*y/x1-54.0/x1-27.0*x1/y));

  FSLpppp=-1.5;

// The minus sign in mppp and ppmp comes from spinor phases.

  FSLmppp=-1.0/8.0*( (x1*x1+1)/y/y*(X1*X1+2.0*I*PI*X1)
                   +0.5*(x1*x1+y*y)*((X1-Y)*(X1-Y)+PI*PI)
                   -4.0*(1.0/y-x1)*(X1+I*PI)
                   +(y*y+1)/x1/x1*(Y*Y+2.0*I*PI*Y)
                   +0.5*(x1*x1+y*y)*((X1-Y)*(X1-Y)+PI*PI)
                   -4.0*(1.0/x1-y)*(Y+I*PI));


  FSLppmp=FSLmppp;

  FSLmmpp=-2.0*(x1*x1+y*y)*Li4(-x1)-(x1-y)*Li4(-x1/y)+2.0*x1*x1*X1*(Li3(-x1)+Li3(-y))
         +PI*PI/6.0*(x1-y)*Li2(-x1)-x1*x1*(1.0/6.0*pow(X1,4)-2.0/3.0*X1*X1*X1*Y
            +PI*PI/6.0*X1*Y-4.0/45.0*pow(PI,4))
         -x1*(2.0*Li3(-x1)-2.0*X1*Li2(-x1)-Li3(-x1/y)+(X1-Y)*Li2(-x1/y)-3.0*Zeta3
         -2.0/3.0*X1*(X1*X1+PI*PI)+5.0/12.0*(X1-Y)*((X1-Y)*(X1-Y)+PI*PI)+X1*X1*(X1-Y))
         -0.25*(2.0-3.0*x1*x1+pow(x1,4)/y/y)*X1*X1+0.25*(2.0*x1*y+3.0)*X1*Y
         -1.0/24.0*(6.0*x1*y+7.0)*PI*PI+(0.5/y+x1)*X1-0.25
         +I*PI*(-2.0*x1*x1*(Li3(-y/x1)-Zeta3)+2.0*x1*(Li2(-y/x1)-0.25*(X1-Y)*(X1-Y)
         -PI*PI/12.0)-0.5*(2.0*x1*x1-1.0)/y/y*X1+0.5/y-0.5)
         -2.0*(y*y+x1*x1)*Li4(-y)-(y-x1)*Li4(-y/x1)+2.0*y*y*Y*(Li3(-y)+Li3(-x1))
         +PI*PI/6.0*(y-x1)*Li2(-y)-y*y*(1.0/6.0*pow(Y,4)-2.0/3.0*Y*Y*Y*X1
            +PI*PI/6.0*Y*X1-4.0/45.0*pow(PI,4))
         -y*(2.0*Li3(-y)-2.0*Y*Li2(-y)-Li3(-y/x1)+(Y-X1)*Li2(-y/x1)-3.0*Zeta3
         -2.0/3.0*Y*(Y*Y+PI*PI)+5.0/12.0*(Y-X1)*((Y-X1)*(Y-X1)+PI*PI)+Y*Y*(Y-X1))
         -0.25*(2.0-3.0*y*y+pow(y,4)/x1/x1)*Y*Y+0.25*(2.0*y*x1+3.0)*Y*X1
         -1.0/24.0*(6.0*y*x1+7.0)*PI*PI+(0.5/x1+y)*Y-0.25
         +I*PI*(-2.0*y*y*(Li3(-x1/y)-Zeta3)+2.0*y*(Li2(-x1/y)-0.25*(Y-X1)*(Y-X1)
         -PI*PI/12.0)-0.5*(2.0*y*y-1.0)/x1/x1*Y+0.5/x1-0.5);

  FSLmpmp=-2.0*(x1*x1+1.0)/y/y*(Li4(-x1/y)-Li4(-y)-Y*(Li3(-x1)-Zeta3)
          +0.5*X1*(Li3(-x1)-Zeta3)+1.0/24.0*pow(X1,4)-1.0/6.0*X1*Y*Y*Y
          +1.0/24.0*pow(Y,4)+PI*PI/12.0*Y*Y+7.0/360.0*pow(PI,4))
          -2.0*(x1-1.0)/y*(Li4(-x1)-0.5*X1*(Li3(-x1)-Zeta3)
          +PI*PI/6.0*(Li2(-x1)-0.5*X1*X1)
          -1.0/48.0*pow(X1,4)-7.0/180.0*pow(PI,4))
          +(2.0*x1/y-1.0)*(Li3(-x1)-X1*Li2(-x1)+Zeta3-1.0/6.0*X1*X1*X1
          -PI*PI/3.0*(X1+Y))+2.0*(2.0*x1/y+1.0)*(Li3(-y)+Y*Li2(-x1)-Zeta3
          +0.25*X1*(2.0*Y*Y+PI*PI)-1.0/8.0*X1*X1*X1)-0.25*(2.0*x1*x1-y*y)*
          ((X1-Y)*(X1-Y)+PI*PI)-0.25*(3.0+2.0*x1/y/y)*X1*X1
          -(2.0-y*y)/4.0/x1/x1*Y*Y+PI*PI/6.0+0.5*(2.0*x1+y*y)*(X1/y+Y/x1)-0.5
          +I*PI*(2.0/y/y*(Li3(-x1)-Zeta3)-(x1*x1+1.0)/6.0/y/y*X1*X1*X1
          -(2.0/y-1.0)*Li2(-x1)-3.0/4.0*(x1-1.0)/y*X1*X1
          -(x1/y/y+1.5)*X1-(2.0-y*y)/2.0/x1/x1*Y-1.0/y-y/2.0/x1);




  FSLmppm=-2.0*(y*y+1.0)/x1/x1*(Li4(-y/x1)-Li4(-x1)-X1*(Li3(-y)-Zeta3)
          +0.5*Y*(Li3(-y)-Zeta3)+1.0/24.0*pow(Y,4)-1.0/6.0*Y*X1*X1*X1
          +1.0/24.0*pow(X1,4)+PI*PI/12.0*X1*X1+7.0/360.0*pow(PI,4))
          -2.0*(y-1.0)/x1*(Li4(-y)-0.5*Y*(Li3(-y)-Zeta3)
          +PI*PI/6.0*(Li2(-y)-0.5*Y*Y)
          -1.0/48.0*pow(Y,4)-7.0/180.0*pow(PI,4))
          +(2.0*y/x1-1.0)*(Li3(-y)-Y*Li2(-y)+Zeta3-1.0/6.0*Y*Y*Y
          -PI*PI/3.0*(Y+X1))+2.0*(2.0*y/x1+1.0)*(Li3(-x1)+X1*Li2(-y)-Zeta3
          +0.25*Y*(2.0*X1*X1+PI*PI)-1.0/8.0*Y*Y*Y)-0.25*(2.0*y*y-x1*x1)*
          ((Y-X1)*(Y-X1)+PI*PI)-0.25*(3.0+2.0*y/x1/x1)*Y*Y
          -(2.0-x1*x1)/4.0/y/y*X1*X1+PI*PI/6.0+0.5*(2.0*y+x1*x1)*(Y/x1+X1/y)-0.5
          +I*PI*(2.0/x1/x1*(Li3(-y)-Zeta3)-(y*y+1.0)/6.0/x1/x1*Y*Y*Y
          -(2.0/x1-1.0)*Li2(-y)-3.0/4.0*(y-1.0)/x1*Y*Y
          -(y/x1/x1+1.5)*Y-(2.0-x1*x1)/2.0/y/y*X1-1.0/x1-x1/2.0/y);


  return as2pi*fgg0*(sigLO*nloV+sigma0*4.0*real(

  // factor 2 from 2*Re and factor of 2 from + <-> - everywhere:

                 conj(Mpppp)*(Nc*FLpppp-FSLpppp/Nc)
            +2.0*conj(Mmppp)*(Nc*(FLmppp+FLppmp)-(FSLmppp+FSLppmp)/Nc)
                +conj(Mmmpp)*(Nc*FLmmpp-FSLmmpp/Nc)
                +conj(Mmpmp)*(Nc*FLmpmp-FSLmpmp/Nc)
                +conj(Mmppm)*(Nc*FLmppm-FSLmppm/Nc) ));


  // FSL terms initially had wrong sign!!!  Now Correct.
  // Corrected on 7/30/02.

}

double ggNLO::wgtR(void)
{ 

complex<double> M;
double Msq=0.0;

//  (1+2+3+4+5+)   (Checked out versus Zvi.)

M = Mppppp(1,2,3,4,5)+Mppppp(1,4,5,2,3)+Mppppp(1,2,4,5,3)
   +Mppppp(1,2,3,5,4)+Mppppp(1,5,4,2,3)+Mppppp(1,2,5,4,3)
   +Mppppp(1,2,4,3,5)+Mppppp(1,5,2,3,4)+Mppppp(1,4,2,5,3)
   +Mppppp(1,2,5,3,4)+Mppppp(1,4,2,3,5)+Mppppp(1,5,2,4,3);

Msq+=norm(M);

//  (1-2+3+4+5+)   (Checked out versus Zvi.)

M = Mmpppp(1,2,3,4,5)+Mmpppp(1,4,5,2,3)+Mmpppp(1,2,4,5,3)
   +Mmpppp(1,2,3,5,4)+Mmpppp(1,5,4,2,3)+Mmpppp(1,2,5,4,3)
   +Mmpppp(1,2,4,3,5)+Mmpppp(1,5,2,3,4)+Mmpppp(1,4,2,5,3)
   +Mmpppp(1,2,5,3,4)+Mmpppp(1,4,2,3,5)+Mmpppp(1,5,2,4,3);

Msq+=norm(M);

//  (1+2-3+4+5+)   (Checked out versus Zvi.)

M = Mmpppp(2,3,4,5,1)+Mmpppp(2,3,1,4,5)+Mmpppp(2,4,5,3,1)
   +Mmpppp(2,3,5,4,1)+Mmpppp(2,3,1,5,4)+Mmpppp(2,5,4,3,1)
   +Mmpppp(2,4,3,5,1)+Mmpppp(2,3,4,1,5)+Mmpppp(2,5,3,1,4)
   +Mmpppp(2,5,3,4,1)+Mmpppp(2,3,5,1,4)+Mmpppp(2,4,3,1,5);

Msq+=norm(M);

//  (1+2+3-4+5+)   (Checked out versus Zvi.)

M = Mmpppp(3,4,5,1,2)+Mmpppp(3,1,4,5,2)+Mmpppp(3,1,2,4,5)
   +Mmpppp(3,5,4,1,2)+Mmpppp(3,1,5,4,2)+Mmpppp(3,1,2,5,4)
   +Mmpppp(3,5,1,2,4)+Mmpppp(3,4,1,5,2)+Mmpppp(3,1,4,2,5)
   +Mmpppp(3,4,1,2,5)+Mmpppp(3,5,1,4,2)+Mmpppp(3,1,5,2,4);

Msq+=norm(M);

//  (1+2+3+4-5+)   (Checked out versus Zvi.)

M = Mmpppp(4,5,1,2,3)+Mmpppp(4,5,2,3,1)+Mmpppp(4,5,3,1,2)
   +Mmpppp(4,1,2,3,5)+Mmpppp(4,2,3,1,5)+Mmpppp(4,3,1,2,5)
   +Mmpppp(4,3,5,1,2)+Mmpppp(4,1,5,2,3)+Mmpppp(4,2,5,3,1)
   +Mmpppp(4,1,2,5,3)+Mmpppp(4,2,3,5,1)+Mmpppp(4,3,1,5,2);

Msq+=norm(M);

//  (1+2+3+4+5-)   (Checked out versus Zvi.)

M = Mmpppp(5,1,2,3,4)+Mmpppp(5,2,3,1,4)+Mmpppp(5,3,1,2,4)
   +Mmpppp(5,4,1,2,3)+Mmpppp(5,4,2,3,1)+Mmpppp(5,4,3,1,2)
   +Mmpppp(5,1,2,4,3)+Mmpppp(5,2,3,4,1)+Mmpppp(5,3,1,4,2)
   +Mmpppp(5,3,4,1,2)+Mmpppp(5,1,4,2,3)+Mmpppp(5,2,4,3,1);

Msq+=norm(M);

//  (1-2-3+4+5+)   (Checked out versus Zvi.)

M = Mmmppp(1,2,3,4,5)+Mmpmpp(2,3,1,4,5)+Mmmppp(1,2,4,5,3)
   +Mmmppp(1,2,3,5,4)+Mmpmpp(2,3,1,5,4)+Mmmppp(1,2,5,4,3)
   +Mmmppp(1,2,4,3,5)+Mmpmpp(1,5,2,3,4)+Mmpmpp(1,4,2,5,3)
   +Mmmppp(1,2,5,3,4)+Mmpmpp(1,4,2,3,5)+Mmpmpp(1,5,2,4,3);

Msq+=norm(M);

//  (1-2+3-4+5+)   (Checked out versus Zvi.)

M = Mmpmpp(1,2,3,4,5)+Mmmppp(3,1,4,5,2)+Mmmppp(3,1,2,4,5)
   +Mmpmpp(1,2,3,5,4)+Mmmppp(3,1,5,4,2)+Mmmppp(3,1,2,5,4)
   +Mmpmpp(3,5,1,2,4)+Mmpmpp(3,4,1,5,2)+Mmmppp(3,1,4,2,5)
   +Mmpmpp(3,4,1,2,5)+Mmpmpp(3,5,1,4,2)+Mmmppp(3,1,5,2,4);

Msq+=norm(M);

//  (1-2+3+4-5+)   (Checked out versus Zvi.)

M = Mmpmpp(4,5,1,2,3)+Mmmppp(1,4,5,2,3)+Mmpmpp(1,2,4,5,3)
   +Mmmppp(4,1,2,3,5)+Mmpmpp(1,5,4,2,3)+Mmpmpp(4,3,1,2,5)
   +Mmpmpp(1,2,4,3,5)+Mmmppp(4,1,5,2,3)+Mmmppp(1,4,2,5,3)
   +Mmmppp(4,1,2,5,3)+Mmmppp(1,4,2,3,5)+Mmpmpp(4,3,1,5,2);

Msq+=norm(M);

//  (1-2+3+4+5-)   (Checked out versus Zvi.)

M = Mmmppp(5,1,2,3,4)+Mmpmpp(1,4,5,2,3)+Mmpmpp(5,3,1,2,4)
   +Mmpmpp(5,4,1,2,3)+Mmmppp(1,5,4,2,3)+Mmpmpp(1,2,5,4,3)
   +Mmmppp(5,1,2,4,3)+Mmmppp(1,5,2,3,4)+Mmpmpp(5,3,1,4,2)
   +Mmpmpp(1,2,5,3,4)+Mmmppp(5,1,4,2,3)+Mmmppp(1,5,2,4,3);

Msq+=norm(M);

//  (1+2-3-4+5+)   (Checked out versus Zvi.)

M = Mmmppp(2,3,4,5,1)+Mmmppp(2,3,1,4,5)+Mmpmpp(3,1,2,4,5)
   +Mmmppp(2,3,5,4,1)+Mmmppp(2,3,1,5,4)+Mmpmpp(3,1,2,5,4)
   +Mmpmpp(2,4,3,5,1)+Mmmppp(2,3,4,1,5)+Mmpmpp(2,5,3,1,4)
   +Mmpmpp(2,5,3,4,1)+Mmmppp(2,3,5,1,4)+Mmpmpp(2,4,3,1,5);

Msq+=norm(M);

//  (1+2-3+4-5+)   (Checked out versus Zvi.)

M = Mmpmpp(2,3,4,5,1)+Mmpmpp(4,5,2,3,1)+Mmmppp(2,4,5,3,1)
   +Mmpmpp(4,1,2,3,5)+Mmmppp(4,2,3,1,5)+Mmpmpp(2,5,4,3,1)
   +Mmmppp(2,4,3,5,1)+Mmpmpp(2,3,4,1,5)+Mmmppp(4,2,5,3,1)
   +Mmpmpp(4,1,2,5,3)+Mmmppp(4,2,3,5,1)+Mmmppp(2,4,3,1,5);

Msq+=norm(M);

//  (1+2-3+4+5-)   (Checked out versus Zvi.)

M = Mmpmpp(5,1,2,3,4)+Mmmppp(5,2,3,1,4)+Mmpmpp(2,4,5,3,1)
   +Mmpmpp(2,3,5,4,1)+Mmpmpp(5,4,2,3,1)+Mmmppp(2,5,4,3,1)
   +Mmpmpp(5,1,2,4,3)+Mmmppp(5,2,3,4,1)+Mmmppp(2,5,3,1,4)
   +Mmmppp(2,5,3,4,1)+Mmpmpp(2,3,5,1,4)+Mmmppp(5,2,4,3,1);

Msq+=norm(M);

//  (1+2+3-4-5+)   (Checked out versus Zvi.)

M = Mmmppp(3,4,5,1,2)+Mmpmpp(3,1,4,5,2)+Mmpmpp(4,5,3,1,2)
   +Mmpmpp(3,5,4,1,2)+Mmpmpp(4,2,3,1,5)+Mmmppp(4,3,1,2,5)
   +Mmmppp(4,3,5,1,2)+Mmmppp(3,4,1,5,2)+Mmpmpp(3,1,4,2,5)
   +Mmmppp(3,4,1,2,5)+Mmpmpp(4,2,3,5,1)+Mmmppp(4,3,1,5,2);

Msq+=norm(M);

//  (1+2+3-4+5-)   (Checked out versus Zvi.)

M = Mmpmpp(3,4,5,1,2)+Mmpmpp(5,2,3,1,4)+Mmmppp(5,3,1,2,4)
   +Mmmppp(3,5,4,1,2)+Mmpmpp(3,1,5,4,2)+Mmpmpp(5,4,3,1,2)
   +Mmmppp(3,5,1,2,4)+Mmpmpp(5,2,3,4,1)+Mmmppp(5,3,1,4,2)
   +Mmmppp(5,3,4,1,2)+Mmmppp(3,5,1,4,2)+Mmpmpp(3,1,5,2,4);

Msq+=norm(M);

//  (1+2+3+4-5-)   (Checked out versus Zvi.)

M = Mmmppp(4,5,1,2,3)+Mmmppp(4,5,2,3,1)+Mmmppp(4,5,3,1,2)
   +Mmmppp(5,4,1,2,3)+Mmmppp(5,4,2,3,1)+Mmmppp(5,4,3,1,2)
   +Mmpmpp(4,3,5,1,2)+Mmpmpp(4,1,5,2,3)+Mmpmpp(4,2,5,3,1)
   +Mmpmpp(5,3,4,1,2)+Mmpmpp(5,1,4,2,3)+Mmpmpp(5,2,4,3,1);

Msq+=norm(M);

//  (1-2-3-4-5-)   (Checked out versus Zvi.)

M = Mmmmmm(1,2,3,4,5)+Mmmmmm(1,4,5,2,3)+Mmmmmm(1,2,4,5,3)
   +Mmmmmm(1,2,3,5,4)+Mmmmmm(1,5,4,2,3)+Mmmmmm(1,2,5,4,3)
   +Mmmmmm(1,2,4,3,5)+Mmmmmm(1,5,2,3,4)+Mmmmmm(1,4,2,5,3)
   +Mmmmmm(1,2,5,3,4)+Mmmmmm(1,4,2,3,5)+Mmmmmm(1,5,2,4,3);

Msq+=norm(M);

//  (1+2-3-4-5-)   (Checked out versus Zvi.)

M = Mpmmmm(1,2,3,4,5)+Mpmmmm(1,4,5,2,3)+Mpmmmm(1,2,4,5,3)
   +Mpmmmm(1,2,3,5,4)+Mpmmmm(1,5,4,2,3)+Mpmmmm(1,2,5,4,3)
   +Mpmmmm(1,2,4,3,5)+Mpmmmm(1,5,2,3,4)+Mpmmmm(1,4,2,5,3)
   +Mpmmmm(1,2,5,3,4)+Mpmmmm(1,4,2,3,5)+Mpmmmm(1,5,2,4,3);

Msq+=norm(M);

//  (1-2+3-4-5-)   (Checked out versus Zvi.)

M = Mpmmmm(2,3,4,5,1)+Mpmmmm(2,3,1,4,5)+Mpmmmm(2,4,5,3,1)
   +Mpmmmm(2,3,5,4,1)+Mpmmmm(2,3,1,5,4)+Mpmmmm(2,5,4,3,1)
   +Mpmmmm(2,4,3,5,1)+Mpmmmm(2,3,4,1,5)+Mpmmmm(2,5,3,1,4)
   +Mpmmmm(2,5,3,4,1)+Mpmmmm(2,3,5,1,4)+Mpmmmm(2,4,3,1,5);

Msq+=norm(M);

//  (1-2-3+4-5-)   (Checked out versus Zvi.)

M = Mpmmmm(3,4,5,1,2)+Mpmmmm(3,1,4,5,2)+Mpmmmm(3,1,2,4,5)
   +Mpmmmm(3,5,4,1,2)+Mpmmmm(3,1,5,4,2)+Mpmmmm(3,1,2,5,4)
   +Mpmmmm(3,5,1,2,4)+Mpmmmm(3,4,1,5,2)+Mpmmmm(3,1,4,2,5)
   +Mpmmmm(3,4,1,2,5)+Mpmmmm(3,5,1,4,2)+Mpmmmm(3,1,5,2,4);

Msq+=norm(M);

//  (1-2-3-4+5-)   (Checked out versus Zvi.)

M = Mpmmmm(4,5,1,2,3)+Mpmmmm(4,5,2,3,1)+Mpmmmm(4,5,3,1,2)
   +Mpmmmm(4,1,2,3,5)+Mpmmmm(4,2,3,1,5)+Mpmmmm(4,3,1,2,5)
   +Mpmmmm(4,3,5,1,2)+Mpmmmm(4,1,5,2,3)+Mpmmmm(4,2,5,3,1)
   +Mpmmmm(4,1,2,5,3)+Mpmmmm(4,2,3,5,1)+Mpmmmm(4,3,1,5,2);

Msq+=norm(M);

//  (1-2-3-4-5+)   (Checked out versus Zvi.)

M = Mpmmmm(5,1,2,3,4)+Mpmmmm(5,2,3,1,4)+Mpmmmm(5,3,1,2,4)
   +Mpmmmm(5,4,1,2,3)+Mpmmmm(5,4,2,3,1)+Mpmmmm(5,4,3,1,2)
   +Mpmmmm(5,1,2,4,3)+Mpmmmm(5,2,3,4,1)+Mpmmmm(5,3,1,4,2)
   +Mpmmmm(5,3,4,1,2)+Mpmmmm(5,1,4,2,3)+Mpmmmm(5,2,4,3,1);

Msq+=norm(M);

//  (1+2+3-4-5-)   (Checked out versus Zvi.)

M = Mppmmm(1,2,3,4,5)+Mpmpmm(2,3,1,4,5)+Mppmmm(1,2,4,5,3)
   +Mppmmm(1,2,3,5,4)+Mpmpmm(2,3,1,5,4)+Mppmmm(1,2,5,4,3)
   +Mppmmm(1,2,4,3,5)+Mpmpmm(1,5,2,3,4)+Mpmpmm(1,4,2,5,3)
   +Mppmmm(1,2,5,3,4)+Mpmpmm(1,4,2,3,5)+Mpmpmm(1,5,2,4,3);

Msq+=norm(M);

//  (1+2-3+4-5-)   (Checked out versus Zvi.)

M = Mpmpmm(1,2,3,4,5)+Mppmmm(3,1,4,5,2)+Mppmmm(3,1,2,4,5)
   +Mpmpmm(1,2,3,5,4)+Mppmmm(3,1,5,4,2)+Mppmmm(3,1,2,5,4)
   +Mpmpmm(3,5,1,2,4)+Mpmpmm(3,4,1,5,2)+Mppmmm(3,1,4,2,5)
   +Mpmpmm(3,4,1,2,5)+Mpmpmm(3,5,1,4,2)+Mppmmm(3,1,5,2,4);

Msq+=norm(M);

//  (1+2-3-4+5-)   (Checked out versus Zvi.)

M = Mpmpmm(4,5,1,2,3)+Mppmmm(1,4,5,2,3)+Mpmpmm(1,2,4,5,3)
   +Mppmmm(4,1,2,3,5)+Mpmpmm(1,5,4,2,3)+Mpmpmm(4,3,1,2,5)
   +Mpmpmm(1,2,4,3,5)+Mppmmm(4,1,5,2,3)+Mppmmm(1,4,2,5,3)
   +Mppmmm(4,1,2,5,3)+Mppmmm(1,4,2,3,5)+Mpmpmm(4,3,1,5,2);

Msq+=norm(M);

//  (1+2-3-4-5+)   (Checked out versus Zvi.)

M = Mppmmm(5,1,2,3,4)+Mpmpmm(1,4,5,2,3)+Mpmpmm(5,3,1,2,4)
   +Mpmpmm(5,4,1,2,3)+Mppmmm(1,5,4,2,3)+Mpmpmm(1,2,5,4,3)
   +Mppmmm(5,1,2,4,3)+Mppmmm(1,5,2,3,4)+Mpmpmm(5,3,1,4,2)
   +Mpmpmm(1,2,5,3,4)+Mppmmm(5,1,4,2,3)+Mppmmm(1,5,2,4,3);

Msq+=norm(M);

//  (1-2+3+4-5-)   (Checked out versus Zvi.)

M = Mppmmm(2,3,4,5,1)+Mppmmm(2,3,1,4,5)+Mpmpmm(3,1,2,4,5)
   +Mppmmm(2,3,5,4,1)+Mppmmm(2,3,1,5,4)+Mpmpmm(3,1,2,5,4)
   +Mpmpmm(2,4,3,5,1)+Mppmmm(2,3,4,1,5)+Mpmpmm(2,5,3,1,4)
   +Mpmpmm(2,5,3,4,1)+Mppmmm(2,3,5,1,4)+Mpmpmm(2,4,3,1,5);

Msq+=norm(M);

//  (1-2+3-4+5-)   (Checked out versus Zvi.)

M = Mpmpmm(2,3,4,5,1)+Mpmpmm(4,5,2,3,1)+Mppmmm(2,4,5,3,1)
   +Mpmpmm(4,1,2,3,5)+Mppmmm(4,2,3,1,5)+Mpmpmm(2,5,4,3,1)
   +Mppmmm(2,4,3,5,1)+Mpmpmm(2,3,4,1,5)+Mppmmm(4,2,5,3,1)
   +Mpmpmm(4,1,2,5,3)+Mppmmm(4,2,3,5,1)+Mppmmm(2,4,3,1,5);

Msq+=norm(M);

//  (1-2+3-4-5+)   (Checked out versus Zvi.)

M = Mpmpmm(5,1,2,3,4)+Mppmmm(5,2,3,1,4)+Mpmpmm(2,4,5,3,1)
   +Mpmpmm(2,3,5,4,1)+Mpmpmm(5,4,2,3,1)+Mppmmm(2,5,4,3,1)
   +Mpmpmm(5,1,2,4,3)+Mppmmm(5,2,3,4,1)+Mppmmm(2,5,3,1,4)
   +Mppmmm(2,5,3,4,1)+Mpmpmm(2,3,5,1,4)+Mppmmm(5,2,4,3,1);

Msq+=norm(M);

//  (1-2-3+4+5-)   (Checked out versus Zvi.)

M = Mppmmm(3,4,5,1,2)+Mpmpmm(3,1,4,5,2)+Mpmpmm(4,5,3,1,2)
   +Mpmpmm(3,5,4,1,2)+Mpmpmm(4,2,3,1,5)+Mppmmm(4,3,1,2,5)
   +Mppmmm(4,3,5,1,2)+Mppmmm(3,4,1,5,2)+Mpmpmm(3,1,4,2,5)
   +Mppmmm(3,4,1,2,5)+Mpmpmm(4,2,3,5,1)+Mppmmm(4,3,1,5,2);

Msq+=norm(M);

//  (1-2-3+4-5+)   (Checked out versus Zvi.)

M = Mpmpmm(3,4,5,1,2)+Mpmpmm(5,2,3,1,4)+Mppmmm(5,3,1,2,4)
   +Mppmmm(3,5,4,1,2)+Mpmpmm(3,1,5,4,2)+Mpmpmm(5,4,3,1,2)
   +Mppmmm(3,5,1,2,4)+Mpmpmm(5,2,3,4,1)+Mppmmm(5,3,1,4,2)
   +Mppmmm(5,3,4,1,2)+Mppmmm(3,5,1,4,2)+Mpmpmm(3,1,5,2,4);

Msq+=norm(M);

//  (1-2-3-4+5+)   (Checked out versus Zvi.)

M = Mppmmm(4,5,1,2,3)+Mppmmm(4,5,2,3,1)+Mppmmm(4,5,3,1,2)
   +Mppmmm(5,4,1,2,3)+Mppmmm(5,4,2,3,1)+Mppmmm(5,4,3,1,2)
   +Mpmpmm(4,3,5,1,2)+Mpmpmm(4,1,5,2,3)+Mpmpmm(4,2,5,3,1)
   +Mpmpmm(5,3,4,1,2)+Mpmpmm(5,1,4,2,3)+Mpmpmm(5,2,4,3,1);

Msq+=norm(M);

/// Matrix elements are multiplied by (4*pi)^2*sqrt(shat).

  double pcrossGG=Nc*zm*z*Msq/4.0;

  return sigma0*as2pi*pcrossGG*fgg;
}

double ggNLO::wgtVplus(void)
{ 
  double pcrossGG=4.0*Nc*(2.0*logzm+logHF)/zm;
  return -sigLO*as2pi*z*pcrossGG*fgg0;
}

double ggNLO::wgtVa(void)
{ 
  double L=2.0*logzm+logHF-logz;
  double pcrossGG=pgg*L;
  return sigLO*as2pi*pcrossGG*fgg;
}

double ggNLO::wgtVb(void)
{ 
  double L=2.0*logzm+logHF-logz;
  double pcrossGG=pgg*L;
  return sigLO*as2pi*pcrossGG*fgg;
}


double ggNLO::wgtSa(void)
{ 

double Dpmp=z*z;
double Dmpm=Dpmp;
double Dppp=0.0;
double Dmmm=Dppp;
double Dpmm=1.0;
double Dmpp=Dpmm;
complex<double> Dppm=zm*zm*exp(-2.0*I*(phi-phi1));
complex<double> Dmmp=conj(Dppm);

complex<double> M;
double Msq=0.0;

//  (1+2+3+4+5+)

M = Dppp*Mmppp+Dpmp*Mpppp;

Msq+=norm(M);

//  (1-2+3+4+5+)

M = Dmpp*Mmppp+Dmmp*Mpppp;

Msq+=norm(M);

//  (1+2-3+4+5+)

M = Dppp*Mmmpp+Dpmp*Mmppp;

Msq+=norm(M);

//  (1+2+3-4+5+)

M = Dppm*Mmppp+Dpmm*Mpppp;

Msq+=norm(M);

//  (1+2+3+4-5+)

M = Dppp*Mmpmp+Dpmp*Mmppp;

Msq+=norm(M);

//  (1+2+3+4+5-)

M = Dppp*Mmppm+Dpmp*Mmppp;

Msq+=norm(M);

//  (1-2-3+4+5+)

M = Dmpp*Mmmpp+Dmmp*Mmppp;

Msq+=norm(M);

//  (1-2+3-4+5+)

M = Dmpm*Mmppp+Dmmm*Mpppp;

Msq+=norm(M);

//  (1-2+3+4-5+)

M = Dmpp*Mmpmp+Dmmp*Mmppp;

Msq+=norm(M);

//  (1-2+3+4+5-)

M = Dmpp*Mmppm+Dmmp*Mmppp;

Msq+=norm(M);

//  (1+2-3-4+5+)

M = Dppm*Mmmpp+Dpmm*Mmppp;

Msq+=norm(M);

//  (1+2-3+4-5+)

M = Dppp*Mmppp+Dpmp*Mmppm;

Msq+=norm(M);

//  (1+2-3+4+5-)

M = Dppp*Mmppp+Dpmp*Mmpmp;

Msq+=norm(M);

//  (1+2+3-4-5+)

M = Dppm*Mmpmp+Dpmm*Mmppp;

Msq+=norm(M);

//  (1+2+3-4+5-)

M = Dppm*Mmppm+Dpmm*Mmppp;

Msq+=norm(M);

//  (1+2+3+4-5-)

M = Dppp*Mmppp+Dpmp*Mmmpp;

Msq+=norm(M);

//  (1-2-3-4-5-)

M = Dmmm*Mmppp+Dmpm*Mpppp;

Msq+=norm(M);

//  (1+2-3-4-5-)

M = Dpmm*Mmppp+Dppm*Mpppp;

Msq+=norm(M);

//  (1-2+3-4-5-)

M = Dmmm*Mmmpp+Dmpm*Mmppp;

Msq+=norm(M);

//  (1-2-3+4-5-)

M = Dmmp*Mmppp+Dmpp*Mpppp;

Msq+=norm(M);

//  (1-2-3-4+5-)

M = Dmmm*Mmpmp+Dmpm*Mmppp;

Msq+=norm(M);

//  (1-2-3-4-5+)

M = Dmmm*Mmppm+Dmpm*Mmppp;

Msq+=norm(M);

//  (1+2+3-4-5-)

M = Dpmm*Mmmpp+Dppm*Mmppp;

Msq+=norm(M);

//  (1+2-3+4-5-)

M = Dpmp*Mmppp+Dppp*Mpppp;

Msq+=norm(M);

//  (1+2-3-4+5-)

M = Dpmm*Mmpmp+Dppm*Mmppp;

Msq+=norm(M);

//  (1+2-3-4-5+)

M = Dpmm*Mmppm+Dppm*Mmppp;

Msq+=norm(M);

//  (1-2+3+4-5-)

M = Dmmp*Mmmpp+Dmpp*Mmppp;

Msq+=norm(M);

//  (1-2+3-4+5-)

M = Dmmm*Mmppp+Dmpm*Mmppm;

Msq+=norm(M);

//  (1-2+3-4-5+)

M = Dmmm*Mmppp+Dmpm*Mmpmp;

Msq+=norm(M);

//  (1-2-3+4+5-)

M = Dmmp*Mmpmp+Dmpp*Mmppp;

Msq+=norm(M);

//  (1-2-3+4-5+)

M = Dmmp*Mmppm+Dmpp*Mmppp;

Msq+=norm(M);

//  (1-2-3-4+5+)

M = Dmmm*Mmppp+Dmpm*Mmmpp;

Msq+=norm(M);

  double pcrossGG=Nc/w/z/zm*Msq;
  return -sigma0*as2pi*pcrossGG*fgg;
}

double ggNLO::wgtSb(void)
{ 

double Dpmp=z*z;
double Dmpm=Dpmp;
double Dppp=0.0;
double Dmmm=Dppp;
double Dpmm=1.0;
double Dmpp=Dpmm;
complex<double> Dppm=zm*zm*exp(2.0*I*(phi-phi1));
complex<double> Dmmp=conj(Dppm);

complex<double> M;
double Msq=0.0;

//  (1+2+3+4+5+)

M = Dppp*Mmppp+Dpmp*Mpppp;

Msq+=norm(M);

//  (1-2+3+4+5+)

M = Dppp*Mmmpp+Dpmp*Mmppp;

Msq+=norm(M);

//  (1+2-3+4+5+)

M = Dmpp*Mmppp+Dmmp*Mpppp;

Msq+=norm(M);

//  (1+2+3-4+5+)

M = Dppm*Mmppp+Dpmm*Mpppp;

Msq+=norm(M);

//  (1+2+3+4-5+)

M = Dppp*Mmppm+Dpmp*Mmppp;

Msq+=norm(M);

//  (1+2+3+4+5-)

M = Dppp*Mmpmp+Dpmp*Mmppp;

Msq+=norm(M);

//  (1-2-3+4+5+)

M = Dmpp*Mmmpp+Dmmp*Mmppp;

Msq+=norm(M);

//  (1-2+3-4+5+)

M = Dppm*Mmmpp+Dpmm*Mmppp;

Msq+=norm(M);

//  (1-2+3+4-5+)

M = Dppp*Mmppp+Dpmp*Mmpmp;

Msq+=norm(M);

//  (1-2+3+4+5-)

M = Dppp*Mmppp+Dpmp*Mmppm;

Msq+=norm(M);

//  (1+2-3-4+5+)

M = Dmpm*Mmppp+Dmmm*Mpppp;

Msq+=norm(M);

//  (1+2-3+4-5+)

M = Dmpp*Mmppm+Dmmp*Mmppp;

Msq+=norm(M);

//  (1+2-3+4+5-)

M = Dmpp*Mmpmp+Dmmp*Mmppp;

Msq+=norm(M);

//  (1+2+3-4-5+)

M = Dppm*Mmppm+Dpmm*Mmppp;

Msq+=norm(M);

//  (1+2+3-4+5-)

M = Dppm*Mmpmp+Dpmm*Mmppp;

Msq+=norm(M);

//  (1+2+3+4-5-)

M = Dppp*Mmppp+Dpmp*Mmmpp;

Msq+=norm(M);

//  (1-2-3-4-5-)

M = Dmmm*Mmppp+Dmpm*Mpppp;

Msq+=norm(M);

//  (1+2-3-4-5-)

M = Dmmm*Mmmpp+Dmpm*Mmppp;

Msq+=norm(M);

//  (1-2+3-4-5-)

M = Dpmm*Mmppp+Dppm*Mpppp;

Msq+=norm(M);

//  (1-2-3+4-5-)

M = Dmmp*Mmppp+Dmpp*Mpppp;

Msq+=norm(M);

//  (1-2-3-4+5-)

M = Dmmm*Mmppm+Dmpm*Mmppp;

Msq+=norm(M);

//  (1-2-3-4-5+)

M = Dmmm*Mmpmp+Dmpm*Mmppp;

Msq+=norm(M);

//  (1+2+3-4-5-)

M = Dpmm*Mmmpp+Dppm*Mmppp;

Msq+=norm(M);

//  (1+2-3+4-5-)

M = Dmmp*Mmmpp+Dmpp*Mmppp;

Msq+=norm(M);

//  (1+2-3-4+5-)

M = Dmmm*Mmppp+Dmpm*Mmpmp;

Msq+=norm(M);

//  (1+2-3-4-5+)

M = Dmmm*Mmppp+Dmpm*Mmppm;

Msq+=norm(M);

//  (1-2+3+4-5-)

M = Dpmp*Mmppp+Dppp*Mpppp;

Msq+=norm(M);

//  (1-2+3-4+5-)

M = Dpmm*Mmppm+Dppm*Mmppp;

Msq+=norm(M);

//  (1-2+3-4-5+)

M = Dpmm*Mmpmp+Dppm*Mmppp;

Msq+=norm(M);

//  (1-2-3+4+5-)

M = Dmmp*Mmppm+Dmpp*Mmppp;

Msq+=norm(M);

//  (1-2-3+4-5+)

M = Dmmp*Mmpmp+Dmpp*Mmppp;

Msq+=norm(M);

//  (1-2-3-4+5+)

M = Dmmm*Mmppp+Dmpm*Mmmpp;

Msq+=norm(M);

  double pcrossGG=Nc/wm/z/zm*Msq;
  return -sigma0*as2pi*pcrossGG*fgg;
}

complex<double> ggNLO::Mppppp(int i1, int i2, int i3, int i4, int i5)
{
return -(S[i1][i2]*S[i2][i3]+S[i2][i3]*S[i3][i4]+S[i3][i4]*S[i4][i5]
       +S[i4][i5]*S[i5][i1]+S[i5][i1]*S[i1][i2]
       +Spb[i1][i2]*Spa[i2][i3]*Spb[i3][i4]*Spa[i4][i1]
       -Spa[i1][i2]*Spb[i2][i3]*Spa[i3][i4]*Spb[i4][i1])
        /Spa[i1][i2]/Spa[i2][i3]/Spa[i3][i4]/Spa[i4][i5]/Spa[i5][i1]/6.0;
}

complex<double> ggNLO::Mmmmmm(int i1, int i2, int i3, int i4, int i5)
{
return (S[i1][i2]*S[i2][i3]+S[i2][i3]*S[i3][i4]+S[i3][i4]*S[i4][i5]
       +S[i4][i5]*S[i5][i1]+S[i5][i1]*S[i1][i2]
       +Spa[i1][i2]*Spb[i2][i3]*Spa[i3][i4]*Spb[i4][i1]
       -Spb[i1][i2]*Spa[i2][i3]*Spb[i3][i4]*Spa[i4][i1])
        /Spb[i1][i2]/Spb[i2][i3]/Spb[i3][i4]/Spb[i4][i5]/Spb[i5][i1]/6.0;
}

complex<double> ggNLO::Mmpppp(int i1, int i2, int i3, int i4, int i5)
{
return -((S[i2][i3]+S[i3][i4]+S[i4][i5])*Spb[i2][i5]*Spb[i2][i5]
        -Spb[i2][i4]*Spa[i4][i3]*Spb[i3][i5]*Spb[i2][i5]
        -Spb[i1][i2]*Spb[i1][i5]/Spa[i1][i2]/Spa[i1][i5]*(
            Spa[i1][i2]*Spa[i1][i2]*Spa[i1][i3]*Spa[i1][i3]
                     *Spb[i2][i3]/Spa[i2][i3]
           +Spa[i1][i3]*Spa[i1][i3]*Spa[i1][i4]*Spa[i1][i4]
                     *Spb[i3][i4]/Spa[i3][i4]
           +Spa[i1][i4]*Spa[i1][i4]*Spa[i1][i5]*Spa[i1][i5]
                     *Spb[i4][i5]/Spa[i4][i5]))
        /Spb[i1][i2]/Spa[i2][i3]/Spa[i3][i4]/Spa[i4][i5]/Spb[i5][i1]/3.0;
}

complex<double> ggNLO::Mpmmmm(int i1, int i2, int i3, int i4, int i5)
{
return ((S[i2][i3]+S[i3][i4]+S[i4][i5])*Spa[i2][i5]*Spa[i2][i5]
        -Spa[i2][i4]*Spb[i4][i3]*Spa[i3][i5]*Spa[i2][i5]
        -Spa[i1][i2]*Spa[i1][i5]/Spb[i1][i2]/Spb[i1][i5]*(
            Spb[i1][i2]*Spb[i1][i2]*Spb[i1][i3]*Spb[i1][i3]
                     *Spa[i2][i3]/Spb[i2][i3]
           +Spb[i1][i3]*Spb[i1][i3]*Spb[i1][i4]*Spb[i1][i4]
                     *Spa[i3][i4]/Spb[i3][i4]
           +Spb[i1][i4]*Spb[i1][i4]*Spb[i1][i5]*Spb[i1][i5]
                     *Spa[i4][i5]/Spb[i4][i5]))
        /Spa[i1][i2]/Spb[i2][i3]/Spb[i3][i4]/Spb[i4][i5]/Spa[i5][i1]/3.0;
}

complex<double> ggNLO::Mmmppp(int i1, int i2, int i3, int i4, int i5)
{

  //  Note: nothing should depend on mm.

  double mm=1.0;

  complex<double> l23,l51;

  if (S[i2][i3]<0) {
     l23=log(mm/(-S[i2][i3]));
  } else {
     l23=log(mm/S[i2][i3])+I*PI;
  }

  if (S[i5][i1]<0) {
     l51=log(mm/(-S[i5][i1]));
  } else {
     l51=log(mm/S[i5][i1])+I*PI;
  }


  complex<double> Vf=-0.5*(l23+l51);
  complex<double> Vs=-1.0/3.0*Vf;

  complex<double> tree=pow(Spa[i1][i2],4)/Spa[i1][i2]/Spa[i2][i3]/Spa[i3][i4]
                      /Spa[i4][i5]/Spa[i5][i1];

  complex<double> Ff=-0.5*Spa[i1][i2]*Spa[i1][i2]*(
                    Spa[i2][i3]*Spb[i3][i4]*Spa[i4][i1]
                   +Spa[i2][i4]*Spb[i4][i5]*Spa[i5][i1])
              *L0(-S[i2][i3],-S[i5][i1])/S[i5][i1]
           /Spa[i2][i3]/Spa[i3][i4]/Spa[i4][i5]/Spa[i5][i1];

  complex<double> Fs=-1.0/3.0*Spb[i3][i4]*Spa[i4][i1]*Spa[i2][i4]*Spb[i4][i5]*(
                    Spa[i2][i3]*Spb[i3][i4]*Spa[i4][i1]
                   +Spa[i2][i4]*Spb[i4][i5]*Spa[i5][i1])
              *L2(-S[i2][i3],-S[i5][i1])/S[i5][i1]/S[i5][i1]/S[i5][i1]
           /Spa[i3][i4]/Spa[i4][i5]
             -1.0/3.0*Ff
             -1.0/3.0*Spa[i3][i5]*Spb[i3][i5]*Spb[i3][i5]*Spb[i3][i5]
                   /Spb[i1][i2]/Spb[i2][i3]/Spa[i3][i4]
                      /Spa[i4][i5]/Spb[i5][i1]
             +1.0/3.0*Spa[i1][i2]*Spb[i3][i5]*Spb[i3][i5]
                   /Spb[i2][i3]/Spa[i3][i4]
                      /Spa[i4][i5]/Spb[i5][i1]
             +1.0/6.0*Spa[i1][i2]*Spb[i3][i4]
                      *Spa[i4][i1]*Spa[i2][i4]*Spb[i4][i5]
                    /S[i2][i3]/Spa[i3][i4]/Spa[i4][i5]/S[i5][i1];

return -( (Vf+Vs)*tree +Ff+Fs );
}
complex<double> ggNLO::Mppmmm(int i1, int i2, int i3, int i4, int i5)
{

  //  Note: nothing should depend on mm.

  double mm=1.0;

  complex<double> l23,l51;

  if (S[i2][i3]<0) {
     l23=log(mm/(-S[i2][i3]));
  } else {
     l23=log(mm/S[i2][i3])+I*PI;
  }

  if (S[i5][i1]<0) {
     l51=log(mm/(-S[i5][i1]));
  } else {
     l51=log(mm/S[i5][i1])+I*PI;
  }


  complex<double> Vf=-0.5*(l23+l51);
  complex<double> Vs=-1.0/3.0*Vf;

  complex<double> tree=pow(Spb[i1][i2],4)/Spb[i1][i2]/Spb[i2][i3]/Spb[i3][i4]
                      /Spb[i4][i5]/Spb[i5][i1];

  complex<double> Ff=-0.5*Spb[i1][i2]*Spb[i1][i2]*(
                    Spb[i2][i3]*Spa[i3][i4]*Spb[i4][i1]
                   +Spb[i2][i4]*Spa[i4][i5]*Spb[i5][i1])
              *L0(-S[i2][i3],-S[i5][i1])/S[i5][i1]
           /Spb[i2][i3]/Spb[i3][i4]/Spb[i4][i5]/Spb[i5][i1];

  complex<double> Fs=-1.0/3.0*Spa[i3][i4]*Spb[i4][i1]*Spb[i2][i4]*Spa[i4][i5]*(
                    Spb[i2][i3]*Spa[i3][i4]*Spb[i4][i1]
                   +Spb[i2][i4]*Spa[i4][i5]*Spb[i5][i1])
              *L2(-S[i2][i3],-S[i5][i1])/S[i5][i1]/S[i5][i1]/S[i5][i1]
           /Spb[i3][i4]/Spb[i4][i5]
             -1.0/3.0*Ff
             -1.0/3.0*Spb[i3][i5]*Spa[i3][i5]*Spa[i3][i5]*Spa[i3][i5]
                   /Spa[i1][i2]/Spa[i2][i3]/Spb[i3][i4]
                      /Spb[i4][i5]/Spa[i5][i1]
             +1.0/3.0*Spb[i1][i2]*Spa[i3][i5]*Spa[i3][i5]
                   /Spa[i2][i3]/Spb[i3][i4]
                      /Spb[i4][i5]/Spa[i5][i1]
             +1.0/6.0*Spb[i1][i2]*Spa[i3][i4]
                      *Spb[i4][i1]*Spb[i2][i4]*Spa[i4][i5]
                    /S[i2][i3]/Spb[i3][i4]/Spb[i4][i5]/S[i5][i1];

return ( (Vf+Vs)*tree +Ff+Fs );
}

complex<double> ggNLO::Mmpmpp(int i1, int i2, int i3, int i4, int i5)
{

  //  Note: nothing should depend on mm.

  double mm=1.0;

  complex<double> l34,l51;

  if (S[i3][i4]<0) {
     l34=log(mm/(-S[i3][i4]));
  } else {
     l34=log(mm/S[i3][i4])+I*PI;
  }

  if (S[i5][i1]<0) {
     l51=log(mm/(-S[i5][i1]));
  } else {
     l51=log(mm/S[i5][i1])+I*PI;
  }

  complex<double> Vf=-0.5*(l34+l51);
  complex<double> Vs=-1.0/3.0*Vf;

  complex<double> tree=pow(Spa[i1][i3],4)/Spa[i1][i2]/Spa[i2][i3]/Spa[i3][i4]
                      /Spa[i4][i5]/Spa[i5][i1];

  complex<double> Ff=-Spa[i1][i3]*Spa[i1][i3]*Spa[i4][i1]*Spb[i2][i4]*Spb[i2][i4]
                  *Ls1(-S[i2][i3],-S[i5][i1],-S[i3][i4],-S[i5][i1])
                   /Spa[i4][i5]/Spa[i5][i1]/S[i5][i1]/S[i5][i1]
             +Spa[i1][i3]*Spa[i1][i3]*Spa[i5][i3]*Spb[i2][i5]*Spb[i2][i5]
                  *Ls1(-S[i1][i2],-S[i3][i4],-S[i5][i1],-S[i3][i4])
                   /Spa[i3][i4]/Spa[i4][i5]/S[i3][i4]/S[i3][i4]
             -0.5*Spa[i1][i3]*Spa[i1][i3]*Spa[i1][i3]*(
                    Spa[i1][i5]*Spb[i5][i2]*Spa[i2][i3]
                   -Spa[i3][i4]*Spb[i4][i2]*Spa[i2][i1])
              *L0(-S[i3][i4],-S[i5][i1])/S[i5][i1]
           /Spa[i1][i2]/Spa[i2][i3]/Spa[i3][i4]/Spa[i4][i5]/Spa[i5][i1];

  complex<double> Fs=-Spa[i1][i2]*Spa[i2][i3]*Spa[i3][i4]*Spa[i4][i1]
                      *Spa[i4][i1]*Spb[i2][i4]*Spb[i2][i4]*(
                2.0*Ls1(-S[i2][i3],-S[i5][i1],-S[i3][i4],-S[i5][i1])
                   +L1(-S[i2][i3],-S[i5][i1])
                   +L1(-S[i3][i4],-S[i5][i1]))
                  /Spa[i4][i5]/Spa[i5][i1]/Spa[i2][i4]/Spa[i2][i4]
                  /S[i5][i1]/S[i5][i1]
             +Spa[i3][i2]*Spa[i2][i1]*Spa[i1][i5]*Spa[i5][i3]
                      *Spa[i5][i3]*Spb[i2][i5]*Spb[i2][i5]*(
                2.0*Ls1(-S[i1][i2],-S[i3][i4],-S[i5][i1],-S[i3][i4])
                   +L1(-S[i1][i2],-S[i3][i4])
                   +L1(-S[i5][i1],-S[i3][i4]))
                  /Spa[i5][i4]/Spa[i4][i3]/Spa[i2][i5]/Spa[i2][i5]
                  /S[i3][i4]/S[i3][i4]
             +2.0/3.0*Spa[i2][i3]*Spa[i2][i3]
                     *Spa[i4][i1]*Spa[i4][i1]*Spa[i4][i1]
                     *Spb[i2][i4]*Spb[i2][i4]*Spb[i2][i4]
                     *L2(-S[i2][i3],-S[i5][i1])
                     /Spa[i4][i5]/Spa[i5][i1]/Spa[i2][i4]
                     /S[i5][i1]/S[i5][i1]/S[i5][i1]
             -2.0/3.0*Spa[i2][i1]*Spa[i2][i1]
                     *Spa[i5][i3]*Spa[i5][i3]*Spa[i5][i3]
                     *Spb[i2][i5]*Spb[i2][i5]*Spb[i2][i5]
                     *L2(-S[i1][i2],-S[i3][i4])
                     /Spa[i5][i4]/Spa[i4][i3]/Spa[i2][i5]
                     /S[i3][i4]/S[i3][i4]/S[i3][i4]
             +L2(-S[i3][i4],-S[i5][i1])/S[i5][i1]/S[i5][i1]/S[i5][i1]*(
	        1.0/3.0*Spa[i1][i3]*Spb[i2][i4]*Spb[i2][i5]*(
                    Spa[i1][i5]*Spb[i5][i2]*Spa[i2][i3]
                   -Spa[i3][i4]*Spb[i4][i2]*Spa[i2][i1])/Spa[i4][i5]
               +2.0/3.0*Spa[i1][i2]*Spa[i1][i2]*Spa[i3][i4]*Spa[i3][i4]
                    *Spa[i4][i1]*Spb[i2][i4]*Spb[i2][i4]*Spb[i2][i4]
                    /Spa[i4][i5]/Spa[i5][i1]/Spa[i2][i4]
               -2.0/3.0*Spa[i3][i2]*Spa[i3][i2]*Spa[i1][i5]*Spa[i1][i5]
                    *Spa[i5][i3]*Spb[i2][i5]*Spb[i2][i5]*Spb[i2][i5]
                    /Spa[i5][i4]/Spa[i4][i3]/Spa[i2][i5])
            +1.0/6.0*Spa[i1][i3]*Spa[i1][i3]*Spa[i1][i3]*(
                    Spa[i1][i5]*Spb[i5][i2]*Spa[i2][i3]
                   -Spa[i3][i4]*Spb[i4][i2]*Spa[i2][i1])
              *L0(-S[i3][i4],-S[i5][i1])/S[i5][i1]
           /Spa[i1][i2]/Spa[i2][i3]/Spa[i3][i4]/Spa[i4][i5]/Spa[i5][i1]
             +1.0/3.0*Spb[i2][i4]*Spb[i2][i4]*Spb[i2][i5]*Spb[i2][i5]
           /Spb[i1][i2]/Spb[i2][i3]/Spb[i3][i4]/Spa[i4][i5]/Spb[i5][i1]
             -1.0/3.0*Spa[i1][i2]*Spa[i4][i1]*Spa[i4][i1]
                     *Spb[i2][i4]*Spb[i2][i4]*Spb[i2][i4]
                   /Spa[i4][i5]/Spa[i5][i1]/Spa[i2][i4]
                      /Spb[i2][i3]/Spb[i3][i4]/S[i5][i1]
             +1.0/3.0*Spa[i3][i2]*Spa[i5][i3]*Spa[i5][i3]
                     *Spb[i2][i5]*Spb[i2][i5]*Spb[i2][i5]
                   /Spa[i5][i4]/Spa[i4][i3]/Spa[i2][i5]
                      /Spb[i2][i1]/Spb[i1][i5]/S[i3][i4]
             +1.0/6.0*Spa[i1][i3]*Spa[i1][i3]*Spb[i2][i4]*Spb[i2][i5]
                      /S[i3][i4]/Spa[i4][i5]/S[i5][i1];

return -( (Vf+Vs)*tree +Ff+Fs );
}

complex<double> ggNLO::Mpmpmm(int i1, int i2, int i3, int i4, int i5)
{

  //  Note: nothing should depend on mm.

  double mm=1.0;

  complex<double> l34,l51;

  if (S[i3][i4]<0) {
     l34=log(mm/(-S[i3][i4]));
  } else {
     l34=log(mm/S[i3][i4])+I*PI;
  }

  if (S[i5][i1]<0) {
     l51=log(mm/(-S[i5][i1]));
  } else {
     l51=log(mm/S[i5][i1])+I*PI;
  }

  complex<double> Vf=-0.5*(l34+l51);
  complex<double> Vs=-1.0/3.0*Vf;

  complex<double> tree=pow(Spb[i1][i3],4)/Spb[i1][i2]/Spb[i2][i3]/Spb[i3][i4]
                      /Spb[i4][i5]/Spb[i5][i1];

  complex<double> Ff=-Spb[i1][i3]*Spb[i1][i3]*Spb[i4][i1]*Spa[i2][i4]*Spa[i2][i4]
                  *Ls1(-S[i2][i3],-S[i5][i1],-S[i3][i4],-S[i5][i1])
                   /Spb[i4][i5]/Spb[i5][i1]/S[i5][i1]/S[i5][i1]
             +Spb[i1][i3]*Spb[i1][i3]*Spb[i5][i3]*Spa[i2][i5]*Spa[i2][i5]
                  *Ls1(-S[i1][i2],-S[i3][i4],-S[i5][i1],-S[i3][i4])
                   /Spb[i3][i4]/Spb[i4][i5]/S[i3][i4]/S[i3][i4]
             -0.5*Spb[i1][i3]*Spb[i1][i3]*Spb[i1][i3]*(
                    Spb[i1][i5]*Spa[i5][i2]*Spb[i2][i3]
                   -Spb[i3][i4]*Spa[i4][i2]*Spb[i2][i1])
              *L0(-S[i3][i4],-S[i5][i1])/S[i5][i1]
           /Spb[i1][i2]/Spb[i2][i3]/Spb[i3][i4]/Spb[i4][i5]/Spb[i5][i1];

  complex<double> Fs=-Spb[i1][i2]*Spb[i2][i3]*Spb[i3][i4]*Spb[i4][i1]
                      *Spb[i4][i1]*Spa[i2][i4]*Spa[i2][i4]*(
                2.0*Ls1(-S[i2][i3],-S[i5][i1],-S[i3][i4],-S[i5][i1])
                   +L1(-S[i2][i3],-S[i5][i1])
                   +L1(-S[i3][i4],-S[i5][i1]))
                  /Spb[i4][i5]/Spb[i5][i1]/Spb[i2][i4]/Spb[i2][i4]
                  /S[i5][i1]/S[i5][i1]
             +Spb[i3][i2]*Spb[i2][i1]*Spb[i1][i5]*Spb[i5][i3]
                      *Spb[i5][i3]*Spa[i2][i5]*Spa[i2][i5]*(
                2.0*Ls1(-S[i1][i2],-S[i3][i4],-S[i5][i1],-S[i3][i4])
                   +L1(-S[i1][i2],-S[i3][i4])
                   +L1(-S[i5][i1],-S[i3][i4]))
                  /Spb[i5][i4]/Spb[i4][i3]/Spb[i2][i5]/Spb[i2][i5]
                  /S[i3][i4]/S[i3][i4]
             +2.0/3.0*Spb[i2][i3]*Spb[i2][i3]
                     *Spb[i4][i1]*Spb[i4][i1]*Spb[i4][i1]
                     *Spa[i2][i4]*Spa[i2][i4]*Spa[i2][i4]
                     *L2(-S[i2][i3],-S[i5][i1])
                     /Spb[i4][i5]/Spb[i5][i1]/Spb[i2][i4]
                     /S[i5][i1]/S[i5][i1]/S[i5][i1]
             -2.0/3.0*Spb[i2][i1]*Spb[i2][i1]
                     *Spb[i5][i3]*Spb[i5][i3]*Spb[i5][i3]
                     *Spa[i2][i5]*Spa[i2][i5]*Spa[i2][i5]
                     *L2(-S[i1][i2],-S[i3][i4])
                     /Spb[i5][i4]/Spb[i4][i3]/Spb[i2][i5]
                     /S[i3][i4]/S[i3][i4]/S[i3][i4]
             +L2(-S[i3][i4],-S[i5][i1])/S[i5][i1]/S[i5][i1]/S[i5][i1]*(
	        1.0/3.0*Spb[i1][i3]*Spa[i2][i4]*Spa[i2][i5]*(
                    Spb[i1][i5]*Spa[i5][i2]*Spb[i2][i3]
                   -Spb[i3][i4]*Spa[i4][i2]*Spb[i2][i1])/Spb[i4][i5]
               +2.0/3.0*Spb[i1][i2]*Spb[i1][i2]*Spb[i3][i4]*Spb[i3][i4]
                    *Spb[i4][i1]*Spa[i2][i4]*Spa[i2][i4]*Spa[i2][i4]
                    /Spb[i4][i5]/Spb[i5][i1]/Spb[i2][i4]
               -2.0/3.0*Spb[i3][i2]*Spb[i3][i2]*Spb[i1][i5]*Spb[i1][i5]
                    *Spb[i5][i3]*Spa[i2][i5]*Spa[i2][i5]*Spa[i2][i5]
                    /Spb[i5][i4]/Spb[i4][i3]/Spb[i2][i5])
            +1.0/6.0*Spb[i1][i3]*Spb[i1][i3]*Spb[i1][i3]*(
                    Spb[i1][i5]*Spa[i5][i2]*Spb[i2][i3]
                   -Spb[i3][i4]*Spa[i4][i2]*Spb[i2][i1])
              *L0(-S[i3][i4],-S[i5][i1])/S[i5][i1]
           /Spb[i1][i2]/Spb[i2][i3]/Spb[i3][i4]/Spb[i4][i5]/Spb[i5][i1]
             +1.0/3.0*Spa[i2][i4]*Spa[i2][i4]*Spa[i2][i5]*Spa[i2][i5]
           /Spa[i1][i2]/Spa[i2][i3]/Spa[i3][i4]/Spb[i4][i5]/Spa[i5][i1]
             -1.0/3.0*Spb[i1][i2]*Spb[i4][i1]*Spb[i4][i1]
                     *Spa[i2][i4]*Spa[i2][i4]*Spa[i2][i4]
                   /Spb[i4][i5]/Spb[i5][i1]/Spb[i2][i4]
                      /Spa[i2][i3]/Spa[i3][i4]/S[i5][i1]
             +1.0/3.0*Spb[i3][i2]*Spb[i5][i3]*Spb[i5][i3]
                     *Spa[i2][i5]*Spa[i2][i5]*Spa[i2][i5]
                   /Spb[i5][i4]/Spb[i4][i3]/Spb[i2][i5]
                      /Spa[i2][i1]/Spa[i1][i5]/S[i3][i4]
             +1.0/6.0*Spb[i1][i3]*Spb[i1][i3]*Spa[i2][i4]*Spa[i2][i5]
                      /S[i3][i4]/Spb[i4][i5]/S[i5][i1];

return ( (Vf+Vs)*tree +Ff+Fs );
}


complex<double> L0(double r1, double r2)
{

  double r=r1/r2;

  complex<double> lr;

  if (r<0) {
      if (r1<0) {
          lr=log(-r)-I*PI;
      } else {
          lr=log(-r)+I*PI;
      }
  } else {
      lr=log(r);
  }

  return lr/(1.0-r);
}

complex<double> L1(double r1, double r2)
{
  double r=r1/r2;

  complex<double> lr;

  if (r<0) {
      if (r1<0) {
          lr=log(-r)-I*PI;
      } else {
          lr=log(-r)+I*PI;
      }
  } else {
      lr=log(r);
  }

  double rm=1.0-r;

  return (lr+rm)/rm/rm;
}

complex<double> L2(double r1, double r2)
{

  double r=r1/r2;
  double rm=1.0-r;

  complex<double> lr;

  if (r<0) {
      if (r1<0) {
          lr=log(-r)-I*PI;
      } else {
          lr=log(-r)+I*PI;
      }
  } else {
      lr=log(r);
  }

  return (lr-(r-1.0/r)/2.0)/rm/rm/rm;
}

complex<double> Ls1(double r1, double r2, double r3, double r4)
{
  double R1=r1/r2;
  double R2=r3/r4;

  double r1m=1.0-R1;
  double r2m=1.0-R2;
  double rm=1.0-R1-R2;

  complex<double> lr1, lr2;

  complex<double> cLi21, cLi22;

  if (R1<0) {
      if (r1<0) {
          lr1=log(-R1)-I*PI;
          cLi21=I*PI*log(r1m);
      } else {
          lr1=log(-R1)+I*PI;
          cLi21=-I*PI*log(r1m);
      }
  } else {
      lr1=log(R1);
      cLi21=0.0;
  }

  if (R2<0) {
      if (r3<0) {
          lr2=log(-R2)-I*PI;
          cLi22=I*PI*log(r2m);
      } else {
          lr2=log(-R2)+I*PI;
          cLi22=-I*PI*log(r2m);
      }
  } else {
      lr2=log(R2);
      cLi22=0.0;
  }

  return (Li2(r1m)+cLi21+Li2(r2m)+cLi22
      +lr1*lr2-PI*PI/6.0+rm*(L0(r1,r2)+L0(r3,r4)))
             /rm/rm;
}

