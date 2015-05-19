#include "gamma2NLO.h"

void gamma2NLO::setParameters(void)
{ 
  Q=Q0+Qbin*(x[7]-0.5);
  Qsq=Q*Q;
  tau=Qsq/ss; 
  muR=chiR*Q;
  muF=chiF*Q;
  mufr=chifr*Q;
  logHF=log(Qsq/muF/muF);
  logHfr=log(Qsq/mufr/mufr);
  logRF=log(muR*muR/muF/muF);
  as2pi=alphas(muR)/2.0/PI;
  setSigma0();
  xA0=xA=pow(tau,x[1]);
  w1=x[2];
  phi1=2.0*PI*x[3];
  z=pow(tau/xA0,x[4]);
  w=x[5];
  phi=phi1+2.0*PI*x[6];
  xB0=tau/xA;
  xB=xB0/z;

  logzmin=log(1.0-xB0);
  jacob0=(2.0*Q)*tau*(-log(tau));
  jacob=jacob0*(-log(tau/xA));
  zm=1.0-z;
  wm=1.0-w;
  w1m=1.0-w1;
  logz=log(z);
  logzm=log(zm);
}

void gamma2NLO::buildProducts(void)
{ 
  double zr=sqrt(z);
  double zmr=sqrt(zm);
  double wr=sqrt(w);
  double wmr=sqrt(wm);
  double w1r=sqrt(w1);
  double w1mr=sqrt(w1m);
  complex<double> eplus=exp(I*phi/2.0);
  complex<double> eminus=exp(-I*phi/2.0);
  complex<double> eplus1=exp(I*(phi-phi1)/2.0);
  complex<double> eminus1=exp(-I*(phi-phi1)/2.0);

  //  1=pA, 2=pB, 3=p1, 4=k1, 5=k2  (4 and 5 are photons)


  //  Note that we should be able to remove the overall eplus
  //  and eminus phases by changing the phases of |1> and |2>.

  // The Spa[i][j] and Spb[i][j] are normalized by rshat=sqrt(ss*xA*xB).
  // The S[i][j] are normalized by shat=ss*xA*xB.

  Spa[1][2]=1.0;
  Spa[1][3]=zmr*wr*eplus;
  Spa[2][3]=-zmr*wmr*eminus;
  Spa[1][4]=eplus*(wr*wmr*w1mr*(zr-1.0)*eplus1+(w*zr+wm)*w1r*eminus1);
  Spa[2][4]=-eminus*(wr*wmr*w1r*(zr-1.0)*eminus1+(wm*zr+w)*w1mr*eplus1);
  Spa[3][4]=zmr*(-w1mr*wr*eplus1+wmr*w1r*eminus1);
  Spa[1][5]=eplus*(-wr*wmr*w1r*(zr-1.0)*eplus1+(w*zr+wm)*w1mr*eminus1);
  Spa[2][5]=-eminus*(wr*wmr*w1mr*(zr-1.0)*eminus1-(wm*zr+w)*w1r*eplus1);
  Spa[3][5]=zmr*(w1r*wr*eplus1+wmr*w1mr*eminus1);
  Spa[4][5]=zr;

  for(int i=1;i<6;i++)
    for (int j=i+1;j<6;j++) {
       Spa[j][i]=-Spa[i][j];
       if ((i<3)&&(j!=2)) {
            Spb[j][i]=-conj(Spa[i][j]);
       } else {
            Spb[j][i]=conj(Spa[i][j]);
       }            
       Spb[i][j]=-Spb[j][i];
       S[i][j]=S[j][i]=real(Spa[i][j]*Spb[j][i]);
    }
}

void gamma2NLO::buildParticles(void)
{ 
  k10=k20=photon;         
  double cosTheta1=1.0-2.0*w1;
  k10.setMom(Q/2.0,cosTheta1,phi1);
  k20.setMom(Q/2.0,-cosTheta1,phi1+PI);
  k11=k1A=k1B=k10;
  k21=k2A=k2B=k20;
  double b=(xA0-xB0)/(xA0+xB0);
  k10.boostZ(b);
  k20.boostZ(b);

  b=(xA*z-xB)/(xA*z+xB);
  k1A.boostZ(b);
  k2A.boostZ(b);
  b=(xA-xB*z)/(xA+xB*z);
  k1B.boostZ(b);
  k2B.boostZ(b);
  double cosTheta=1.0-2.0*w;
  double rshat=sqrt(ss*xA*xB);
  double energy=(1.0-z)*rshat/2.0;
  p1=particle(0.0,0.0,     " parton   ");
  p1.setMom(energy,cosTheta,phi);

  double theta=acos(cosTheta);
  k11.invBoost(0.0,theta,phi);
  k21.invBoost(0.0,theta,phi);
  b=(z-1.0)/(1.0+z);
  k11.boost(b,theta,phi);
  k21.boost(b,theta,phi);
  b=(xA-xB)/(xA+xB);
  p1.boostZ(b);
  k11.boostZ(b);
  k21.boostZ(b);
}

double gamma2NLO::fxn(void)
{
  initialize();
  double xwgt=0.0;
  switch ( evt ) {
      case NLO :
             xwgtLO=wgtLO()*jacob0;
             xwgtV=wgtV()*jacob0;
             xwgtVplus=wgtVplus()*jacob;
             xwgtVa=wgtVa()*jacob;
             xwgtVb=wgtVb()*jacob;
             xwgtSa=wgtSa()*jacob;
             xwgtSb=wgtSb()*jacob;
             xwgtR=wgtR()*jacob;
             if (!cc->cut(k10,k20))   { xwgt+= xwgtLO+xwgtV+xwgtVplus; }
             if (!cc->cut(k1A,k2A))   { xwgt+= xwgtSa+xwgtVa; }
             if (!cc->cut(k1B,k2B))   { xwgt+= xwgtSb+xwgtVb; }
             if (!cc->cut(k11,k21,p1)) { xwgt+= xwgtR; }
             break;
      case HARDNLO :
             xwgtSa=wgtSa()*jacob;
             xwgtSb=wgtSb()*jacob;
             xwgtR=wgtR()*jacob;
             if (!cc->cut(k1A,k2A))   { xwgt+= xwgtSa; }
             if (!cc->cut(k1B,k2B))   { xwgt+= xwgtSb; }
             if (!cc->cut(k11,k21,p1)) { xwgt+= xwgtR; }
             break;
      case LO :
             xwgtLO=wgtLO()*jacob0;
             if (!cc->cut(k10,k20))   { xwgt+= xwgtLO; }
             break;
      case EASYNLO :
             xwgtLO=wgtLO()*jacob0;
             xwgtV=wgtV()*jacob0;
             xwgtVplus=wgtVplus()*jacob;
             xwgtVa=wgtVa()*jacob;
             xwgtVb=wgtVb()*jacob;
             if (!cc->cut(k10,k20))   { xwgt+= xwgtLO+xwgtV+xwgtVplus; }
             if (!cc->cut(k1A,k2A))   { xwgt+= xwgtVa; }
             if (!cc->cut(k1B,k2B))   { xwgt+= xwgtVb; }
             break;
  };
  return wgtnorm()*xwgt;    
}



void gamma2NLO::fillBins(void)
{
  int n=0;
  double cts[4];
  double wgt[4];
  double wn=wgtnorm()*vwgt;
  switch ( evt ) {
      case NLO :
             if (!cc->cut(k10,k20))   { 
	       cts[n] = dist->Param(k10,k20);
               wgt[n] = wn*(xwgtLO+xwgtV+xwgtVplus);
               n++;
             } 
             if (!cc->cut(k1A,k2A))   { 
	       cts[n] = dist->Param(k1A,k2A);
               wgt[n] = wn*(xwgtSa+xwgtVa);
               n++;
             }
             if (!cc->cut(k1B,k2B))   { 
	       cts[n] = dist->Param(k1B,k2B);
               wgt[n] = wn*(xwgtSb+xwgtVb); 
               n++;
             }
             if (!cc->cut(k11,k21,p1)) { 
	       cts[n] = dist->Param(k11,k21);
               wgt[n] = wn*xwgtR; 
               n++;
             }
             break;
      case HARDNLO :
             if (!cc->cut(k1A,k2A))   { 
	       cts[n] = dist->Param(k1A,k2A);
               wgt[n] = wn*xwgtSa; 
               n++;
             }
             if (!cc->cut(k1B,k2B))   { 
	       cts[n] = dist->Param(k1B,k2B);
               wgt[n] = wn*xwgtSb; 
               n++;
             }
             if (!cc->cut(k11,k21,p1)) { 
	       cts[n] = dist->Param(k11,k21);
               wgt[n] = wn*xwgtR; 
               n++;
             }
             break;
      case LO :
             if (!cc->cut(k10,k20))   { 
	       cts[n] = dist->Param(k10,k20);
               wgt[n] = wn*xwgtLO; 
               n++;
             }
             break;
      case EASYNLO :
             if (!cc->cut(k10,k20))   { 
	       cts[n] = dist->Param(k10,k20);
               wgt[n] = wn*(xwgtLO+xwgtV+xwgtVplus); 
               n++;
             }
             if (!cc->cut(k1A,k2A))   { 
	       cts[n] = dist->Param(k1A,k2A);
               wgt[n] = wn*xwgtVa; 
               n++;
             }
             if (!cc->cut(k1B,k2B))   { 
	       cts[n] = dist->Param(k1B,k2B);
               wgt[n] = wn*xwgtVb; 
               n++;
             }
             break;
  };
  dist->bin(n,cts,wgt);
}


