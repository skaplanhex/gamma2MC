#include "HiggsNLO.h"


void HiggsNLO::initializeSubprocess(void)
{ 
  pdfpairs(col,xA,xB,muF,fgg,fgq,fqg,fqiqj,fqiqi,fqiai);
  long int iparton=0;
  fgg0=pdf(iparton,xA0,muF)*pdf(iparton,xB0,muF);
  pgg=Nc*(1.0+pow(zm,4)+pow(z,4))/z/zm;
  pgq=CF*(1.0+zm*zm)/z;

}


double HiggsNLO::wgtLO(void)
{ 
  return sigLO*fgg0;
}

double HiggsNLO::wgtV(void)
{ 
  double nloV=2.0/3.0*PI*PI*Nc+11.0+beta0*logRF
              +4.0*Nc*logzmin*(logzmin+logHF);
  return wgtLO()*as2pi*nloV;
}


double HiggsNLO::wgtR(void)
{ 
  double pcrossGG=Nc*(1.0+pow(z,4)+pow(zm*w,4)+pow(zm*wm,4))/w/wm/z/zm;
  double pcrossQG=CF*(1.0+zm*zm*wm*wm)/z/w;
  double pcrossGQ=CF*(1.0+zm*zm*w*w)/z/wm;
  double pcrossQA=32.0/9.0*zm*zm*zm/z*(w*w+wm*wm);
  return sigLO*as2pi*
              (pcrossGG*fgg+pcrossQA*fqiai+pcrossQG*fqg+pcrossGQ*fgq);
}

double HiggsNLO::wgtVplus(void)
{ 
  double pcrossGG=4.0*Nc*(2.0*logzm+logHF)/zm;
  return -sigLO*as2pi*z*
              (pcrossGG*fgg0);
}

double HiggsNLO::wgtVa(void)
{ 
  double L=2.0*logzm+logHF-logz;
  double pcrossGG=pgg*L;
  double pcrossQG=pgq*L+CF*z;
  return sigLO*as2pi*
              (pcrossGG*fgg+pcrossQG*fqg);
}

double HiggsNLO::wgtVb(void)
{ 
  double L=2.0*logzm+logHF-logz;
  double pcrossGG=pgg*L;
  double pcrossGQ=pgq*L+CF*z;
  return sigLO*as2pi*
              (pcrossGG*fgg+pcrossGQ*fgq);
}

double HiggsNLO::wgtSa(void)
{ 
  double pcrossGG=pgg/w;
  double pcrossQG=pgq/w;
  return -sigLO*as2pi*
              (pcrossGG*fgg+pcrossQG*fqg);
}

double HiggsNLO::wgtSb(void)
{ 
  double pcrossGG=pgg/wm;
  double pcrossGQ=pgq/wm;
  return -sigLO*as2pi*
              (pcrossGG*fgg+pcrossGQ*fgq);
}

double HiggsNLO::A1sq(void)
{
  double mf=MT;
  double mfsq=mf*mf;
  return norm(mfsq/Qsq*(4.0-W2(Qsq)*(1.0-4.0*mfsq/Qsq)));
}

complex<double> HiggsNLO::W2(double s)
{
  double mf=MT;
  complex<double> WW;
  if (s<0.0) {
     WW=4.0*pow(asinh(sqrt(-s)/2.0/mf),2);
  } else if (s<4.0*mf*mf) {
     WW=-4.0*pow(asin(sqrt(s)/2.0/mf),2);
  } else {
     WW=4.0*acosh(sqrt(s)/2.0/mf)*(acosh(sqrt(s)/2.0/mf)-I*PI)-PI*PI;
  }
  return WW;
}

double muF0;
double tau0;
double Qsq0;
double zmin(double xa);
double zmax(double xa);
double dXsecLO(double xa);
double dXsecNLO(double xa, double z);
collider coll;

double zmin(double zp)
{
return tau0/zp;
}

double zmax(double zp)
{
return 1.0;
}

double dXsecLO(double zp)
{
  double xa=zp;
  double xb=tau0/xa;
  long int iparton=0;
  return pdf(iparton,xa,muF0)*pdf(iparton,xb,muF0)/zp;
}

double dXsecNLO(double zp, double z)
{
  double fgg,fgq,fqg,fqiqj,fqiqi,fqiai;
  double xa=zp;
  double xb=tau0/xa/z;
  double xb0=tau0/xa;
  pdfpairs(coll,xa,xb,muF0,fgg,fgq,fqg,fqiqj,fqiqi,fqiai);
  long int iparton=0;
  double fgg0=pdf(iparton,xa,muF0)*pdf(iparton,xb0,muF0);
  double zm=1.0-z;
  double pgg=Nc*(1.0+pow(zm,4)+pow(z,4))/z/zm;
  double pgq=CF*(1.0+zm*zm)/z;
  double lh=log(Qsq0/muF0/muF0);
  double zmmin=1.0-zmin(zp);
  double lzmmin=log(zmmin);
  double L=2.0*log(zm)+lh-log(z);
  double L0=2.0*log(zm)+lh;
  double pcrossGG=2.0*pgg*L-11.0/3.0*Nc*zm*zm*zm/z;
  double pcrossQA=64.0/27.0*zm*zm*zm/z;
  double pcrossQG=pgq*L+CF*z-CF*3.0/2.0*zm*zm/z;
  double pcrossGG0=4.0*Nc*(L0/zm-lzmmin/zmmin*(lzmmin+lh));
  //The last term contains the integral z=0..zmin of the
  //subtraction term in the plus function.  It has then been
  //included as a constant in z by int 1/(1-zmin)*int(z=zmin..1).
  return (pcrossGG*fgg+pcrossQA*fqiai+pcrossQG*(fqg+fgq))/zp/z
         -pcrossGG0*fgg0/zp;


}

double HiggsNLO::XsecLO(void)
{
  muF0=muF;
  tau0=tau;
  double one=1.0;
  return sigLO*tau*gauss(dXsecLO,tau,one);
}

double HiggsNLO::XsecNLO(void)
{
  tau0=tau;
  muF0=muF;
  Qsq0=Qsq;
  coll=col;
  double one=1.0;
  double nloV=2.0/3.0*PI*PI*Nc+11.0+beta0*logRF;
  return sigLO*tau*( gauss(dXsecLO,tau,one)*(1.0+as2pi*nloV)
                +as2pi*digauss(dXsecNLO,tau,one,zmin,zmax) );
}

double HiggsNLO::XsecLOlargeMt(void)
{
  return XsecLO()*sigma0/sigLO;
}

double HiggsNLO::XsecNLOlargeMt(void)
{
  return XsecNLO()*sigma0/sigLO;
}

