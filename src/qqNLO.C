#include "qqNLO.h"

void qqNLO::initializeSubprocess(void)
{ 
  long int iparton,aparton;
  double e4, pdfiA, pdfiB, pdfaA, pdfaB;  
  fqq0=fqq=fgq=fqg=0.0;

  for (iparton=1;iparton<6;iparton++) {
    if ((iparton==1)||(iparton==4)) {
         e4=16.0/81.0;
    } else {
         e4=1.0/81.0;
    }
    aparton=-iparton;
    fqq0+=e4*(pdf(iparton,xA0,muF)*pdf(aparton,xB0,muF)
             +pdf(aparton,xA0,muF)*pdf(iparton,xB0,muF));
    pdfiA=pdf(iparton,xA,muF);
    pdfiB=pdf(iparton,xB,muF);
    pdfaA=pdf(aparton,xA,muF);
    pdfaB=pdf(aparton,xB,muF);
    fqq+=e4*(pdfiA*pdfaB
             +pdfaA*pdfiB);
    fqg+=e4*(pdfiA+pdfaA);
    fgq+=e4*(pdfiB+pdfaB);
  }
  iparton=0;
  fqg*=pdf(iparton,xB,muF);
  fgq*=pdf(iparton,xA,muF);

  sigLO=sigma0*(w1m*w1m+w1*w1)/w1/w1m;
  sigLOqg=sigma0*(1.0+wm*wm)/wm;
  sigLOgq=sigma0*(1.0+w*w)/w;

  // Note that I use w rather than w1 to generate
  // the LO kinematics in the fragmentation term.

  pqq=CF*(1.0+z*z)/zm;
  pqg=0.5*(z*z+zm*zm);
  pgammaq=(1.0+zm*zm)/z;

  x1=-w1m;
  y=-w1;
  X1=log(-x1);
  Y=log(-y);

  if (qqflag==QGONLY) {  fqq0=fqq=0.0;}
  if (qqflag==QQONLY) {  fgq=fqg=0.0;}

  double b=(xA-xB)/(xA+xB);

  k1f=k2f=photon;       
  pf=p1;
  double cosTheta=1.0-2.0*w;
  double rshat=sqrt(ss*xA*xB);
  k1f.setMom(z*rshat/2.0,cosTheta,phi);
  k2f.setMom(rshat/2.0,-cosTheta,phi+PI);
  k1f.boostZ(b);
  k2f.boostZ(b);


}


double qqNLO::wgtLO(void)
{ 
  return sigLO*fqq0;
}

double qqNLO::wgtV(void)
{ 
  double nloV=CF*(2.0/3.0*PI*PI+3.0*logHF
		  -7.0+1.0/(x1*x1+y*y)*(
		  x1*x1*( X1*X1 +(1.0-2.0/x1)*Y +Y*Y/x1/x1 )
		 +y*y*  ( Y*Y   +(1.0-2.0/y)*X1 +X1*X1/y/y ) )
              +4.0*logzmin*(logzmin+logHF));
  return wgtLO()*as2pi*nloV;
}


double qqNLO::wgtR(void)
{ 
  double pcrossQQ=CF*S[1][2]/(S[1][3]*S[1][4]*S[1][5]
			     *S[2][3]*S[2][4]*S[2][5])
	          *(S[1][3]*S[2][3]*S[2][3]*S[2][3]
	           +S[1][4]*S[2][4]*S[2][4]*S[2][4]
	           +S[1][5]*S[2][5]*S[2][5]*S[2][5]
	           +S[2][3]*S[1][3]*S[1][3]*S[1][3]
	           +S[2][4]*S[1][4]*S[1][4]*S[1][4]
	           +S[2][5]*S[1][5]*S[1][5]*S[1][5]);

  double pcrossQG=-0.5*S[1][3]/(S[1][2]*S[1][4]*S[1][5]
                               *S[2][3]*S[3][4]*S[3][5])
		  *(S[1][2]*S[2][3]*S[2][3]*S[2][3]
                   +S[1][4]*S[3][4]*S[3][4]*S[3][4]
                   +S[1][5]*S[3][5]*S[3][5]*S[3][5]
		   +S[2][3]*S[1][2]*S[1][2]*S[1][2]
                   +S[3][4]*S[1][4]*S[1][4]*S[1][4]
                   +S[3][5]*S[1][5]*S[1][5]*S[1][5]);

  double pcrossGQ=-0.5*S[2][3]/(S[1][2]*S[2][4]*S[2][5]
                               *S[1][3]*S[3][4]*S[3][5])
		  *(S[1][2]*S[1][3]*S[1][3]*S[1][3]
                   +S[2][4]*S[3][4]*S[3][4]*S[3][4]
                   +S[2][5]*S[3][5]*S[3][5]*S[3][5]
		   +S[1][3]*S[1][2]*S[1][2]*S[1][2]
                   +S[3][4]*S[2][4]*S[2][4]*S[2][4]
                   +S[3][5]*S[2][5]*S[2][5]*S[2][5]);

  return sigma0*as2pi*z*zm*
              (pcrossQQ*fqq+pcrossQG*fqg+pcrossGQ*fgq);
}

double qqNLO::wgtVplus(void)
{ 
  double pcrossQQ=4.0*CF*(2.0*logzm+logHF)/zm;
  return -sigLO*as2pi*z*
              (pcrossQQ*fqq0);
}

double qqNLO::wgtVa(void)
{ 
  double L=2.0*logzm+logHF-logz;
  double pcrossQQ=pqq*L+CF*zm;
  double pcrossGQ=pqg*L+z*zm;
  return sigLO*as2pi*
              (pcrossQQ*fqq+pcrossGQ*fgq);
}

double qqNLO::wgtVb(void)
{ 
  double L=2.0*logzm+logHF-logz;
  double pcrossQQ=pqq*L+CF*zm;
  double pcrossQG=pqg*L+z*zm;
  return sigLO*as2pi*
              (pcrossQQ*fqq+pcrossQG*fqg);
}

double qqNLO::wgtSa(void)
{ 
  double pcrossQQ=pqq/w;
  double pcrossGQ=pqg/w;
  return -sigLO*as2pi*
              (pcrossQQ*fqq+pcrossGQ*fgq);
}

double qqNLO::wgtSb(void)
{ 
  double pcrossQQ=pqq/wm;
  double pcrossQG=pqg/wm;
  return -sigLO*as2pi*
              (pcrossQQ*fqq+pcrossQG*fqg);
}

double qqNLO::wgtSf(void)
{ 

  // Note that this contains the collinear singularities
  // with both photon 1 and 2 and includes the 1/2 compared
  // to the qqbar->gamma gamma matrix element, in addition
  // to the 1/2 Bose factor, contained in Sigma0.

  double ww1 = S[3][4]/zm;

  return -as2pi*z*pgammaq*
              (sigLOqg*fqg+sigLOgq*fgq)/2.0/ww1/(1.0-ww1);
}

double qqNLO::wgtVf(void)
{ 

  double L=2.0*logzm+logHfr;

  return as2pi*z*(pgammaq*L+z)*
              (sigLOqg*fqg+sigLOgq*fgq);

  // The factor of z comes from the fact that the flux 
  // is 1/shat, not 1/(shat*z) for fragmentation.
  // Note that the matrix element has a 1/2 compared to
  // the qqbar->gamma gamma matrix element.  We have then
  // multiplied by this by 2, since the identification of photon 1
  // with the fragmented photon breaks the Bose symmetry
  // (which had been included in sigma0).

}



double qqNLO::fxn(void)
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
             xwgtSf=wgtSf()*jacob;
             xwgtVf=wgtVf()*jacob;
             if (!cc->cut(k10,k20))   { xwgt+= xwgtLO+xwgtV+xwgtVplus; }
             if (!cc->cut(k1A,k2A))   { xwgt+= xwgtSa+xwgtVa; }
             if (!cc->cut(k1B,k2B))   { xwgt+= xwgtSb+xwgtVb; }
             if (!cc->cut(k11,k21,p1)) { xwgt+= xwgtR; }
             if (!cc->cut(k1f,k2f,pf)) { xwgt+= xwgtSf+xwgtVf; }

             break;
      case HARDNLO :
             xwgtSa=wgtSa()*jacob;
             xwgtSb=wgtSb()*jacob;
             xwgtR=wgtR()*jacob;
             xwgtSf=wgtSf()*jacob;
             if (!cc->cut(k1A,k2A))   { xwgt+= xwgtSa; }
             if (!cc->cut(k1B,k2B))   { xwgt+= xwgtSb; }
             if (!cc->cut(k11,k21,p1)) { xwgt+= xwgtR; }
             if (!cc->cut(k1f,k2f,pf)) { xwgt+= xwgtSf; }
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
             xwgtVf=wgtVf()*jacob;
             if (!cc->cut(k10,k20))   { xwgt+= xwgtLO+xwgtV+xwgtVplus; }
             if (!cc->cut(k1A,k2A))   { xwgt+= xwgtVa; }
             if (!cc->cut(k1B,k2B))   { xwgt+= xwgtVb; }
             if (!cc->cut(k1f,k2f,pf)) { xwgt+= xwgtVf; }
             break;
  };
  return wgtnorm()*xwgt;    
}



void qqNLO::fillBins(void)
{
  int n=0;
  double cts[5];
  double wgt[5];
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
             if (!cc->cut(k1f,k2f,pf)) { 
	       cts[n] = dist->Param(k1f,k2f);
               wgt[n] = wn*(xwgtSf+xwgtVf); 
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
             if (!cc->cut(k1f,k2f,pf)) { 
	       cts[n] = dist->Param(k1f,k2f);
               wgt[n] = wn*xwgtSf; 
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
             if (!cc->cut(k1f,k2f,pf)) { 
	       cts[n] = dist->Param(k1f,k2f);
               wgt[n] = wn*xwgtVf; 
               n++;
             }
             break;
  };
  dist->bin(n,cts,wgt);
}

