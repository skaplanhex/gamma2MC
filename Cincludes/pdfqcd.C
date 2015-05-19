#include <string>
#include "pdfqcd.h"
#include "LHAPDF/LHAPDF.h"
using namespace std;

const double mb=4.5;

void pdfinit(const string &pdf_name,long int iset) 
{
  LHAPDF::setVerbosity(LHAPDF::SILENT);
  LHAPDF::initLHAPDF();
  LHAPDF::initPDFSetByName(pdf_name);
  LHAPDF::initPDF(iset);
}

double pdf(long int iparton, double x, double Q)
{
  return LHAPDF::xfx(x,Q,iparton)/x;
}

double pdfeff(double x, double Q)
{
    double af=9.0/4.0;
    double pd=0.0;
    for (long int iparton=-5;iparton<6;iparton++) {
        if (iparton==0) {pd+= af*pdf(iparton,x,Q);
        } else {pd+=pdf(iparton,x,Q);}
    }
    return pd/af;
}

void pdfpairs(collider col, double x1, double x2, double Q, double &fgg, 
    double &fgq, double &fqg, double &fqiqj, double &fqiqi, double &fqiai)
{
    long int iparton, aparton;
    double fq1, fq2, fa1, fa2, fg1, fg2;
    double fqt1=0.0;
    double fqt2=0.0;
    fqiqi=0.0;
    fqiai=0.0;
    for (iparton=1;iparton<6;iparton++) {
        aparton=-iparton;
        fq1=pdf(iparton,x1,Q);
        fa1=pdf(aparton,x1,Q);
        if (col==pp) {
           fa2=pdf(aparton,x2,Q);  // hadron 2 is p.
           fq2=pdf(iparton,x2,Q);  // hadron 2 is p.
        } else {
           fa2=pdf(iparton,x2,Q);  // hadron 2 is pbar.
           fq2=pdf(aparton,x2,Q);  // hadron 2 is pbar.
        }
    fqt1+=fq1+fa1;
    fqt2+=fq2+fa2;
    fqiqi+=fq1*fq2+fa1*fa2;
    fqiai+=fq1*fa2+fa1*fq2;
    }
    fqiqj=fqt1*fqt2 - fqiqi - fqiai;
    iparton=0;
    fg1=pdf(iparton,x1,Q);
    fg2=pdf(iparton,x2,Q);
    fgg = fg1*fg2;
    fgq = fg1*fqt2;
    fqg = fqt1*fg2;
}

double pdfquarks(double x, double Q)
{
    long int iparton, aparton;
    double fq=0.0;
    for (iparton=1;iparton<6;iparton++) {
        aparton=-iparton;
        fq+= pdf(iparton,x,Q)
            +pdf(aparton,x,Q);
    }
    return fq;
}

double alphas(double mu)
{
return LHAPDF::alphasPDF(mu);
}

const double alphasMZ = 0.1273;
double alphasLO(double QQ) {
  double AA=1.0/alphasMZ + (23.0/6.0/PI)*log(QQ/MZ);
  if ( QQ<=mb ) {AA = AA + log(QQ/mb)/(3.0*PI);}
  if ( AA<=0.0  ) {  
       cerr<<" Value of Q given was "<<QQ<<"\n";
       cerr<<" You have reached the Landau singularity.\n";
       cerr<<" All hope is lost!\n";
       exit(1);
       }
  return 1.0/AA;
}

double nflavors(double QQ) {
  double nfl;
  if ( QQ<=mb ) {nfl=4.0;}
  else {nfl=5.0;}
  return nfl;
}

