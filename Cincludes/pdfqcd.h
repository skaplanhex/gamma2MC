#ifndef PDFQCDHDR
#define PDFQCDHDR

/*

Updated January 26, 2011.  Now uses LHAPDF's.

*/

#include <math.h>
#include <iostream>
#include <string>
#include "particle.h"
using namespace std;

enum collider { pp = 0, ppbar = 1 } ;


void pdfinit(const string &pdft, long int iset);
           // Initialize PDF's, where pdft is the PDFset name 
           // (eq MRST2001nlo.LHpdf) ("lhapdf-config" script finds
           // the correct path) and iset is the PDF member number.
           // (best fit set is always iset=0)

double pdf(long int iparton, double x, double Q);
           // Call PDF's.
           // -6 tbar
           // -5 bbar
           // -4 cbar
           // -3 sbar
           // -2 ubar
           // -1 dbar
           //  0 g
           //  1 d
           //  2 u
           //  3 s
           //  4 c
           //  5 b
           //  6 t



void pdfpairs(collider coll, double x1, double x2, double Q, double &fgg, 
        double &fgq, double &fqg, double &fqiqj, double &fqiqi, double &fqiai);

//  NOTE!!! pdfpairs is modified from earlier version.  fgq and fqg
//  are now separated.  In pdfqcd.h they were added together and called
//  simply fgq.
//  8/2/01:  Second correction.  Original pdfpairs was for ppbar.
//  Now separated into a pp and a ppbar choice.


double pdfquarks(double x, double Q);

//  Just the sum over all quark AND antiquark pdfs.

double pdfeff(double x, double Q);
           // This is the effective BFKL pdf.

double alphas(double mu);
           // Returns appropriate alpha_s, corresponding to PDF set.

//  The following are legacy:

double alphasLO(double mu);
           // LO running with alphasMZ=.1273 and mb=4.5 GeV.

double nflavors(double mu);
           // 4 flavors if mu<mb=4.5 GeV, else 5.  

#endif

