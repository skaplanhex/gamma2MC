#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include "HiggsNLO.h"
#include "ggNLO.h"
#include "qqNLO.h"
#include "Cuts.h"
#define MXDIM 50
int ndim = 7;
int init =0;
double ll[MXDIM+1] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double ul[MXDIM+1] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

enum process_type {HIGGS=0, GG=1, QQ=2};
// event_types are {LO, NLO, EASYNLO, HARDNLO}
//enum qq_type { QQTOTAL = 0, QQONLY = 1, QGONLY = 2 };
enum distribution_type { COSTHETASTAR=0, YSTAR=1, TANHYSTAR=2, YGG=3, PHIGG=4, QT=5, MGG=6, USERDEFINED=7 };
//enum cut_type {none=0, standard=1, ptcut=2, hcut=3, frixione=4,
//               noIsolation=5, jetveto=6, annulus=7, D0=8, userdefined=9};
enum norm_type {scaled, unscaled};	           
int main()
{
  gamma2NLO *mc;

  ifstream infile ("gamma2MC.indat");
  ofstream outfile;
  streambuf *psbuf, *backup;
  string inputstring;
  // Read Name of Output file.  Re-direct all output to it.
  getline(infile, inputstring,' ');
  const char *outfilename=0;
  outfilename=inputstring.data();
  outfile.open(outfilename);
  backup=cout.rdbuf();
  psbuf=outfile.rdbuf();
  cout.rdbuf(psbuf);
  infile.ignore(1024,'\n');
  // Read Collider Type
  collider coll;
  getline(infile, inputstring,' ');
  if (inputstring=="pp") {coll=pp;}
  else if (inputstring=="ppbar") {coll=ppbar;}
  else { 
       cerr<<" Invalid Collider type. "<<endl;
       return 0;
  };
  cout<<"Collider type: " <<inputstring<<endl;
  infile.ignore(1024,'\n');
  // Read Collider Energy
  double energy;
  infile>>energy;
  cout<<"COM Energy: "<<energy<<" GeV"<<endl;
  infile.ignore(1024,'\n');
  // Read Parton distributions
  long int iset;
  getline(infile, inputstring,' ');
  infile.ignore(1024,'\n');
  infile>>iset;
  cout<<"PDF set name: "<<inputstring<<" with member number="<<iset<<endl;
  infile.ignore(1024,'\n');
  pdfinit(inputstring,iset);
  // Read Process Type
  process_type proc;
  getline(infile, inputstring,' ');
  if (inputstring=="GG") {proc=GG; mc = new ggNLO(energy,coll);
  cout<<"Process: gg --> gamma gamma X"<<endl;}
  else if (inputstring=="HIGGS") {proc=HIGGS; mc = new HiggsNLO(energy,coll); 
  cout<<"Process: gg --> H X,  H --> gamma gamma"<<endl;}
  else if (inputstring=="QQ") {proc=QQ; mc = new qqNLO(energy,coll);
  cout<<"Process: qq --> gamma gamma X"<<endl;}
  else { 
       cerr<<" Invalid Process type. "<<endl;
       return 0;
  };
  infile.ignore(1024,'\n');
  //
  if (proc==QQ) {
     qq_type qq;
     getline(infile, inputstring,' ');
     if (inputstring=="QQTOTAL") {qq=QQTOTAL; 
     cout<<"    QQ+QG+GQ included"<<endl;}
     else if (inputstring=="QQONLY") {qq=QQONLY;
     cout<<"    QQ only included"<<endl;}
     else if (inputstring=="QQTOTAL") {qq=QGONLY;
     cout<<"    QG+GQ only included"<<endl;}
     else { 
       cerr<<" Invalid QQ type. "<<endl;
       return 0;
     };
     mc->setQQType(qq);
  };
  infile.ignore(1024,'\n');
  // Read Order of Calculation
  event_type evt;
  getline(infile, inputstring,' ');
  if (inputstring=="LO") {evt=LO;}
  else if (inputstring=="NLO") {evt=NLO;}
  else if (inputstring=="EASYNLO") {evt=EASYNLO;}
  else if (inputstring=="HARDNLO") {evt=HARDNLO;}
  else { 
       cerr<<" Invalid Order of Calculation. "<<endl;
       return 0;
  };
  cout<<"Order of Calculation: "<<inputstring<<endl;
  mc->setEventType(evt);
  infile.ignore(1024,'\n');
  // Read Cut types
  cutStandard cs;
  cutFrixione cf;
  cutD0 cd;
  cutNone cn;
  cutHiggs ch;
  cutPt cp;
  cutPhoton cni;
  cutJet cj;
  cutAnnulus ca;
  cutUser cu;
  getline(infile, inputstring,' ');
  infile.ignore(1024,'\n');
  if (inputstring=="none") {
    infile.ignore(1024,'\n');
    mc->setCuts(&cn);
    cout<<"No cuts. "<<endl;
  } else if (inputstring=="standard") {
    double ycut, pt2, pt1;
    double Rcut, Etcut;
    infile>>ycut>>pt2>>pt1;
    infile.ignore(1024,'\n');
    infile>>Rcut>>Etcut;
    cs.setPhotonCuts(ycut,pt2,pt1);
    cs.setIsolation(Rcut,Etcut);
    mc->setCuts(&cs);
    cout<<"Standard cuts with Standard isolation "<<endl;
    cout<<"Photon Cuts: |y_gam| < "<<ycut<<" , pt2 > "<<pt2 
        <<"GeV, pt1 > "<<pt1<<" GeV"<<endl;
    cout<<"Isolation Cuts: R = "<<Rcut<<" , Et = "<<Etcut<<" GeV "<<endl;
  } else if (inputstring=="ptcut") {
    double ycut, pt2, pt1;
    double Rcut, Etcut, ptg;
    infile>>ycut>>pt2>>pt1;
    infile.ignore(1024,'\n');
    infile>>Rcut>>Etcut>>ptg;
    cp.setPhotonCuts(ycut,pt2,pt1);
    cp.setIsolation(Rcut,Etcut);
    cp.setPtCut(ptg);
    mc->setCuts(&cp);
    cout<<"Standard cuts with Standard isolation plus observed parton pt cut "
	<<endl;
    cout<<"Photon Cuts: |y_gam| < "<<ycut<<" , pt2 > "<<pt2 
        <<"GeV, pt1 > "<<pt1<<" GeV"<<endl;
    cout<<"Isolation Cuts: R = "<<Rcut<<" , Et = "<<Etcut<<" GeV "<<endl;
    cout<<"Parton Pt Cut: pt_parton > "<<ptg<<" GeV "<<endl;
  } else if (inputstring=="hcut") {
    double yh;
    infile.ignore(1024,'\n');
    infile>>yh;
    ch.setYcut(yh);
    mc->setCuts(&ch);
    cout<<"Simple rapidity cut on the gamma-gamma (Higgs) system "<<endl;
    cout<<"|y_gamgam| < "<<yh<<endl;
  } else if (inputstring=="frixione") {
    double ycut, pt2, pt1;
    double Rcut, epsilon;
    infile>>ycut>>pt2>>pt1;
    infile.ignore(1024,'\n');
    infile>>Rcut>>epsilon;
    cf.setPhotonCuts(ycut,pt2,pt1);
    cf.setIsolation(Rcut,epsilon);
    mc->setCuts(&cf);
    cout<<"Standard cuts with Frixione isolation "<<endl;
    cout<<"Photon Cuts: |y_gam| < "<<ycut<<" , pt2 > "<<pt2 
        <<"GeV, pt1 > "<<pt1<<" GeV"<<endl;
    cout<<"Isolation Cuts: R = "<<Rcut<<" , epsilon = "<<epsilon<<endl;
  } else if (inputstring=="noIsolation") {
    double ycut, pt2, pt1;
    infile>>ycut>>pt2>>pt1;
    infile.ignore(1024,'\n');
    cni.setPhotonCuts(ycut,pt2,pt1);
    mc->setCuts(&cni);
    cout<<"Standard cuts with no isolation "<<endl;
    cout<<"Photon Cuts: |y_gam| < "<<ycut<<" , pt2 > "<<pt2 
        <<"GeV, pt1 > "<<pt1<<" GeV"<<endl;
  } else if (inputstring=="jetveto") {
    double ycut, pt2, pt1;
    double Rcut, Etcut;
    double Rjet, Etjet;
    infile>>ycut>>pt2>>pt1;
    infile.ignore(1024,'\n');
    infile>>Rcut>>Etcut>>Rjet>>Etjet;
    cj.setPhotonCuts(ycut,pt2,pt1);
    cj.setIsolation(Rcut,Etcut);
    cj.setJetVeto(Rjet, Etjet);
    mc->setCuts(&cj);
    cout<<"Standard cuts with Standard isolation and Jet Veto "
	<<endl;
    cout<<"Photon Cuts: |y_gam| < "<<ycut<<" , pt2 > "<<pt2 
        <<"GeV, pt1 > "<<pt1<<" GeV"<<endl;
    cout<<"Isolation Cuts: R = "<<Rcut<<" , Et = "<<Etcut<<" GeV "<<endl;
    cout<<"Jet Veto Cut: Rj = "<<Rjet<<" , Etj = "<<Etjet<<" GeV "<<endl;
 } else if (inputstring=="annulus") {
    double ycut, pt2, pt1;
    double Rcut, Etcut;
    double Rjet, Etjet;
    infile>>ycut>>pt2>>pt1;
    infile.ignore(1024,'\n');
    infile>>Rcut>>Etcut>>Rjet>>Etjet;
    ca.setPhotonCuts(ycut,pt2,pt1);
    ca.setIsolation(Rcut,Etcut);
    ca.setJetVeto(Rjet, Etjet);
    mc->setCuts(&ca);
    cout<<"Standard cuts with Standard isolation and an observed jet in annulus around a photon "<<endl;
    cout<<"Photon Cuts: |y_gam| < "<<ycut<<" , pt2 > "<<pt2 
        <<"GeV, pt1 > "<<pt1<<" GeV"<<endl;
    cout<<"Isolation Cuts: R = "<<Rcut<<" , Et = "<<Etcut<<" GeV "<<endl;
    cout<<"Jet Annulus: R1 = "<<Rcut<<" , R2 = "<<Rjet<<endl;
    cout<<"Minimum jet energy: Etj = "<<Etjet<<" GeV "<<endl;
 } else if (inputstring=="D0") {
    double ycut, pt2, pt1;
    double Rcut, epsilon;
    infile>>ycut>>pt2>>pt1;
    infile.ignore(1024,'\n');
    infile>>Rcut>>epsilon;
    cd.setPhotonCuts(ycut,pt2,pt1);
    cd.setIsolation(Rcut,epsilon);
    mc->setCuts(&cd);
    cout<<"Standard cuts with D0 isolation "<<endl;
    cout<<"Photon Cuts: |y_gam| < "<<ycut<<" , pt2 > "<<pt2 
        <<"GeV, pt1 > "<<pt1<<" GeV"<<endl;
    cout<<"Isolation Cuts: R = "<<Rcut<<" , epsilon = "<<epsilon<<endl;
 } else if (inputstring=="userdefined") {
    double ycut, pt2, pt1;
    double P1, P2, P3, P4;
    infile>>ycut>>pt2>>pt1;
    infile.ignore(1024,'\n');
    infile>>P1>>P2>>P3>>P4;
    cu.setPhotonCuts(ycut,pt2,pt1);
    cu.setIsolation(P1,P2,P3,P4);
    mc->setCuts(&cu);
    cout<<"Standard Cuts with User-defined Isolation "<<endl;
    cout<<"Photon Cuts: |y_gam| < "<<ycut<<" , pt2 > "<<pt2 
        <<"GeV, pt1 > "<<pt1<<" GeV"<<endl;
    cout<<"Isolation Cut Parameters: P1 = "<<P1<<" , P2 = "<<P2<<" , P3 = "<<P3<<" , P4 = "<<P4<<endl;
 } else { 
       cerr<<" Invalid Cut Type. "<<endl;
       return 0;
  };
  infile.ignore(1024,'\n');
  // Read upper and lower values of M_gamgam
  double Qlow, Qhigh, Q, Qbin;
  infile>>Qlow>>Qhigh;
  if (proc==HIGGS) {
    Q=Qlow;
    Qbin=0.0;
    cout<<"Higgs Mass = M_gamgam = "<<Q<<" GeV"<<endl;
  } else {
    Q=(Qlow+Qhigh)/2.0;
    Qbin=Qhigh-Qlow;
    cout<<"Diphoton invariant mass : "<<Qlow<<" < M_gamgam < "<<Qhigh<<" GeV"<<endl;
  };
  mc->newQ(Q);
  mc->newQbin(Qbin);
  infile.ignore(1024,'\n');
  // Read Scales
  double chiR, chiF, chifr;
  infile>>chiR>>chiF>>chifr;
  cout<<"alpha_s( "<<chiR*Q<<" GeV ) = "<<alphas(chiR*Q)<<endl; 
  cout<<"chiR = "<<chiR<<endl;
  cout<<"chiF = "<<chiF<<endl;
  cout<<"chifr = "<<chifr<<endl;
  mc->newChiF(chiF);
  mc->newChiR(chiR);
  mc->newChifr(chifr);
  infile.ignore(1024,'\n');
  // Read Grid Adaptation variables
  int itmx;
  unsigned long ncallAdapt;
  infile>>itmx>>ncallAdapt;
  cout<<"Grid Variables : itmx = "<<itmx<<", ncallAdapt = "<<ncallAdapt<<endl;
  infile.ignore(1024,'\n');
  // Read Final number of Calls
  unsigned long ncallFinal;
  infile>>ncallFinal;
  cout<<"Final number of calls = "<<ncallFinal<<endl;
  infile.ignore(1024,'\n');
  // Read Vegas Seed
  long iseed;
  infile>>iseed;
  cout<<"Vegas Seed = "<<iseed<<endl;
  infile.ignore(1024,'\n');
  // Read Distribution types
  int nbins;
  double hmin, hmax;
  distribution_type dtype;
  getline(infile, inputstring,' ');
  infile.ignore(1024,'\n');
  infile>>nbins>>hmin>>hmax;
  CosThetaStar d1(nbins,hmin,hmax);
  YStar d2(nbins,hmin,hmax);
  TanhYStar d3(nbins,hmin,hmax);
  Ygg d4(nbins,hmin,hmax);
  Phigg d5(nbins,hmin,hmax);
  qT d6(nbins,hmin,hmax);
  Mgamgam d7(nbins,hmin,hmax);
  distUser d8(nbins,hmin,hmax);     
  if (inputstring=="COSTHETASTAR") {
    dtype=COSTHETASTAR;
    mc->setDistribution(&d1);
  } else if (inputstring=="YSTAR") {
    dtype=YSTAR;
    mc->setDistribution(&d2);
  } else if (inputstring=="TANHYSTAR") {
    dtype=TANHYSTAR;
    mc->setDistribution(&d3);
  } else if (inputstring=="YGG") {
    dtype=YGG;
    mc->setDistribution(&d4);
  } else if (inputstring=="PHIGG") {
    dtype=PHIGG;
    mc->setDistribution(&d5);
  } else if (inputstring=="QT") {
    dtype=QT;
    mc->setDistribution(&d6);
  } else if (inputstring=="MGG") {
    dtype=MGG;
    mc->setDistribution(&d7);
  } else if (inputstring=="USERDEFINED") {
    dtype=USERDEFINED;
    mc->setDistribution(&d8);
 } else { 
       cerr<<" Invalid Distribution Type. "<<endl;
       return 0;
  };
  infile.ignore(1024,'\n');
  // Read Normalization type
  norm_type dnorm;
  getline(infile, inputstring,' ');
  if (inputstring=="scaled") {
    dnorm=scaled;
  } else if (inputstring=="unscaled") {
    dnorm=unscaled;
  } else { 
       cerr<<" Scaling of Distribution is undetermined. "<<endl;
       return 0;
  };
  infile.ignore(1024,'\n');
  infile.close();
  mc->limits(ll,ul,ndim);
  mc->seed(iseed);
  cout<<endl;
  mc->integrate(init,ncallAdapt,itmx);
  mc->multiEvents(ncallFinal);
  cout<<endl;
  cout<<"cross section = "<<mc->xsec.mean()<<" +/- "<<setprecision(2)<<mc->xsec.stanDev()<<" fb"<<setprecision(7)<<endl<<endl;

  switch(dtype) {
  case COSTHETASTAR :
     cout<<"Cos(thetaStar) Distribution:";
     cout<<" ( nbins = "<<nbins<<" , hmin = "<<hmin<<" , hmax = "<<hmax<<" )"<<endl<<endl;
     cout<<"Cos(thetsStar) : d sigma/d Cos(thetaStar)(fb): Error Estimate (fb)"<<endl;
     break;
  case YSTAR :
     cout<<"ystar Distribution:";
     cout<<" ( nbins = "<<nbins<<" , hmin = "<<hmin<<" , hmax = "<<hmax<<" )"<<endl<<endl;
     cout<<"ystar          : d sigma/d ystar (fb)        : Error Estimate (fb)"<<endl;
     break;
  case TANHYSTAR :
     cout<<"tanh(ystar) Distribution:";
     cout<<" ( nbins = "<<nbins<<" , hmin = "<<hmin<<" , hmax = "<<hmax<<" )"<<endl<<endl;
     cout<<"tanh(ystar)    : d sigma/d tanh(ystar) (fb)  : Error Estimate (fb)"<<endl;
     break;
  case YGG :
     cout<<"y_{gamma gamma} Distribution:";
     cout<<" ( nbins = "<<nbins<<" , hmin = "<<hmin<<" , hmax = "<<hmax<<" )"<<endl<<endl;
     cout<<"y_gamgam       : d sigma/d y_gamgam (fb)     : Error Estimate (fb)"<<endl;
     break;
  case PHIGG :
     cout<<"phi_{gamma gamma} Distribution:";
     cout<<" ( nbins = "<<nbins<<" , hmin = "<<hmin<<" , hmax = "<<hmax<<" )"<<endl<<endl;
     cout<<"phi_gamgam     : d sigma/d phi_gamgam (fb)   : Error Estimate (fb)"<<endl;
     break;
  case QT :
     cout<<"qT Distribution:";
     cout<<" ( nbins = "<<nbins<<" , hmin = "<<hmin<<" , hmax = "<<hmax<<" )"<<endl<<endl;
     cout<<"qT (GeV)       : d sigma/d qT (fb/GeV)       : Error Estimate (fb/GeV)"<<endl;
     break;
  case MGG :
     cout<<"M_{gamma gamma} Distribution:";
     cout<<" ( nbins = "<<nbins<<" , hmin = "<<hmin<<" , hmax = "<<hmax<<" )"<<endl<<endl;
     cout<<"M_gamgam (GeV) : d sigma/d M_gamgam (fb/GeV) : Error Estimate (fb/GeV)"<<endl;
     break;
  case USERDEFINED :
     cout<<"User-Defined Distribution:";
     cout<<" ( nbins = "<<nbins<<" , hmin = "<<hmin<<" , hmax = "<<hmax<<" )"<<endl<<endl;
     cout<<"X (unit)       : d sigma/d X (fb/unit)       : Error Estimate (fb/unit)"<<endl;
     break;
  };
  double hnorm=1.0;
  switch(dnorm) {
  case scaled :
    hnorm=mc->xsec.mean();
    cout<<"               : *(1/sigma)                  :"<<endl;
    break;
  case unscaled :
    break;
  };
  mc->printDistribution(hnorm);

  cout.rdbuf(backup);
  outfile.close();

  return 0;
}


