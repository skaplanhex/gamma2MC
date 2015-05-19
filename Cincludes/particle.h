#ifndef PARTICLEHDR
#define PARTICLEHDR

#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <string>

using namespace std;

/*
  Defines a four vector momenta class, with all the expected 
  algebraic properties, as well as a particle class which is
  essentially just a four momentum with a well-defined mass.
  Also defines a finalState class which is just an indefinite
  array of particles.
*/


class fourMom {
 protected:
    double px; double py; double pz; double pe;
    inline friend ostream& operator<<(ostream &os, const fourMom &k);
 public:
    fourMom( double xx=0.0, double yy=0.0, double zz = 0.0, 
               double ee=0.0 ) {px=xx; py=yy; pz=zz; pe=ee;}
    fourMom( const fourMom &k ) {
                  px=k.px; py=k.py; pz=k.pz; pe=k.pe;}
    fourMom operator=(const fourMom &k) {
                  px=k.px; py=k.py; pz=k.pz; pe=k.pe; return *this;}
    fourMom operator+(const fourMom &k ) {
                  fourMom q; q.px=px+k.px; q.py=py+k.py; 
                  q.pz=pz+k.pz; q.pe=pe+k.pe; return q;}
    fourMom operator-(const fourMom &k ) {
                  fourMom q; q.px=px-k.px; q.py=py-k.py; 
                  q.pz=pz-k.pz; q.pe=pe-k.pe; return q;}
    double operator*(const fourMom &p ) {
                  double w; w=pe*p.pe-px*p.px-py*p.py-pz*p.pz;
                  return w;}
    fourMom operator+=(const fourMom &k) {
                  px+=k.px; py+=k.py; pz+=k.pz; pe+=k.pe; return *this;}
    fourMom operator-=(const fourMom &k) {
                  px-=k.px; py-=k.py; pz-=k.pz; pe-=k.pe; return *this;}
    fourMom operator/(const double &a) {
                  fourMom q; q.px=px/a; q.py=py/a; q.pz=pz/a;
                  q.pe=pe/a; return q;}  
    friend fourMom operator-(const fourMom &p ) {
                  fourMom q; q.px=-p.px; q.py=-p.py; q.pz=-p.pz;
                  q.pe=-p.pe;  return q;}  
    friend fourMom operator+(const fourMom &p ) {
                  fourMom q; q.px=p.px; q.py=p.py; q.pz=p.pz;
                  q.pe=p.pe;  return q;}  
    friend fourMom operator*(const double &a, const fourMom &p) {
                  fourMom q; q.px=a*p.px; q.py=a*p.py; q.pz=a*p.pz;
                  q.pe=a*p.pe; return q;}  
    fourMom operator*(const double &a) {
                  fourMom q; q.px=a*px; q.py=a*py; q.pz=a*pz;
                  q.pe=a*pe; return q;}  
    double squared(void) const {
                  double r; r=pe*pe-px*px-py*py-pz*pz;
                  return r;}
    void boostZ(const double &beta) {
                  double gamma = 1.0/sqrt(1.0-beta*beta);
                  double zz=pz; double ee=pe; pe=gamma*(ee+beta*zz); 
                  pz=gamma*(beta*ee+zz);}
    void rotateZX(const double &phi) {
                  double zz=pz; double xx=px; pz=zz*cos(phi) - xx*sin(phi);
                  px= zz*sin(phi) + xx*cos(phi);}
    void rotateXY(const double &phi) {
                  double xx=px; double yy=py; px=xx*cos(phi) - yy*sin(phi);
                  py= xx*sin(phi) + yy*cos(phi);}
    void boost(const double &beta, const double &theta, const double &phi) { 
                  boostZ(beta);
                  rotateZX(theta);
                  rotateXY(phi); }
    void invBoost(const double &beta, const double &theta, const double &phi) 
                 { 
                  rotateXY(-phi); 
                  rotateZX(-theta);
                  boostZ(-beta); }    
    double gamma(void) const {return pe/sqrt(squared());}
    double beta(void) const {return p()/pe;}
    double p(void) const {return sqrt(px*px+py*py+pz*pz);}
    double pt(void) const {return sqrt(px*px+py*py);}
    double rapidity(void) const {return 0.5*log((pe+pz)/(pe-pz));}
    double e(void) const {return pe;}
    double x(void) const {return px;}
    double y(void) const {return py;}
    double z(void) const {return pz;}
    double cosTheta(void) const {
      double pp=p();
      if (pp!=0.0) {return pz/pp;}
       else {return 0.0;}}
    double phi(void) const { return atan2(py,px);}
};


class particle : public fourMom {
 private:
    double mass;
    double qe;  // electric charge
    string nType; // particle name  
    inline friend ostream& operator<<(ostream &os, const particle &k);
 public:
    particle( double m=0.0, double qqe=0.0,
              string n="" ) : fourMom(0.0,0.0,0.0,m) {
              mass = m; qe=qqe; nType=n;}
    particle(const particle &part) : 
              fourMom(part.px,part.py,part.pz,part.pe) {
              mass = part.mass; qe=part.qe; 
              nType=part.nType;}
    particle( fourMom k) : fourMom(k) {
               double msq=k.squared();
               if  (msq>=0.0) 
                    mass = sqrt(msq);
               else if (msq>-1.0e-8) 
                    mass = 0.0;
               else {
                    cerr<<"Attempt to create a particle with\n";
                    cerr<<"negative mass-squared.\n";
                    cerr<<"msq = "<<msq<<"\n";
                    cerr<<"e,x,y,z "<<k.e()<<" "<<k.x()<<" "<<k.y()
                                     <<" "<<k.z()<<"\n";
                    exit(1);}
               qe=0.0; nType="";
               }         

    double charge(void) {return qe;}
    string name(void) {return nType;}
    double getMass(void) {return mass;}

    inline void setMom(const double &energy, const double &cosTheta=1.0,
           const double &phi=0.0);
    inline void setMomYPtPhi(const double &y, const double &pt,
           const double &phi=0.0);
    particle operator=(const particle &part) {
                    px=part.px; py=part.py; pz=part.pz; pe=part.pe;
                    mass=part.mass; qe=part.qe;
                    nType=part.nType; return *this;} 
    void boost(const double &beta, const double &theta, const double &phi)  
                    {fourMom::boost(beta,theta,phi);}
    void boost(const particle &pdad) {
                  double beta = sqrt(1.0-pdad.squared()/pdad.e()/pdad.e());
                  double theta = acos(pdad.cosTheta());
                  double phi = pdad.phi();
                  fourMom::boost(beta,theta,phi);}
    void boostDirect(const particle &pdad) {
                  double beta = sqrt(1.0-pdad.squared()/pdad.e()/pdad.e());
                  double theta = acos(pdad.cosTheta());
                  double phi = pdad.phi();
                  fourMom::boost(beta,theta,phi);
                  beta=0.0;
                  fourMom::invBoost(beta,theta,phi);}
    // boostDirect can be realized as a single boost along the
    // momentum direction of pdad.  
    void invBoost(const double &beta, const double &theta, const double &phi) 
                    {fourMom::invBoost(beta,theta,phi);} 
    void invBoost(const particle &pdad) {
                  double beta = sqrt(1.0-pdad.squared()/pdad.e()/pdad.e());
                  double theta = acos(pdad.cosTheta());
                  double phi = pdad.phi();
                  fourMom::invBoost(beta,theta,phi);}

    void invBoostDirect(const particle &pdad) {
                  double beta = sqrt(1.0-pdad.squared()/pdad.e()/pdad.e());
                  double theta = acos(pdad.cosTheta());
                  double phi = pdad.phi();
                  fourMom::invBoost(beta,theta,phi);
                  beta=0.0;
                  fourMom::boost(beta,theta,phi);}
    // invBoostDirect can be realized as a single boost against the
    // momentum direction of pdad.


};


class finalStates {
 protected:
   int nParticles;  // number of particles produced in the event
                    // ( = 6 at tree level)

   particle *p;    //  the array of final state particles

 public:
   finalStates(void) {nParticles=0; p=0;}
   finalStates( const int &size ) ;
   finalStates( const finalStates &fs);
   ~finalStates(void) {delete [] p;}
   finalStates& operator=(const finalStates &k) ;
   void increment(particle &k);
   void rangeCheck(int index) const {
              if (index<0 || index>=nParticles) { 
                    cerr<<"Attempt to index a non-existent particle.\n";
                    cerr<<"Only "<<nParticles<<
                          " particles in the final state.\n";
                    exit(1);
	      } }
   particle& operator[](int index) const { 
              rangeCheck(index);
              return p[index];}
   int number(void) const {return nParticles;}
   void boostZ(const double &beta) {
                  for (int i=0;i<nParticles;i++) p[i].boostZ(beta);}
   void boost(const double &beta, const double &theta, const double &phi) { 
                  for (int i=0;i<nParticles;i++) p[i].boost(beta,theta,phi);}
   void invBoost(const double &beta, const double &theta, const double &phi) {
                  for (int i=0;i<nParticles;i++) p[i].invBoost(beta,theta,phi);}

};


inline void particle::setMom(const double &energy, const double &cosTheta,
           const double &phi)
{
   double p;
   double sinTheta=sqrt(1.0-cosTheta*cosTheta);
   pe=energy;
   p = sqrt(energy*energy-mass*mass);
   px=p*sinTheta*cos(phi);
   py=p*sinTheta*sin(phi);
   pz=p*cosTheta;
}


inline void particle::setMomYPtPhi(const double &y, const double &pt,
           const double &phi)
{
   double mt = sqrt(mass*mass+pt*pt);
   pe=mt*cosh(y);
   px=pt*cos(phi);
   py=pt*sin(phi);
   pz=mt*sinh(y);
}


inline ostream& operator<<(ostream &os , const particle &k) {
          os.setf(ios::fixed,ios::floatfield);
          os<<k.nType<<setprecision(3)<<setw(10)<<k.pe
                     <<setw(10)<<k.px<<setw(10)<<k.py
                     <<setw(10)<<k.pz;
          os.setf(ios::scientific,ios::floatfield);
          os<<setprecision(5);
          return os;
	}


inline ostream& operator<<(ostream &os , const fourMom &k) {
          os.setf(ios::fixed,ios::floatfield);
          os<<setprecision(3)<<setw(10)<<k.pe
                     <<setw(10)<<k.px<<setw(10)<<k.py
                     <<setw(10)<<k.pz;
          os.setf(ios::scientific,ios::floatfield);
          os<<setprecision(5);
          return os;
	}



// Basic definitions:

#define PI M_PI
#define SQRT2 M_SQRT2

// The standard model parameters:  
  
const double alpha = 1.0/137.0356;
const double alphaMZ = 1.0/127.934;
const double SSqtw = 0.23124;
const double CSqtw = 0.76876;
const double StwCtw = sqrt(SSqtw*CSqtw);
const double hbarCsqPB = 0.389379292e9;  // in GeV^2 pb.
const double hbarCsqNB = 0.389379292e6;  // in GeV^2 nb.

// up quark weak couplings:
const double QU = 2.0/3.0;
const double QUL = (0.5-QU*SSqtw)/StwCtw;
const double QUR = -QU*SSqtw/StwCtw;
const double QUV = 0.5*(QUR+QUL);
const double QUA = 0.5*(QUR-QUL);

// down quark weak couplings;
const double QD = -1.0/3.0;
const double QDL = (-0.5-QD*SSqtw)/StwCtw;
const double QDR = -QD*SSqtw/StwCtw;
const double QDV = 0.5*(QDR+QDL);
const double QDA = 0.5*(QDR-QDL);

// electron weak couplings:
const double QE=-1.0;
const double QEL= (-0.5+SSqtw)/StwCtw;
const double QER= SSqtw/StwCtw;
const double QEV = 0.5*(QER+QEL);
const double QEA = 0.5*(QER-QEL);

// neutrino weak couplings:
const double QN= 0.0;
const double QNL= 0.5/StwCtw;
const double QNR= 0.0;
const double QNV = 0.5*(QNR+QNL);
const double QNA = 0.5*(QNR-QNL);

// lepton masses:
const double ME=0.5110e-3;
const double MMU=0.1057;
const double MTAU=1.777;
const double MNU=0.0;

// quark masses:
const double MU=0.0;
const double MD=0.0;
const double MS=0.2;
const double MC=1.4;
const double MB=4.5;
const double MT=174.3;

// Nucleon masses:
const double MP=0.938272;
const double MN=0.939566;

// boson masses;
const double MZ = 91.187;
const double MZ2= MZ*MZ;
const double MW = 80.41;
const double MW2= MW*MW;

const double GF=1.16639e-5;
//const double vev=sqrt(MW2*SSqtw/PI/alphaMZ);
const double vev=1.0/sqrt(SQRT2*GF);

// generic particles:

const particle eMinus(ME,QE,    " e-       ");
const particle muMinus(MMU,QE,  " mu-      ");
const particle tauMinus(MTAU,QE," tau-     ");
const particle nuE(MNU,QN,      " nuE      ");
const particle nuMu(MNU,QN,     " nuMu     ");
const particle nuTau(MNU,QN,    " nuTau    ");

const particle ePlus(ME,-QE,    " e+       ");
const particle muPlus(MMU,-QE,  " mu+      ");
const particle tauPlus(MTAU,-QE," tau+     ");
const particle nubarE(MNU,-QN,  " nubarE   ");
const particle nubarMu(MNU,-QN, " nubarMu  ");
const particle nubarTau(MNU,-QN," nubarTau ");

const particle up(MU,QU,        " u        ");
const particle charm(MC,QU,     " c        ");
const particle top(MT,QU,       " t        ");
const particle down(MD,QD,      " d        ");
const particle strange(MS,QD,   " s        ");
const particle bottom(MB,QD,    " b        ");

const particle ubar(MU,-QU,     " ubar     ");
const particle cbar(MC,-QU,     " cbar     ");
const particle tbar(MT,-QU,     " tbar     ");
const particle dbar(MD,-QD,     " dbar     ");
const particle sbar(MS,-QD,     " sbar     ");
const particle bbar(MB,-QD,     " bbar     ");

const particle proton(MP,1.0,   " P        ");
const particle neutron(MN,0.0,  " N        ");
const particle pbar(MP,-1.0,    " Pbar     ");
const particle nbar(MN,0.0,     " Nbar     ");

const particle photon(0.0,0.0,  " photon   ");
const particle Zzero(MZ,0.0,    " Z        ");
const particle Wplus(MW,1.0,    " W+       ");
const particle gluon(0.0,0.0,   " gluon    ");
const particle Wminus(MW,-1.0,  " W-       ");

double Li(const double &x);

// Li is the dilogarithm.  It assumes 0 < x < 1.

double Li2(double r, double theta);

// Li2 is the Real part of the dilog function of
// the complex variable Z = r * exp(I*theta), 
// with r>0 and -PI < theta < PI.  
// Li2 = - Real( int(0..Z) log(1.0-z)/z ).
//
// (If theta=0 and r < 1, it reduces to Li(r).)

double Li2(double x);

// This is the same as Li2(r, theta) with (r=x,theta=0) for x>0 an
// (r=-x,theta=PI) for x<0.  (Note that this only gives the Real part.)

double Li3(double x);

double Li4(double x);

#endif

