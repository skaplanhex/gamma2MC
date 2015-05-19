#include "particle.h"

finalStates::finalStates( const int &size ) 
{
   nParticles = size;
   p = new particle[nParticles];
   for (int i=0; i<nParticles; i++)
      p[i] = particle(0.0);
}


finalStates::finalStates(const finalStates &fs) 
{
   nParticles=fs.nParticles;
   p=new particle[nParticles];
   for (int i=0;i<nParticles;i++)
      p[i]=fs.p[i];
}


finalStates& finalStates::operator=(const finalStates &fs) 
{
   if (nParticles!=fs.nParticles) {
      delete [] p;  
      nParticles=fs.nParticles;
      p=new particle[nParticles];
   }
   for (int i=0;i<nParticles;i++)
      p[i]=fs.p[i];
   return *this;
}


void finalStates::increment(particle &k)
{
   nParticles++;
   particle *ptemp=p;
   p=new particle[nParticles];
   for (int i=0;i<nParticles-1;i++) p[i]=ptemp[i];
   p[nParticles-1]=k;
   delete [] ptemp;
}


const int NCUT = 22;

double Li(const double &x)
{
   if (x<=0.5) {
     double L=0.0;
     double n=0.0;
     double r=1.0;
     for(int i=0;i<NCUT;i++) {
       n+=1.0;
       r*=x;
       if (r<=1.0e-10)
          return L;
       L+=r/n/n;
     }
     return L;
   } else {
     double xm = (1.0-x);          
     double L = PI*PI/6.0-log(x)*log(xm);
     double n = 0.0;
     double r = 1.0;
     for(int i=0;i<NCUT;i++) {
       n+=1.0;
       r*=xm;
       if (r<=1.0e-10) 
          return L;
       L-=r/n/n;
     }
     return L;
   }
}

double Li2(double xx)
{
   double x=xx;
   double L=0.0;
   double b=1.0;
   if (x<0.0) {
     L+=PI*PI/6.0-log(1.0-x)*log(-x);
     b=-1.0;
     x=1.0-x;
   }
   if (x>1.0) {
     double lx=log(x);
     L+=b*(PI*PI/3.0-0.5*lx*lx);
     b=-b;
     x=1.0/x;
   }
   return L+b*Li(x);
}

double Li2(double rr, double ttheta)
{
   double r=rr;
   double theta=ttheta;
   if (fabs(r-1.0)<1.0e-10) {
       return pow(PI-fabs(theta),2)/4.0 - PI*PI/12.0;
   }
   double L;
   double b;
   if (r>1.0) {
      r=1.0/r;
      double lr=log(r);
      theta=-theta;
      b=-1.0;
      L= PI*PI/3.0-lr*lr/2.0+theta*theta/2.0-PI*fabs(theta);
   } else {
      b=1.0;
      L=0.0;
   }
   double x=r*cos(theta);
   if (x>0.5) {
      x=1.0-x;
      double y=-r*sin(theta);
      double rold=r;
      double told=theta;
      r=sqrt(x*x+y*y);
      theta=atan2(y,x);
      L+=b*(PI*PI/6.0-log(r)*log(rold)+theta*told);
      b=-b;
   }
   if (r<0.53) {
     double n=0.0;
     double xr=1.0;
     for(int i=0;i<NCUT;i++) {
       n+=1.0;
       xr*=r;
       double zr=xr*cos(n*theta);
       if (fabs(zr)<1.0e-10)
          return L;
       L+=b*zr/n/n;
     }
     return L;
   } else {
     double rm = (1.0-r);          
     double sh=2.0*sin(theta/2.0);
     L += b*(pow(PI-fabs(theta),2)/4.0-PI*PI/12.0-log(r)*log(fabs(sh)));
     double Zn = 0.0;
     double k = 0.0;
     double xr=rm;
     double Sn=1.0;
     for(int i=0;i<NCUT;i++) {
        k+=1.0;
        xr*=rm;
        Sn/=sh;
        Zn+=Sn*cos( (theta-PI)*k/2.0 )/k;
        double zr=Zn*xr;          
        if (fabs(zr)<1.0e-10) 
            return L;
        L-=b*zr/(k+1.0);
     }
     return L;
   }
}

// The following are from Lance Dixon:

const double PISQ6=PI*PI/6.0;
const double ZETA3=1.202056903;

// this version uses 't Hooft and Veltman's change of variable
// good to ~ 10^(-16)

/*
double li2(double x){
 double x_0 = -0.30;
 double x_1 = 0.25;
 double x_2 = 0.51;
 if (x == 1) return PISQ6;
 if (x <= x_0){ 
   double temp = log(abs(1.0-x));
   return -li2(-x/(1.0-x)) - temp*temp/2 ; }
 else if (x < x_1){
   double z = - log(1.0-x);
   double temp = z*(1.0-z/4.0*(1.0-z/9.0*(1.0-z*z/100.0
                  *(1.0-5.0*z*z/294.0*(1.0-7.0*z*z/360.0
                  *(1.0-5.0*z*z/242.0*(1.0-7601.0*z*z/354900.0
                  *(1.0-91.0*z*z/4146.0*(1.0-3617.0*z*z/161840.0)
                   ))))))));
   return temp; }
   else if (x < x_2) return - li2(-x) + li2(x*x)/2.0 ;
   else { return PISQ6 - li2(1.0-x) 
                  - log(abs(x))*log(abs(1.0-x)) ; }
}
*/

double Li3(double x){
 double x_0 = -1.0;
 double x_1 = -0.85;
 double x_2 = 0.25;
 double x_3 = 0.63;
 double x_4 =  1.0;
 if (x == 1) return ZETA3;
 if (x <= x_0){ 
   double lnx = log(-x);
   return Li3(1.0/x) - PISQ6*lnx - lnx*lnx*lnx/6.0; }
 else if (x < x_1){
   return Li3(x*x)/4.0 - Li3(-x); }
   else if (x < x_2){
     double z = - log(1.0-x);
     double temp = z*(1.0-3.0*z/8.0*(1.0-17.0*z/81.0*(1.0-15*z/136.0
                    *(1.0-28.0*z/1875.0*(1.0+5.0*z/8.0*(1.0-304.0*z/7203.0
                    *(1.0+945.0*z/2432.0*(1.0-44.0*z/675.0*(1.0+7.0*z/24.0
                    *(1.0-26104.0*z/307461.0*(1.0+1925.0*z/8023.0
                    *(1.0-53598548.0*z/524808375.0
                    *(1.0+22232925.0*z/107197096.0
                     )))))))))))));
     return temp; }
     else if (x < x_3){
       return Li3(x*x)/4.0 - Li3(-x); }
       else if (x < x_4){
         double ln1x = log(1.0-x); 
         return -Li3(1.0-x) - Li3(-x/(1.0-x)) + ZETA3 + PISQ6*ln1x
           - log(x)*ln1x*ln1x/2.0 + ln1x*ln1x*ln1x/6.0; }
       else { 
         double lnx = log(x);
         return Li3(1/x) + 2.0*PISQ6*lnx - lnx*lnx*lnx/6.0; }
}

// the tetralog, Li4.
// good to ~ 1.7 * 10^(-12) (worst areas:  x = 0.9501, also x = - 0.9701)

double Li4(double x){
 double x_0 = -1.0;
 double x_1 = -0.97;
 double x_2 = 0.25;
 double x_3 = 0.95;
 double x_4 =  1.0;
 if (x == -1) return -0.35 * PISQ6*PISQ6;
 if (x == 1) return 0.4 * PISQ6*PISQ6;
 if (x <= x_0){ 
   double lnx = log(-x);
   return - Li4(1./x) - 0.5 * PISQ6*lnx*lnx
          - 1./24. * lnx*lnx*lnx*lnx - 0.7 * PISQ6*PISQ6; }
 else if (x < x_1){
   return Li4(x*x)/8. - Li4(-x); }
   else if (x < x_2){
     double z = - log(1.-x);
     double temp = z*(1.-7.*z/16.*(1.-151.*z/567.*(1.0-411.*z/2416.
                    *(1.-24986.*z/256875.*(1.-805.*z/49972.
                    *(1.+583406.*z/1159683.*(1.-7455.*z/137272.
                    *(1.+659921.*z/2444175.*(1.-251559.*z/2639684.
                    *(1.+24259894.*z/136410197.*(1.-30625595.*z/218339046.
                    *(1.+2134239258113.*z/16698772722450.
                    *(1.-1640443805715.*z/8536957032452.
                     )))))))))))));
     return temp; }
     else if (x < x_3){
       return Li4(x*x)/8. - Li4(-x); }
       else if (x < x_4){
        double y = 1.-x; 
        double lny = log(y);
        return 0.4*PISQ6*PISQ6 + ZETA3 * log(1.-y) + 0.5*PISQ6 * y*y 
         + (-11./36.+0.5*PISQ6+1./6.*lny) * y*y*y
         + (11./24.*PISQ6+1./4.*lny-19./48.) * y*y*y*y
         + (7./24.*lny+5./12.*PISQ6-599./1440.) * y*y*y*y*y
         + (137./360.*PISQ6+5./16.*lny-79./192.) * y*y*y*y*y*y
         + (7./20.*PISQ6-3343./8400.+29./90.*lny) * y*y*y*y*y*y*y
          + (363./1120*PISQ6-21977./57600.+469./1440.*lny) * y*y*y*y*y*y*y*y;
                        }
       else { 
         double lnx = log(x);
         return - Li4(1/x) + PISQ6*lnx*lnx 
                - 1./24. * lnx*lnx*lnx*lnx + 0.8 * PISQ6*PISQ6; }
}

