#include "integration.h"

extern "C" {
void gset_(double* ax, double* bx, long int* n, double* z, double* w);
}

double gauss(double (*f)(double), double a, double b, int n)
{
  double r=0.0;
  long int nn=n;
  double z[96], w[96];
  gset_(&a,&b,&nn,z,w);
  for (int i=0;i<n;i++) r+=f(z[i])*w[i];
  return r;
}

double digauss(double (*f)(double, double), double a1, double b1, 
               double (*ymin)(double), double (*ymax)(double), int n)
{
  double r=0.0;
  double  a2, b2;
  long int nn=n;
  double z1[96], w1[96], z2[96], w2[96];
  gset_(&a1,&b1,&nn,z1,w1);
  for (int i=0;i<n;i++) {
    a2=ymin(z1[i]);
    b2=ymax(z1[i]);
    gset_(&a2,&b2,&nn,z2,w2);
    for (int j=0;j<n;j++) r+=f(z1[i],z2[j])*w1[i]*w2[j];
  }
  return r;
}



double trapzd(double (*f)(double), double a, double b, int n)
{
  double x,tnm,sum,del;
  static double s;
  int it,j;

  if (n==0) {
      return (s=0.5*(b-a)*(f(a)+f(b)));
  } else {
    for (it=1,j=1;j<n;j++) it <<=1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for(sum=0.0,j=1;j<=it;j++,x+=del) sum += f(x);
    s=0.5*(s+(b-a)*sum/tnm);
    return s;
  }
}


const double EPS=1.0e-6;
const int JMAX=20;

double qromb(double (*f)(double), double a, double b)
{
  double s[JMAX], sold[JMAX], h[JMAX];
  int j, i;

  cout<<" Warning:  The routine 'qgromb' doesn't work so well!"<<endl;

  h[0]=1.0;
  s[0]=trapzd(f,a,b,0); 
  for (j=1;j<JMAX;j++) {
     h[j]=0.25*h[j-1];
     sold[0] = s[0];
     s[0] = trapzd(f,a,b,j);
     for (i=1;i<=j;i++) {
         sold[i] = s[i];
         s[i] = (h[j-i]*s[i-1]-h[j]*sold[i-1])/(h[j-i]-h[j]);
     }
     if (fabs(s[j]-s[j-1])<EPS*fabs(s[j])) return s[j];
  }   
  cerr <<"Too many steps in routine qromb.\n";
  exit(1);

  return 0.0;
}


double qsimp(double (*f)(double), double a, double b)
{
  double s,st,ost,os;

  ost = os = -9.23581e30;
  for (int j=0;j<JMAX;j++) {
     st=trapzd(f,a,b,j);
     s=(4.0*st-ost)/3.0;
     if (fabs(s-os) <EPS*fabs(os)) return s;
     os=s;
     ost = st;
  } 
  cerr <<"Too many steps in routine qsimp.\n";
  exit(1);

  return 0.0;
}

