#include "Distribution.h"


double YStar::Param(particle p1, particle p2)
{
return fabs((p1.rapidity()-p2.rapidity())/2.0);
}

double TanhYStar::Param(particle p1, particle p2)
{
  return tanh(fabs((p1.rapidity()-p2.rapidity())/2.0));
}

double CosThetaStar::Param(particle p1, particle p2)
{

  particle pdad=p1+p2;
  particle k1=p1;
  k1.invBoostDirect(pdad);
  return fabs(k1.cosTheta());
}

double Ygg::Param(particle p1, particle p2)
{

  particle pdad=p1+p2;
  return fabs(pdad.rapidity());
}

double Phigg::Param(particle p1, particle p2)
{
  double dphi=fabs(p1.phi()-p2.phi());
  if (dphi>PI) {dphi=2.0*PI-dphi;}
  return dphi;
}

double qT::Param(particle p1, particle p2)
{
  particle q=p1+p2;
  return q.pt();
}

double Mgamgam::Param(particle p1, particle p2)
{
  particle q=p1+p2;
  return q.getMass();
}

double distUser::Param(particle p1, particle p2)
{
  particle q=p1+p2;
  return q.getMass();
}
