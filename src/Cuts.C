#include "Cuts.h"

double deltaR(particle p1, particle p2);

int cutNone::cut(particle p1, particle p2) 
{
  return 0;
}

int cutNone::cut(particle p1, particle p2, particle p3) 
{
  return 0;
}

int cutPhoton::cut(particle p1, particle p2)     
{  
  return p1.rapidity()<-ycut||p1.rapidity()>ycut
                           ||p2.rapidity()<-ycut||p2.rapidity()>ycut
                           ||p1.pt()<ptcut2||p2.pt()<ptcut2
                           ||(p1.pt()<ptcut1&&p2.pt()<ptcut1);
}

int cutPhoton::cut(particle p1, particle p2, particle p3)
{  
  return (cutPhoton::cut(p1,p2)||cutIsolation(p1,p2,p3));
}

int cutHiggs::cut(particle p1, particle p2)     
{  
  particle p=p1+p2;
  return p.rapidity()<-yHcut||p.rapidity()>yHcut;
}

int cutHiggs::cut(particle p1, particle p2, particle p3)
{  
  return cutHiggs::cut(p1,p2);
}

int cutPt::cut(particle p1, particle p2)     
{  
  return 1;
}

int cutPt::cut(particle p1, particle p2, particle p3)
{  
  //  return (cutPhoton::cut(p1,p2,p3)||(p3.pt()<ptgcut));
  return (cutPhoton::cut(p1,p2,p3)||(p3.pt()<ptgcut)||(p3.rapidity()<-ycut)||
      (p3.rapidity()>ycut));
}

int cutStandard::cutIsolation(particle p1, particle p2, particle p3)
{
  return ((deltaR(p1,p3)<Rcut||deltaR(p2,p3)<Rcut)&&(p3.pt()>Etcut));
}
	           
int cutPhoton::cutIsolation(particle p1, particle p2, particle p3)
{
  return 0;
}
	           
int cutFrixione::cutIsolation(particle p1, particle p2, particle p3)
{
  double R1=deltaR(p1,p3);
  double R2=deltaR(p2,p3);
  return ( ((R1<Rcut)&&(p3.pt()/p1.pt()>epsilon*(1.0-cos(R1))/Ccut)) ||
           ((R2<Rcut)&&(p3.pt()/p2.pt()>epsilon*(1.0-cos(R2))/Ccut)) );
}
	 
int cutD0::cutIsolation(particle p1, particle p2, particle p3)
{
  double R1=deltaR(p1,p3);
  double R2=deltaR(p2,p3);
  return ( ((R1<Rcut)&&(p3.pt()/p1.pt()>epsilon)) ||
           ((R2<Rcut)&&(p3.pt()/p2.pt()>epsilon)) );
}
	 
int cutJet::cutIsolation(particle p1, particle p2, particle p3)
{
  return (cutStandard::cutIsolation(p1,p2,p3)||
((deltaR(p1,p3)<Rjet||deltaR(p2,p3)<Rjet)&&(p3.pt()>Etjet)));
}
          	           
int cutAnnulus::cut(particle p1, particle p2)     
{  
  return 1;
}

int cutAnnulus::cut(particle p1, particle p2, particle p3)
{  
  return (cutPhoton::cut(p1,p2)||(p3.pt()<Etjet)||
       !(   (deltaR(p1,p3)<Rjet&&deltaR(p1,p3)>Rcut)||
         (deltaR(p2,p3)<Rjet&&deltaR(p2,p3)>Rcut)  )  );
}

int cutUser::cutIsolation(particle p1, particle p2, particle p3)
{
  double Rcut=P1;
  double Etcut=P2;
  return ((deltaR(p1,p3)<Rcut||deltaR(p2,p3)<Rcut)&&(p3.pt()>Etcut));
}
	           
double deltaR(particle p1, particle p2)
{
  double dy=p1.rapidity()-p2.rapidity();
  double dphi=fabs(p1.phi()-p2.phi());
  if (dphi>PI) {dphi=2.0*PI-dphi;}
  return sqrt(dy*dy+dphi*dphi);
}

