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
bool cutPhoton::inZoneOne(const particle& p)
{
  double y = fabs( p.rapidity() );
  bool inside = (y > y1absmin) && (y < y1absmax);
  return inside;
}
bool cutPhoton::inZoneTwo(particle p)
{
  double y = fabs( p.rapidity() );
  bool inside = (y > y2absmin) && (y < y2absmax);
  return inside;
}
int cutPhoton::cut(particle p1, particle p2)     
{  
  // return p1.rapidity()<-ycut||p1.rapidity()>ycut
  //                          ||p2.rapidity()<-ycut||p2.rapidity()>ycut
  //                          ||p1.pt()<ptcut2||p2.pt()<ptcut2
  //                          ||(p1.pt()<ptcut1&&p2.pt()<ptcut1);

  // assuming we have access to four cuts: y1absmin, y1absmax, y2absmin, y2absmax
  // we want y1 to be within y1absmin and y1absmax and for y2 to be within y2absmin and y2absmax
  bool p1InZoneOne = inZoneOne(p1);
  bool p1InZoneTwo = inZoneTwo(p1);
  bool p2InZoneOne = inZoneOne(p2);
  bool p2InZoneTwo = inZoneTwo(p2);

  // if etaCheck is true, then eta cuts are satisfied
  bool etaCheck = ( p1InZoneOne && p2InZoneTwo ) || ( p1InZoneTwo || p2InZoneOne );

  // figure out which particle has the highest pt
  double pt1 = p1.pt();
  double pt2 = p2.pt();

  double lowpt,highpt;
  if (pt1 < pt2){
    lowpt = pt1;
    highpt = pt2;
  }
  else if (pt2 < pt1){
    lowpt = pt2;
    highpt = pt1;
  }
  else {
    lowpt = pt2;
    highpt = lowpt;
  }
  bool ptCheck = (highpt > ptcut2) && (lowpt > ptcut2) && (highpt > ptcut1);

  // event is good if both pt and eta cuts are satisfied
  bool eventGood = etaCheck && ptCheck;

  if (eventGood) // don't cut if both eta and pt cuts are met
    return 0;
  else // if one of the cuts aren't met
    return 1;
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
  double ycut = 1.e10; // a placeholder so the code doesn't break, shouldn't affect anything.
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

