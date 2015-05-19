#ifndef VEGASGHDR
#define VEGASGHDR

#include <math.h>
#include <stdio.h>
#include "random.h"
#include "average.h"
#ifdef __GNUG__
#define PI M_PI
#else
#define PI _PI
#endif
#define NDMX 50 // If you change this, make sure to change
                // it in vegas.C also.
#define MXDIM 10

class vegasGrid 
/*
   Uses VEGAS to make a grid for event generation.
   Importance Sampling Only!
*/
{
  protected:
        Random rr;
	int init;        // input flag
	int ndim;        // number of dimensions to be integrated
	double lowlim[MXDIM+1];  // lower limits array
	double uplim[MXDIM+1];   // upper limits array
	int itmx;        // number of iterations
	unsigned long ncall; // approximate number of calls per integration
	int nprn;        // flag which controls diagnostic output
	double tgral;   // Best estimate of integral
	double sd;      // Standard deviation
	double chi2a;   // Chi-squared per degree of freedom
	void rebin(void);
	int i,it,j,k,nd,ndo,ng,npg;
	int ia[MXDIM+1];
	double calls,dv2g,dxg,f,f2,f2b,fb,rc,ti,tsi,wgt,xjac,xn,xnd,xo;
        double volume;
	double vwgt;  // weight for one call.
	double d[NDMX+1][MXDIM+1];    
	double di[NDMX+1][MXDIM+1];   
	double dt[MXDIM+1];
	double dx[MXDIM+1];
	double r[NDMX+1];
	double x[MXDIM+1];
	double xi[MXDIM+1][NDMX+1];  
	double xin[NDMX+1];
	double schi,si,swgt;
	void mainLoop(void);    // Main iteration loop.  
	void rebin(double rrc);
	vegasGrid(int np=0) {nprn=np;}
	virtual ~vegasGrid(void) {}
	
  public:
	virtual double fxn(void)=0;
	virtual void processEvent(void)=0;
	double integral(void) {return tgral;}
	double stanDev(void) {return sd;}
	double chiSqr(void) {return chi2a;}
	void fullPrint(void) {nprn=1;}
	void basicPrint(void) {nprn=0;}
	void noPrint(void) {nprn=-1;}
        void seed(long idum) {rr.seed(idum);}
                     //  The seed for the random number generator
	             //  must be a NEGATIVE integer.
	             //  The default is -1.

	void limits(double ll[], double ul[], int nd) {
					    	ndim=nd;
					    	for (int i=1;i<=ndim;i++) {
					    		lowlim[i]=ll[i];
					    		uplim[i]=ul[i];
					    	}
					    }
	void integrate(int in, unsigned long nc, int it ) {
		                    init=in; ncall=nc; itmx=it;
		                    mainLoop();
				  }

	void multiEvents(unsigned long nc );  // Get nc events with 
                                              // preset grid.

	void singleEvent(void); // Get one event with preset grid.
};

#undef MXDIM
#undef NDMX

#endif
