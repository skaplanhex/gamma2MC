#include "vegasGrid.h"
#define ALPH 1.5
#define NDMX 50 // If you change this, make sure to change
                // it in vegasGrid.h also.
#define TINY 1.0e-30

void vegasGrid::mainLoop(void)
{
	if (init <= 0) {
		ndo=1;  
		for (j=1;j<=ndim;j++) xi[j][1]=1.0;
	}
	if (init <= 1) si=swgt=schi=0.0;
	if (init <= 2) {
		nd=NDMX;
		ng=1;
	 	npg=ncall;
		calls=(double)npg;
		dv2g=1.0/(npg-1.0);
		xnd=nd;
		dxg = xnd;
		xjac=1.0/calls;
                volume=1.0;
		for (j=1;j<=ndim;j++) {
			dx[j]=uplim[j]-lowlim[j];
			xjac *= dx[j];
			volume *= dx[j];
		}
		if (nd != ndo) {     // Do binning if necessary.
	        int ndd;
            if (nd > ndo) {
		       ndd=nd;
  		    } else {
		       ndd=ndo;
            }
		    for (i=1;i<=ndd;i++) r[i]=1.0;
		    for (j=1;j<=ndim;j++) rebin(ndo/xnd);
		    ndo=nd;
		}
		if (nprn >= 0) {
		  cout<<" Input parameters for vegas:  ndim = "<<ndim<<"  ncall = "<<calls<<endl;
		  cout<<"                              it = "<<it<<"    itmx = "<<itmx<<endl;
		  cout<<"                              nprn = "<<nprn<<"  ALPH = "<<ALPH<<endl;
		  cout<<"                              ndo = "<<ndo<<"  nd = "<<nd<<endl;
//			printf("%s:  ndim= %3d  ncall= %8.0f\n",
//				" Input parameters for vegas",ndim,calls);
//			printf("%28s  it=%5d  itmx=%5d\n"," ",it,itmx);
//			printf("%28s  nprn=%3d  ALPH=%5.2f\n"," ",nprn,ALPH);
//			printf("%28s  ndo=%4d  nd=%4d\n"," ",ndo,nd);
			for (j=1;j<=ndim;j++) {
			  cout<<"                              xl["<<j<<"] =  "<<lowlim[j]<<"    xu["<<j<<"] =  "<<uplim[j]<<endl;
//				printf("%30s xl[%2d]= %11.4g xu[%2d]= %11.4g\n",
//					" ",j,lowlim[j],j,uplim[j]);
			}
		}
	}
	for (it=1;it<=itmx;it++) {  // Main loop. Can enter here to do additional
	                  // itmx iterations with all other parameters unchanged.
		ti=tsi=0.0;
		for (j=1;j<=ndim;j++) {
			for (i=1;i<=nd;i++) d[i][j]=di[i][j]=0.0;
		}
		fb=f2b=0.0;
		for (k=1;k<=npg;k++) {
			wgt=xjac;
			for (j=1;j<=ndim;j++) {
				xn=(rr.uniform())*dxg+1.0; // Change to (rr.uniform())
				ia[j]=(int)(xn);
				if (ia[j]>1) {
					xo=xi[j][ia[j]]-xi[j][ia[j]-1];
					rc=xi[j][ia[j]-1]+(xn-ia[j])*xo;
				} else {
					ia[j]=1;
					xo=xi[j][ia[j]];
					rc=(xn-ia[j])*xo;
				}
				x[j]=lowlim[j]+rc*dx[j];
				wgt *= xo*xnd;
			}
			f=wgt*fxn();
			f2=f*f;
			fb += f;
			f2b += f2;
			for (j=1;j<=ndim;j++) {
				di[ia[j]][j] += f;
				d[ia[j]][j] += f2;
			}
		}
		f2b=sqrt(f2b*npg);
		f2b=(f2b-fb)*(f2b+fb);
		if (f2b <= 0.0) f2b=TINY;
		ti = fb;
		tsi = f2b;
		tsi *= dv2g;          // Compute final results for this iteration.
		wgt=1.0/tsi;
		si += wgt*ti;
		schi += wgt*ti*ti;
		swgt += wgt;
		tgral=si/swgt;
		chi2a=(schi-si*tgral)/(it-0.9999);
		if (chi2a < 0.0) chi2a = 0.0;
		sd=sqrt(1.0/swgt);
		tsi=sqrt(tsi);
		if (nprn >= 0) {
		  cout<<" iteration no.   "<<it<<": integral = "<<setprecision(7)<<setw(14)<<ti<<" +/-    "<<setprecision(2)<<setw(9)<<tsi<<endl;
//			printf("%s %3d : integral = %14.7g +/-  %9.2g\n",
//				" iteration no.",it,ti,tsi);
		  cout<<" all iterations: integral = "<<setprecision(7)<<setw(14)<<tgral<<" +/- "<<setprecision(2)<<setw(9)<<sd<<" chi**2/IT n= "<<setw(5)<<chi2a<<setprecision(7)<<endl;
//			printf("%s integral =%14.7g+/-%9.2g chi**2/IT n = %9.2g\n",
//				" all iterations:  ",tgral,sd,chi2a);
			if (nprn) {
				for (j=1;j<=ndim;j++) {
				  cout<<" DATA FOR axis "<<j<<endl;
//					printf(" DATA FOR axis  %2d\n",j);
				  cout<<"     X      delta i          X      delta i          X      delta i"<<endl;    
//					printf("%6s%13s%11s%13s%11s%13s\n",
//						"X","delta i","X","delta i","X","delta i");
					for (i=1+nprn/2;i<=nd;i += nprn+2) {
					  cout<<setprecision(5)<<setw(8)<<xi[j][i]<<setprecision(4)<<setw(12)<<di[i][j];
					  cout<<setprecision(5)<<setw(12)<<xi[j][i+1]<<setprecision(4)<<setw(12)<<di[i+1][j];
					  cout<<setprecision(5)<<setw(12)<<xi[j][i+2]<<setprecision(4)<<setw(12)<<di[i+2][j];
					  cout<<setprecision(7)<<endl;
//						printf("%8.5f%12.4g%12.5f%12.4g%12.5f%12.4g\n",
//							xi[j][i],di[i][j],xi[j][i+1],
//							di[i+1][j],xi[j][i+2],di[i+2][j]);
					}
				}
			}
		}
		for (j=1;j<=ndim;j++) {    // Refine the grid.
			xo=d[1][j];
			xn=d[2][j];
			d[1][j]=(xo+xn)/2.0;
			dt[j]=d[1][j];
			for (i=2;i<nd;i++) {
				rc=xo+xn;
				xo=xn;
				xn=d[i+1][j];
				d[i][j] = (rc+xn)/3.0;
				dt[j] += d[i][j];
			}
			d[nd][j]=(xo+xn)/2.0;
			dt[j] += d[nd][j];
		}
		for (j=1;j<=ndim;j++) {
			rc=0.0;
			for (i=1;i<=nd;i++) {
				if (d[i][j] < TINY) d[i][j]=TINY;
				r[i]=pow((1.0-d[i][j]/dt[j])/
				    (log(dt[j])-log(d[i][j])),ALPH);
				rc += r[i];
			}
			rebin(rc/xnd);
		}
	}
}


void vegasGrid::multiEvents(unsigned long nc)
{
        ncall=nc;
	nd=NDMX;
	dxg = nd;
	if (nd != ndo) {
	  cerr<<"Error!  Must initialize grid first!"<<endl;
	}
	for (k=1;k<=ncall;k++) {
		vwgt=volume;
		for (j=1;j<=ndim;j++) {
			xn=(rr.uniform())*dxg+1.0; // Change to (rr.uniform())
			ia[j]=(int)(xn);
			if (ia[j]>1) {
				xo=xi[j][ia[j]]-xi[j][ia[j]-1];
				rc=xi[j][ia[j]-1]+(xn-ia[j])*xo;
			} else {
				ia[j]=1;
				xo=xi[j][ia[j]];
				rc=(xn-ia[j])*xo;
			}
			x[j]=lowlim[j]+rc*dx[j];
			vwgt *= xo*xnd;
		}
		processEvent();
	}
}


void vegasGrid::singleEvent(void)
{
	wgt=volume;
	for (j=1;j<=ndim;j++) {
		xn=(rr.uniform())*xnd+1.0;
		ia[j]=(int)(xn);
		if (ia[j]>1) {
			xo=xi[j][ia[j]]-xi[j][ia[j]-1];
			rc=xi[j][ia[j]-1]+(xn-ia[j])*xo;
		} else {
			ia[j]=1;
			xo=xi[j][ia[j]];
			rc=(xn-ia[j])*xo;
		}
		x[j]=lowlim[j]+rc*dx[j];
		wgt *= xo*xnd;
	}
	vwgt=wgt;
}


void vegasGrid::rebin(double rrc)
// Utility routine to rebin the vector of densities xi[j] into new bins defined 
// by the vector r.
{
	k=0;
	double dr=0.0;
	xn=0.0;
	
	for (i=1;i<nd;i++) {
		while (rrc > dr) {
			dr += r[++k];
			xo=xn;
			xn=xi[j][k];
		}
		dr -= rrc;
		xin[i]=xn-(xn-xo)*dr/r[k];
	}
	for (i=1;i<nd;i++) xi[j][i]=xin[i];
	xi[j][nd]=1.0;
}

#undef ALPH
#undef NDMX
#undef TINY
