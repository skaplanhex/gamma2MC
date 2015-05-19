#ifndef AVERAGEHDR
#define AVERAGEHDR


#include <math.h>
#include "tdio.h"


class average {
 private:
     int calls;
     double sum, sumSq; 
 public:
     average(void) {calls=0; sum=sumSq=0.0;}
     void add(double wgt=1.0) {
                               calls+=1;
                               sum+=wgt;
                               sumSq+=wgt*wgt;
			      }
     double mean (void) {
                         if (calls<=1) return sum;
                         else return sum/((double) calls);
		        }
     double stanDev (void) {
                            if (calls<=1) return sum;
                            else {
                               double dcalls = (double) calls;
                               double mean = sum/dcalls;
                               return sqrt(fabs(sumSq/dcalls-mean*mean)
                                           /(dcalls-1.0));
			    }
			   }

     void reset(void) {calls=0; sum=sumSq=0.0;}
};

#endif







