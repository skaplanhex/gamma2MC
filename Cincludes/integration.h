#ifndef INTEGRATIONHDR
#define INTEGRATIONHDR

#include <math.h>
#include <iostream>
#include <stdlib.h>

using namespace std;

double gauss(double (*f)(double), double a, double b, int n=96);

double digauss(double (*f)(double, double), double a1, double b1, 
               double (*ymin)(double), double (*ymax)(double), int n=96);


// The following are from Numerical Recipes.  They don't work as well
// as the above routines.

double trapzd(double (*f)(double), double a, double b, int n);

double qromb(double (*f)(double), double a, double b);

double qsimp(double (*f)(double), double a, double b);

#endif
