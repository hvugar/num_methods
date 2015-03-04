#include <stdio.h>
#include "minimum.h"
#include "gradient.h"
#include "methods.h"
#include "penalty.h"

extern double f_rosenbrock(double *x, int n); 
extern double f(double *x, int n); 
extern double h1(double *x, int n); 

void penalty_sample1()
{
	double r = 1.00;
	
	int n = 2;
	int m = 1;
	int p = 0;
	double* x = (double*) malloc( sizeof(double) * n );
    x[0]    = +10;
    x[1]    = +10;
	
	RnFunction *h = (RnFunction*) malloc(sizeof(RnFunction*) * m);
	RnFunction *g = (RnFunction*) malloc(sizeof(RnFunction*) * p);
	h[0] = h1;
	
	penalty_method(f, x, n, h, m, g, p, r);
	
	free(x);
	free(h);
	free(g);
}
