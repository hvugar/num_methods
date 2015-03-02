#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "minimum.h"
#include "gradient.h"
#include "methods.h"
#include "penalty.h"

double f(double *x, int n)
{
	return (x[0]-4.0)*(x[0]-4.0)+(x[1]-4.0)*(x[1]-4.0);
}

double g1(double *x, int n)
{
	return x[0]+x[1]-5.0;
}

int main(int argc, char** argv)
{
	double r = 100000000.00;
	
	int n = 2;
	int m = 0;
	int p = 1;
	double* x = (double*) malloc( sizeof(double) * n );
    x[0]    = +10;
    x[1]    = +10;
	
	RnFunction *h = (RnFunction*) malloc(sizeof(RnFunction*) * m);
	RnFunction *g = (RnFunction*) malloc(sizeof(RnFunction*) * p);
	g[0] = g1;
	
	penalty_method(f, x, n, h, m, g, p, r);
	
	free(x);
	free(h);
	free(g);
	
	return 0;
}