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

double f_rosenbrock(double *x, int n)
{
	double x1 = x[0];
	double x2 = x[1];
    return ((1-x1)*(1-x1)) + 100*(x2-x1*x1)*(x2-x1*x1);
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
	
/*	
    double epsilon	= 0.001;		//dovrun sona catma meyari
	double grad_eps	= 0.005;		//gradient
	double line_eps	= 0.1;			//parcani bolme
	double gold_eps	= 0.0001;		//qizil qayda ucun
    
	int n = 2;
    double* x  = (double*) malloc( sizeof(double) * n );
    
	x[0]    = -1.2;
    x[1]    = +1.0;
	conjugate_gradient_method(f_rosenbrock, x, n, line_eps, gold_eps, grad_eps, epsilon);
	
	puts("");
	
    x[0]    = -1.2;
    x[1]    = +1.0;
	conjugate_gradient_method2(f_rosenbrock, x, n, line_eps, gold_eps, grad_eps, epsilon);

	free(x);
*/	
	return 0;
}