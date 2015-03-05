#include <stdio.h>
#include "penalty.h"

extern double f_rosenbrock(double *x, int n); 

void penalty_sample1()
{
	double f(double *x, int n)  { return (x[0]-4.0)*(x[0]-4.0)+(x[1]-4.0)*(x[1]-4.0); }
	double h1(double *x, int n) { return x[0]+x[1]-5.0; }

	int n = 2;
	int m = 1;
	int p = 0;
	double* x = (double*) malloc( sizeof(double) * n );
    x[0]    = +10;
    x[1]    = +10;
	
	RnFunction *h = (RnFunction*) malloc(sizeof(RnFunction*) * m);
	h[0] = h1;
		
	RnFunction *g = (RnFunction*) malloc(sizeof(RnFunction*) * p);
	
	double r1 = 1.00;
	double r2 = 1.00;
	double epsilon = 0.01;
	penalty_method(f, x, n, h, m, g, p, r1, r2, epsilon);
	
	free(x);
	free(h);
	free(g);
}

void penalty_sample2()
{
	double f(double *x, int n) { return x[0]*x[0] + x[1]*x[1] + 2*x[1]; }

	double h1(double *x, int n) { return x[0]*x[0] + x[1]*x[1] - 1.0; }
	double g1(double *x, int n) { return x[0] + 2*x[1] - 0.5; }
	double g2(double *x, int n) { return x[0]; }
	double g3(double *x, int n) { return x[1]; }
	
	int n = 2;
	int m = 1;
	int p = 3;
	double* x = (double*) malloc( sizeof(double) * n );
    x[0]    = +0.5;
    x[1]    = +0.5;
	
	RnFunction *h = (RnFunction*) malloc(sizeof(RnFunction*) * m);
	h[0] = h1;
	
	RnFunction *g = (RnFunction*) malloc(sizeof(RnFunction*) * p);
	g[0] = g1;
	g[1] = g2;
	g[2] = g3;
	
	double r1 = 1.00;
	double r2 = 1.00;
	double epsilon = 0.1;
	penalty_method(f, x, n, h, m, g, p, r1, r2, epsilon);
	
	free(x);
	free(h);
	free(g);
}

void penalty_sample3()
{
	double f(double *x, int n) { return x[0]*x[0]; }
	
	double g1(double *x, int n) { return x[0]-1.0; }
	double g2(double *x, int n) { return 2.0-x[0]; }
	
	int n = 1;
	int m = 0;
	int p = 2;
	double* x = (double*) malloc( sizeof(double) * n );
    x[0]    = +1.5;
	
	RnFunction *h = (RnFunction*) malloc(sizeof(RnFunction*) * m);
	
	RnFunction *g = (RnFunction*) malloc(sizeof(RnFunction*) * p);
	g[0] = g1;
	g[1] = g2;

	double r1 = 1.00;
	double r2 = 1.00;
	double epsilon = 0.05;
	penalty_method(f, x, n, h, m, g, p, r1, r2, epsilon);
	
	free(x);
	free(h);
	free(g);

}
