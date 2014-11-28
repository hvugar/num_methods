#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "minimum.h"
#include "gradient.h"
#include "methods.h"

double f(double *x, int n);
double g(double alpha);
double f1(double *x, int n);
double g1(double alpha);

double f2(const double x);

double  *X;
int     N;
double  epsilon;
double  delta_x;

int j=0;
int k=0;

int main(int argc, char** argv)
{
//	double a,b;
//	double x0 = -0.2;
//	double m = straight_line_search_metod(f2, x0, 0.1, &a, &b);
//	printf("a: %f, b: %f m: %f\n", a, b, m);
//	double x = golden_section_search_min(f2, a, b, 0.001);
//	printf("a: %f, b: %f x: %f\n", a, b, x);
//	printf("f(a): %f, f(b): %f f(x): %f\n", f2(a), f2(b), f2(x));

//	return 0;
	
    epsilon = 0.0001;
    delta_x = 0.00001;
    
	N = 2;
    X  = (double*) malloc( sizeof(double) * N );

    X[0]    = -1.2;
    X[1]    = +1.0;
	
//	fast_proximal_gradient_method(f, g, X, N, delta_x, epsilon);
	conjugate_gradient_method(f1, g1, X, N, delta_x, epsilon);
	
	free(X);
	
	return 0;
}

double f(double *x, int n)
{
    return pow(x[0],3) + 2*pow(x[1],2) - 3*x[0] - 4*x[1];
}

double g(double alpha)
{
    double* _x = (double*) malloc( sizeof(double) * N );

    double* gr = (double*) malloc(sizeof(double) * N);
    gradient(f, X, N, delta_x, gr);

	int i;
    for (i=0; i<N; i++)
    {
        _x[i] = X[i] - alpha * gr[i];
    }

    double result = f(_x, N);

    free(gr);
    free(_x);

    return result;
}

double f1(double *x, int n)
{
	double _x = x[0];
	double _y = x[1];
    return ((1-_x)*(1-_x)) + 100*(_y-_x*_x)*(_y-_x*_x);
}

double g1(double alpha)
{
    double* _x = (double*) malloc( sizeof(double) * N );

    double* gr = (double*) malloc(sizeof(double) * N);
    gradient(f1, X, N, delta_x, gr);
	
	int i;
    for (i=0; i<N; i++)
    {
        _x[i] = X[i] - alpha * gr[i];
    }
	
    double result = f1(_x, N);

    free(gr);
    free(_x);

    return result;
}

double f2(const double x)
{
	return x*x*x*x - 3*x*x*x + 2*x*x + 2;
}