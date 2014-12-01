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


int main(int argc, char** argv)
{
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
    double *xc = (double*) malloc( sizeof(double) * N );

    double* gr = (double*) malloc(sizeof(double) * N);
    gradient(f, X, N, delta_x, gr);

	int i;
    for (i=0; i<N; i++)
    {
        xc[i] = X[i] - alpha * gr[i];
    }

    double result = f(xc, N);

    free(gr);
    free(xc);

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
    double *xc = (double*) malloc( sizeof(double) * N );

    double* gr = (double*) malloc(sizeof(double) * N);
    gradient(f1, X, N, delta_x, gr);
	
	int i;
    for (i=0; i<N; i++)
    {
        xc[i] = X[i] - alpha * gr[i];
    }
	
    double result = f1(xc, N);

    free(gr);
    free(xc);

    return result;
}

double f2(const double x)
{
	return x*x*x*x - 3*x*x*x + 2*x*x + 2;
}