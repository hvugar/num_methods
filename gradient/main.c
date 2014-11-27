#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "minimum.h"
#include "gradient.h"
#include "methods.h"

double f(double *x, int n);
double g(double alpha);
double g1(double alpha);

double  *x;
int     N;
double  epsilon;
double  delta_x;

double **S;
double **Xj;
int j=0;
int k=0;

int main(int argc, char** argv)
{
    epsilon = 0.001;
    delta_x = 0.0001;
    
	N = 2;
    x  = (double*) malloc( sizeof(double) * N );

    x[0]    = -0.5;
    x[1]    = -1.0;
	
	fast_proximal_gradient_method(f, g, x, N, delta_x, epsilon);

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
    gradient(f, x, N, delta_x, gr);

	int i;
    for (i=0; i<N; i++)
    {
        _x[i] = x[i] - alpha * gr[i];
    }

    double result = f(_x, N);

    free(gr);
    free(_x);

    return result;
}

double g1(double alpha)
{
    double* gr = (double*) malloc( sizeof(double) * N );
    gradient(f, x, N, delta_x, gr);

    double* _x = (double*) malloc( sizeof(double) * N );
	
	int i;
    for (i=0; i<N; i++)
    {
        _x[i] = x[i] - alpha * gr[i];
    }
    double result = f(_x, N);

    free(gr);
    free(_x);

    return result;
}