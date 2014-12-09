#include <stdio.h>f
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "minimum.h"
#include "gradient.h"
#include "methods.h"

double f(double *x, int n);
double g(double alpha);

double f_rosenbrock(double *x, int n);
double g_rosenbrock(double alpha);

double  *X;
int     N;
double  epsilon;
double  delta_x;
double	h;

double f1(double x)
{
	return (x-0.133)*(x-0.133)+0.2;
}

int main(int argc, char** argv)
{
    epsilon = 0.001;
    delta_x = 0.00001;
	h = 0.1;
	
	double a1,b1;
	a1 = b1 = 0.0;
	double a2,b2;
	a2 = b2 = 0.0;
	//straight_line_search_metod(f1, +0.10, 0.01, &a1, &b1);
	//straight_line_search_metod(f1, -1.0, 0.01, &a2, &b2);
	
	//printf("%8.6f, %8.6f\n", a1, b1);
	//printf("%8.6f, %8.6f\n", a2, b2);
	
	//return;
    
	N = 2;
    X  = (double*) malloc( sizeof(double) * N );

    X[0]    = -1.2;
    X[1]    = +1.0;
	
//	fast_proximal_gradient_method(f, g, X, N, delta_x, epsilon);
	conjugate_gradient_method(f_rosenbrock, g_rosenbrock, X, N, delta_x, epsilon);
	
	free(X);
	
	return 0;
}

double f_rosenbrock(double *x, int n)
{
	double x1 = x[0];
	double x2 = x[1];
    return ((1-x1)*(1-x1)) + 100*(x2-x1*x1)*(x2-x1*x1);
}

double g_rosenbrock(double alpha)
{
    double *xc = (double*) malloc( sizeof(double) * N );

    double* gr = (double*) malloc(sizeof(double) * N);
    gradient(f_rosenbrock, X, N, delta_x, gr);
	
	int i;
    for (i=0; i<N; i++)
    {
        xc[i] = X[i] - alpha * gr[i];
    }
	
    double result = f_rosenbrock(xc, N);

    free(gr);
    free(xc);

    return result;
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

