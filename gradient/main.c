#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "minimum.h"
#include "gradient.h"
#include "methods.h"

double f_rosenbrock(double *x, int n);
double f(double *x, int n);
double f1(double x);

double  *x;
int     n;
double  epsilon;
double  delta_x;
double	h;

int count = 0;
                                                                                                  
int main(int argc, char** argv)
{
    epsilon = 0.0001;
    delta_x = 0.00001;
    
	n = 2;
    x  = (double*) malloc( sizeof(double) * n );

    x[0]    = -1.2;
    x[1]    = +1.0;
	
//	fast_proximal_gradient_method(f, g, X, N, delta_x, epsilon);
	conjugate_gradient_method(f_rosenbrock, x, n, delta_x, epsilon);
	
	free(x);
	
	return 0;
}

double f_rosenbrock(double *x, int n)
{
	count++;
	double x1 = x[0];
	double x2 = x[1];
    return ((1-x1)*(1-x1)) + 100*(x2-x1*x1)*(x2-x1*x1);
}

double f(double *x, int n)
{
    return pow(x[0],3) + 2*pow(x[1],2) - 3*x[0] - 4*x[1];
}

double f1(double x)
{
	return (x-0.133)*(x-0.133)+0.2;
}