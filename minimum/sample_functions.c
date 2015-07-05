#include <stdio.h>
#include <math.h>

double f_rosenbrock(double *x, int n)
{
	double x1 = x[0];
	double x2 = x[1];
    return ((1-x1)*(1-x1)) + 100*(x2-x1*x1)*(x2-x1*x1);
}

double f_1(double *x, int n)
{
	return (x[0]-4.0)*(x[0]-4.0)+(x[1]-4.0)*(x[1]-4.0);
}

double h1(double *x, int n)
{
	return x[0]+x[1]-5.0;
}