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

int main(int argc, char** argv)
{
    double epsilon	= 0.005;		//dovrun sona catma meyari
	double grad_eps	= 0.005;		//gradient
	double line_eps	= 0.5;			//parcani bolme
	double gold_eps	= 0.0001;		//qizil qayda ucun
    
	int n = 2;
    double* x  = (double*) malloc( sizeof(double) * n );
    x[0]    = -1.2;
    x[1]    = +1.0;
	conjugate_gradient_method(f_rosenbrock, x, n, line_eps, gold_eps, grad_eps, epsilon);
	free(x);
	return 0;
}

double f_rosenbrock(double *x, int n)
{
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