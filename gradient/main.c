#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "minimum.h"
#include "gradient.h"
#include "methods.h"

int count = 0;

double f_rosenbrock(double *x, int n)
{
	double x1 = x[0];
	double x2 = x[1];
	count++;
    return ((1-x1)*(1-x1)) + 100*(x2-x1*x1)*(x2-x1*x1);
}

int main(int argc, char** argv)
{
    double epsilon	= 0.005;		//dovrun sona catma meyari
	double grad_eps	= 0.005;		//gradient
	double line_eps	= 1.0;			//parcani bolme
	double gold_eps	= 0.0001;		//qizil qayda ucun
    
	int n = 2;
    double* x  = (double*) malloc( sizeof(double) * n );
    x[0]    = -1.2;
    x[1]    = +1.0;
	conjugate_gradient_method(f_rosenbrock, x, n, line_eps, gold_eps, grad_eps, epsilon);
	free(x);
	
	printf("%d\n", count);
	return 0;
}

