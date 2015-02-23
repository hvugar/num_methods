#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "minimum.h"
#include "gradient.h"
#include "methods.h"

int count = 0;

void gradient_test();
void svenn_test();

void lagrange(RnFunction f, double *x, int n, 
			  RnFunction* h, int k, double* u, 
			  RnFunction* g, int p, double* v, 
			  double* w);

int main(int argc, char** argv)
{
	svenn_test();
	return 0;
}

void lagrange(RnFunction f, double *x, int n, RnFunction* h, int k, double* u, RnFunction* g, int p, double* v, double* w)
{
	double h = 0.0001;
	double* grads = (double*) malloc( sizeof(double)*n );
	
	f(x,n,h,grads)+w[i]*
}

double f_rosenbrock(double *x, int n)
{
	double x1 = x[0];
	double x2 = x[1];
	count++;
    return ((1-x1)*(1-x1)) + 100*(x2-x1*x1)*(x2-x1*x1);
}

double f_1(double x)
{
	return (10.0 - x)*(10.0 - x);
}

double f_2(double x)
{
	return 2*x*x+16/x;
}

void svenn_test()
{
	double dx = 0.5;
	double x0 = 3.0;
	double a = 6.0;
	double b = 15.0;

//	search_interval_svenn(f_1, x0, dx, &a, &b);
//	printf("a=%8.2f\nb=%8.2f\n", a, b);
//	halph_interval_method(f_1, 0.001, &a, &b);
//	printf("a=%8.2f\nb=%8.2f\n", a, b);
//	search_method_pauella(f_2, 1.0, 1.0, 0.03, &a, &b);

	newton_raphson(f_2, 1.3333333, 0.0001);
}

void gradient_test()
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
}

