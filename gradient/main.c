#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "minimum.h"
#include "gradient.h"
#include "methods.h"

int count = 0;
double r = 0.0;

void gradient_test();
void straf_test();
void svenn_test();

void lagrange(RnFunction f, double *x, int n, 
			  RnFunction* h, int k, double* u, 
			  RnFunction* g, int p, double* v, 
			  double* w);

			  
double f(double *x, int n)
{
	return (x[0]-4.0)*(x[0]-4.0)+(x[1]-4.0)*(x[1]-4.0);
}

double h(double *x, int n)
{
	return x[0]+x[1]-5.0;
}

double P(double *x, int n)
{
	return f(x,n) + (1.0/r) * h(x,n) * h(x,n);
}

void test1(RnFunction *fs, double* x, int n)
{
	int i=0;
	for ( i = 0; i < n; i++)
	{
		printf("%f\n", fs[i](x, n));
	}
}


int main(int argc, char** argv)
{
	//gradient_test();
	r = 10000000000.00;
	straf_test();
	
/*	double* x = (double*)malloc( sizeof(double) * 2 );
    x[0]    = +5.0;
    x[1]    = +5.0;
	
	RnFunction *fs = (RnFunction*) malloc ( sizeof(RnFunction*) * 3);
	fs[0] = f;
	fs[1] = h;
	fs[2] = P;
	
	test1(fs, x, 3);
	*/
	return 0;
}

void straf_test()
{
    double epsilon	= 0.001;		//dovrun sona catma meyari
	double grad_eps	= 0.005;		//gradient
	double line_eps	= 0.1;			//parcani bolme
	double gold_eps	= 0.0001;		//qizil qayda ucun
	
	int n = 2;
	double* x = (double*)malloc( sizeof(double) * n );
	
    x[0]    = +5.0;
    x[1]    = +5.0;
	while ( r * (h(x,n)) > epsilon ) {
		conjugate_gradient_method(P, x, n, line_eps, gold_eps, grad_eps, epsilon);
		r = r * 0.1;
	}
	
	free(x);
}

void lagrange(RnFunction f, double *x, int n, RnFunction* h, int k, double* u, RnFunction* g, int p, double* v, double* w)
{
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
//	double dx = 0.5;
//	double x0 = 3.0;
//	double a = 6.0;
//	double b = 15.0;

//	search_interval_svenn(f_1, x0, dx, &a, &b);
//	printf("a=%8.2f\nb=%8.2f\n", a, b);
//	halph_interval_method(f_1, 0.001, &a, &b);
//	printf("a=%8.2f\nb=%8.2f\n", a, b);
//	search_method_pauella(f_2, 1.0, 1.0, 0.03, &a, &b);

	newton_raphson(f_2, 1.3333333, 0.0001);
}

void gradient_test()
{
    double epsilon	= 0.001;		//dovrun sona catma meyari
	double grad_eps	= 0.005;		//gradient
	double line_eps	= 0.1;			//parcani bolme
	double gold_eps	= 0.0001;		//qizil qayda ucun
    
	int n = 2;
    double* x  = (double*) malloc( sizeof(double) * n );
    x[0]    = -1.2;
    x[1]    = +1.0;
	conjugate_gradient_method(f_rosenbrock, x, n, line_eps, gold_eps, grad_eps, epsilon);
	free(x);
	
	printf("%d\n", count);
}
