#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "methods.h"

extern void penalty_sample1();
extern void penalty_sample2();
extern void penalty_sample4();
extern void sample_gradient1();

double u = 1.1;

double func(double x, double y)
{
	return  (x/y);
}

double x1(double t, double* x, int n)
{
	return 2*x[0];
}

double x2(double t, double *x, int n)
{
	return x[0]+x[1]+1.0;
}

void sample_runga_kutta2()
{
	int n = 2;
	RmFunction *x = (RmFunction*)malloc(sizeof(RmFunction*)*n);
	x[0] = x1;
	x[1] = x2;
	double *x0 = (double*)malloc(sizeof(double*)*n);
	x0[0] = 1;
	x0[1] = 2;
	double t0 = 0.0;
	double t = 1.0;

	runga_kutta2(x, x0, n, t0, t, 0.00000001);
	printf("%f %f\n", x0[0], x0[1]);
}

int main(int argc, char** argv)
{
//	sample_gradient1();
//	penalty_sample4();

	

	return 0;
}