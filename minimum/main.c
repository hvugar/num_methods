#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "methods.h"
#include "print.h"
#include "optimal.h"

double _f1(double t, double *x, int n)
{
	return 3.0*x[1]*x[1];// - t;
}

double _f2(double t, double *x, int n)
{
	double u = 0.45;
	return x[0] + x[1] - 2*u - t*t*t + 1.0;
}
/*
double smp1_Dpsi1(double t, double *psi, int n, double u)
{
	double u = 1.0;
	return 2.0 * (x[0] - t*t*t) - psi[1];
}

double smp1_Dpsi2(double t, double *x, double *psi, int n, double u)
{
	return 2.0 * x[1] - 6.0*x[1]*psi[0] - psi[1];
}
*/
int main(int argc, char** argv)
{
	
	double t0 = 0.9;
	double t1 = 1.0;
	double x0[2] = {0.8560880482, 1.0005610153};
	double x1[2] = {0.0, 0.0};
	RmFunction f[] = { _f1, _f2 };
	RungaKuttaSystem(f, t0, x0, t1, x1, 2, 0.000001);
	printf("%.10f %.10f %.10f %.10f\n", x0[0], x0[1], x1[0], x1[1]);

    smp1_control();
    return 0;
}
