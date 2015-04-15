#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "methods.h"
#include "print.h"
#include "optimal.h"

double _f1(double t, double *x, int n)
{
	//printf("_f1 %.10f\n", t);
	return 3.0*x[1]*x[1];// - t;
}

double _f2(double t, double *x, int n)
{
	//printf("_f2 %.10f\n", t);
	return x[0] + x[1] + 1.0;
}

int main(int argc, char** argv)
{
	
	double t0 = 0.0;
	double t1 = 0.1;
	double x0[2] = {0.0, 0.0};
	double x1[2] = {0.0, 0.0};
	RmFunction f[] = { _f1, _f2 };
	RungaKuttaSystem(f, t0, x0, t1, x1, 2, 0.000001);
	printf("%.8f %.8f %.8f %.8f\n", x0[0], x0[1], x1[0], x1[1]);

    smp1_control();
    return 0;
}
