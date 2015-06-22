#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "method_grid.h"

//double u(double x, double t) { return 2.0*sin(t) + x*x*x*x + 5.0; }
//double fi(double x) { return x*x*x*x + 5.0; };
//double m1(double t) { return 2.0 * sin(t) + 5.0; }
//double m2(double t) { return 2.0 * sin(t) + 1.0 + 5.0; }
//double f1(double x, double t) { return x*x*x*x + 5.0; }

double u(double x, double t) { return x*x*x*x + t*t*t + 2.0; }
double fi(double x) { return x*x*x*x + 2.0; };
double m1(double t) { return t*t*t + 2.0; }
double m2(double t) { return t*t*t + 3.0; }
double f1(double x, double t) { return 3.0*t*t - 12.0*x*x; }

void calculate_grid()
{
	double alpha = 1.0;
    double dx = 0.01;
    double dt = 0.01;
    double x0 = 0.0;
    double x1 = 1.0;
    double t0 = 0.0;
    double t1 = 1.0;
	
	grid g;
    implicit_difference_scheme(f1, fi, m1, m2, alpha, dx, dt, x0, x1, t0, t1, &g);
	
    int i,j;
    printf("n=%d m=%d\n", g.n, g.m);
    for (j=0; j<g.m; j++)
    {
        for (i=0; i<g.n; i++)
        {
            if (i%(g.n/10)==0)
                printf("%12.8f", g.u[j][i]);
        }
		puts("");
		for (i=0; i<g.n; i++)
        {
            if (i%(g.n/10)==0)
                printf("%12.8f", u(dx*i, dt*j));
        }
		puts("\n---");
    }
    return 0;
}