#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "diff_method.h"

double f1(double x, double y)
{
	return x*y/2;
}


double f(double, double);

//Initial condition
double fi(double x);
//Boundary conditions
double m1(double t);
double m2(double t);

void tomas_algorithm(double **a, double* b, double* x, int size);

int main(int argc, char** argv)
{
    double a  = 1.00;
    double dx = 0.10;
    double dt = 0.10;
	
	//implicit_difference_scheme(f, fi, m1, m2, a, dx, dt, 1.0, 1.0);
	puts("");
	//explicit_difference_scheme(f, fi, m1, m2, a, dx, dt, 1.0, 1.0);
	diff_euler(f1, 0.0, 1.0, 0.000001, 2.0);
	diff_runga_kutta(f1, 0.0, 1.0, 0.000001, 2.0);
	printf("\t\t\t%.20f\n", 2.7182818284590452353602874713527);

    return 0;
}

double f(double x, double t) { return -1.0; }

double fi(double x) { return x*x; }

double m1(double t) { return t; }

double m2(double t) { return 1.0+t; }
