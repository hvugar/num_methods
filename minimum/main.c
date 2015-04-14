#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "methods.h"
#include "print.h"
#include "optimal.h"

double smp1_JSum_1(double *t, double *x1, double *x2, double n, double *u, int N)
{
	double sum = 0.0;
	int i=0;
	for (i=0; i<(N-1); i++)
	{
		int j=i+1;
		double fj = (x1[j]-t[j]*t[j]*t[j])*(x1[j]-t[j]*t[j]*t[j]) + x2[j]*x2[j] - t[j]*t[j] + (2*u[j]*u[j] - 1.0)*(2*u[j]*u[j] - 1.0);
		double fi = (x1[i]-t[i]*t[i]*t[i])*(x1[i]-t[i]*t[i]*t[i]) + x2[i]*x2[i] - t[i]*t[i] + (2*u[i]*u[i] - 1.0)*(2*u[i]*u[i] - 1.0);
		sum = sum + 0.5 * (fj+fi) * (t[j]-t[i]);
	}
	sum = sum + (x2[N-1]-1.0) * (x2[N-1]-1.0);
	return sum;
}

double argmin_1(double alpha)
{
	double t[] =    {  0.00000000,   0.10000000,   0.20000000,   0.30000000,   0.40000000,   0.50000000,   0.60000000,   0.70000000,   0.80000000,   0.90000000,   1.00000000};
double u[] =    {  0.82775325,   0.77854243,   0.72393913,   0.66953625,   0.61898677,   0.57434462,   0.53644838,   0.50530151,   0.48044630,   0.46134179,   0.44780773};
double x1[] =   {  0.00000000,   0.00046356,   0.00365694,   0.01197212,   0.02739002,   0.05170808,   0.08695876,   0.13597712,   0.20312210,   0.29520904,   0.42273154};
double x2[] =   {  0.00000000,  -0.06895336,  -0.13500339,  -0.19723767,  -0.25622457,  -0.31379644,  -0.37284390,  -0.43712060,  -0.51103623,  -0.59942282,  -0.70724304};
double gr[] =   {  6.81337833,   6.51907135,   6.50400955,   6.65395836,   6.87821258,   7.10833282,   7.29506020,   7.40392751,   7.41021778,   7.29423137,   7.03774121};
    int N = 11;
    int n = 2;

    int i;
    double u1[N];
    for (i=0; i<N; i++) u1[i] = u[i] - alpha * gr[i];
    return smp1_JSum_1(t, x1, x2, n, u1, N);
}

int main(int argc, char** argv)
{
/*
    double a,b;
    double alpha0 = 0.0;
    straight_line_search_metod(argmin_1, alpha0, 0.01, &a, &b);
    double alpha = golden_section_search_min(argmin_1, a, b, 0.00000001);
    printf("alpha: %.10f %.10f\n", alpha, argmin_1(alpha));
*/
    smp1_control();
    return 0;
}
