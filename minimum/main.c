#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "methods.h"

typedef double (*RfFunction)(double t, double *x, int n, double *u, int r);

//2.7182818284590452353602874713527

double uk = 0.0;
double *x0 = NULL;
double *xT = NULL;
double *x1 = NULL;
double *u  = NULL;
double *p0 = NULL;
double *p  = NULL;

double f0(double t, double *x, int n)
{
    return x[0]*x[0] + x[1]*x[1] + uk*uk;
}

double f1(double t, double* x, int n)
{
    return 2 * x[0];
}

double f2(double t, double* x, int n)
{
    return x[0] + x[1] + uk;
}

double pf1(double t, double *p, int n)
{
    return 2*x1[0] - p[0] - 2*p[1];
}

double pf2(double t, double *p, int n)
{
    return 2*x1[1] - 2*p[1];
}

double HamiltonPantragen(RmFunction f0, RmFunction *f, double t, double *x, int n, double *u, int r, double *p)
{
    double a = -f0(t, x, n);
    int i=0;
    for (i=0; i<n; i++)
    {
        a = a + f[i](t, x, n) * p[i];
    }
    return a;
}

int main(int argc, char** argv)
{
    int n = 2;
    int r = 1;

    x0 = (double*) malloc( sizeof(double) * n );
    xT = (double*) malloc( sizeof(double) * n );
    x1 = (double*) malloc( sizeof(double) * n );
    u  = (double*) malloc( sizeof(double) * r );
    p  = (double*) malloc( sizeof(double) * n );
    p0 = (double*) malloc( sizeof(double) * n );

    uk = 0.0;
    double h  = 0.00001;
    double t0 = 0.0;
    double T  = +1.0;
    double t1 = 0.1;

    x0[0] = 1.0;
    x0[1] = 2.0;

    RmFunction *f = (RmFunction*) malloc( sizeof(RmFunction*) * n );
    f[0] = f1;
    f[1] = f2;

    runga_kutta2(f, x0, t0, xT, T,  n, h);
    runga_kutta2(f, x0, t0, x1, t1, n, h);

    RmFunction *p_f = (RmFunction*) malloc( sizeof(RmFunction*) * n );
    p_f[0] = pf1;
    p_f[1] = pf2;

    p0[0]  = 2*xT[0];
    p0[1]  = 2*xT[1];
    runga_kutta2(p_f, p0, T, p, t1, n, h);

    printf("x1(1.0) = %.10f x2(1.0) = %.10f\n", xT[0], xT[1]);
    printf("x1(0.1) = %.10f x2(0.1) = %.10f\n", x1[0], x1[1]);
    printf("p1(0.1) = %.10f p2(0.1) = %.10f\n", p[0], p[1]);

    double a = HamiltonPantragen(f0, p_f, t1, x1, n, u, r, p);
    printf("%f\n", a);

    return 0;
}
