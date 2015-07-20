#include "methods.h"

void gradient(RnFunction f, double *x, int n, double dx, DoubleVector& gr)
{
    gradient2(f, x, n, dx, gr);
}

void gradient1(RnFunction f, double *x, int n, double dx, DoubleVector& gr)
{
    double f0 = f(x, n);
    int i;
    for (i=0; i<n; i++)
    {
        x[i] = x[i] + dx;
        double f1 = f(x, n);
        gr[i] = (f1 - f0) / dx;
        x[i] = x[i] - dx;
    }
}

void gradient2(RnFunction f, double *x, int n, double dx, DoubleVector& gr)
{
    int i = 0;
    for (i=0; i<n; i++)
    {
        x[i] = x[i] - dx;
        double f1 = f(x, n);
        x[i] = x[i] + 2*dx;
        double f2 = f(x, n);
        x[i] = x[i] - dx;
        gr[i] = (f2 - f1) / (2 * dx);
    }
}

double derivative1(R1Function f, double x, double h)
{
    return (f(x+h)-f(x-h)) / (2*h);
}


double derivative2(R1Function f, double x, double h)
{
    return (f(x+2*h)-2*f(x)+f(x-2*h)) / (4*h*h);
}

double vertor_norm(double *vctr, int n)
{
    double sum = 0.0;
    int i;
    for (i=0; i<n; i++)
    {
        sum += vctr[i] * vctr[i];
    }
    return sqrt(sum);
}

double grad_module(double *gr, int n)
{
    return vertor_norm(gr, n);
}

double distance(double *x1, double *x2, int n)
{
    double dist = 0.0;
    int i;
    for (i=0; i<n; i++)
    {
        double dx = x1[i] - x2[i];
        dist += dx*dx;
    }
    return sqrt(dist);
}
