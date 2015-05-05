#include "methods.h"

void gradient(RnFunction *f, double *x, int n, double dx, double *gradients)
{
    gradient2(f, x, n, dx, gradients);
}

void gradient1(RnFunction *f, double *x, int n, double dx, double *gradients)
{
    double f0 = f->fx(x, n);
    int i;
    for (i=0; i<n; i++)
    {
        x[i] = x[i] + dx;
        double f1 = f->fx(x, n);
        gradients[i] = (f1 - f0) / dx;
        x[i] = x[i] - dx;
    }
}

void gradient2(RnFunction *f, double *x, int n, double dx, double *gradients)
{
    int i = 0;
    for (i=0; i<n; i++)
    {
        x[i] = x[i] - dx;
        double f1 = f->fx(x, n);
        x[i] = x[i] + 2*dx;
        double f2 = f->fx(x, n);
        x[i] = x[i] - dx;
        gradients[i] = (f2 - f1) / (2 * dx);
    }
}
