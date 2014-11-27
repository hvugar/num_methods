#include "gradient.h"
#include <math.h>

void gradient(RnFunction f, double *x, int n, double dx, double* gradients)
{
    double f0 = f(x, n);

	int i;
    for (i=0; i<n; i++)
    {
        x[i] = x[i] + dx;
        double f1 = f(x, n);

        gradients[i] = (f1 - f0) / dx;
        x[i] = x[i] - dx;
    }
}

double grad_module(double* grads, int n)
{
    double result = 0.0;

	int i;
    for (i=0; i<n; i++)
    {
        result += grads[i] * grads[i];
    }

    return sqrt(result);
}
