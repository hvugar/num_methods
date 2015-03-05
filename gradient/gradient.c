#include "methods.h"

void gradient(RnFunction f, double *x, int n, double dx, double *gradients)
{
	int i = 0;
	for (i=0; i<n; i++)
	{
		x[i] = x[i] - dx;

		double f1 = f(x, n);

		x[i] = x[i] + 2*dx;
		
		double f2 = f(x, n);

		x[i] = x[i] - dx;

		gradients[i] = (f2 - f1) / (2 * dx);
	}
}

void gradient1(RnFunction f, double *x, int n, double dx, double *gradients)
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

double grad_module(double *grads, int n)
{
    double result = 0.0;

	int i;
    for (i=0; i<n; i++)
    {
        result += grads[i] * grads[i];
    }

    return sqrt(result);
}

double vertor_norm(double *x, int n)
{
    double result = 0.0;

	int i;
    for (i=0; i<n; i++)
    {
        result = result + x[i] * x[i];
    }

    return sqrt(result);
}

double distance(double *x1, double *x2, int n)
{
	double dist = 0.0;

	int i;
    for (i=0; i<n; i++)
    {
        dist = dist + (x1[i] - x2[i]) * (x1[i] - x2[i]);
    }

    return sqrt(dist);
}