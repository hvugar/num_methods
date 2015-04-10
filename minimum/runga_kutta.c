#include "methods.h"

double runga_kutta1(R2Function f, double y0, double x0, double x, double h)
{
    while (x0 <= x)
    {
        double k1 = f(x0, y0);
        double k2 = f(x0+h/2.0, y0+(h/2.0)*k1);
        double k3 = f(x0+h/2.0, y0+(h/2.0)*k2);
        double k4 = f(x0+h, y0+h*k3);

        y0 = y0 + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
        x0 = x0 + h;
    }
    return y0;
}

void runga_kutta2(RmFunction *f, double *y0, double x0, double *y1, double x1, int n, double h)
{
    double *k1 = (double*) malloc( sizeof(double) * n );
    double *k2 = (double*) malloc( sizeof(double) * n );
    double *k3 = (double*) malloc( sizeof(double) * n );
    double *k4 = (double*) malloc( sizeof(double) * n );
    double *y  = (double*) malloc( sizeof(double) * n );

    if (fabs(x0-x1) > h)
    {
        if ( x0 < x1 ) h = +fabs(h);
        if ( x0 > x1 ) h = -fabs(h);

        double x = x0;
        memcpy(y1, y0, sizeof(double)*n);
        while (1)
        {
            int i=0;
            // Calculating k1 vector
            for (i = 0; i<n; i++) k1[i] = f[i](x, y1, n);

            // Calculating k2 vector
            for (i = 0; i<n; i++) y[i]  = y1[i] + (h/2.0) * k1[i];
            for (i = 0; i<n; i++) k2[i] = f[i](x+h/2.0, y, n);

            // Calculating k3 vector
            for (i = 0; i<n; i++) y[i]  = y1[i] + (h/2.0) * k2[i];
            for (i = 0; i<n; i++) k3[i] = f[i](x+h/2.0, y, n);

            // Calculating k1 vector
            for (i = 0; i<n; i++) y[i]  = y1[i] + h * k3[i];
            for (i = 0; i<n; i++) k4[i] = f[i](x+h, y, n);

            // Calculating y
            for (i = 0; i<n; i++) y1[i] = y1[i] + (h/6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);

            x = x + h;

            if ( x0 < x1 && x > x1 ) break;
            if ( x0 > x1 && x < x1 ) break;
        }
    }
	
    free(y);
    free(k1);
    free(k2);
    free(k3);
    free(k4);
}

double integeral_trapezoidal_rule1(double *y, double *x, int n)
{
	int i=0;
	double sum = 0.0;
	
	for (i=0; i<(n-1); i++)
	{
		sum += (( y[i+1] + y[i] ) * ( x[i+1] - x[i] )) / 2.0;
	}
	return sum;
}

double integeral_trapezoidal_rule2(R1Function f, double *x, int n)
{
	int i=0;
	double sum = 0.0;
	
	for (i=0; i<(n-1); i++)
	{
		sum += (( f(x[i+1]) + f(x[i]) ) * ( x[i+1] - x[i] )) / 2.0;
	}
	return sum;
}

