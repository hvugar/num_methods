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

void runga_kutta2(RmFunction *f, double *y0, int n, double x0, double x, double h)
{
	double *k1 = (double*) malloc( sizeof(double) * n );
	double *k2 = (double*) malloc( sizeof(double) * n );
	double *k3 = (double*) malloc( sizeof(double) * n );
	double *k4 = (double*) malloc( sizeof(double) * n );
	double *y1 = (double*) malloc( sizeof(double) * n );
	
	int i=0;
	while (x0 <= x)
	{
		for (i = 0; i<n; i++)
		{
			k1[i] = f[i](x0, y0, n);
			y1[i] = y0[i]+(h/2.0)*k1[i];
		}
		
		
		for (i = 0; i<n; i++)
		{
			k2[i] = f[i](x0+h/2.0, y1, n);
			y1[i] = y0[i]+(h/2.0)*k2[i];
		}
		
		for (i = 0; i<n; i++)
		{
			k3[i] = f[i](x0+h/2.0, y1, n);
			y1[i] = y0[i]+h*k3[i];
		}
		
		for (i = 0; i<n; i++)
		{
			k4[i] = f[i](x0+h, y1, n);
			y1[i] = y0[i]+h*k3[i];
		}
		
		for (i = 0; i<n; i++)
		{
			y0[i] = y0[i] + (h/6) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
		}
		
		x0 = x0 + h;
	}
	
	free(y1);
	free(k1);
	free(k2);
	free(k3);
	free(k4);
}
