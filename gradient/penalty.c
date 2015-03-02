#include "penalty.h"

void penalty_method(RnFunction f, double *x, int n, RnFunction* h, int m, RnFunction* g, int p, double r)
{
	double epsilon_p = 0.00001;
	
	double P(double *x, int n)
	{
		int i,j;
		double a = 0.0;
	
		a = f(x,n);

		for (i=0; i<m; i++)
		{
			a += (1.0 / r) * h[i](x,n) * h[i](x,n);
			printf("%f\n", a);
		}
		for (j=0; j<p; j++)
		{
			a += (1.0 / r) * g[j](x,n) * g[j](x,n);
		}
		
		return a;
	}
	
	double R(double *x, int n)
	{
		int i,j;
		double a = 0.0;
		for (i=0; i<m; i++)
		{
			a += (1.0 / r) * h[i](x,n) * h[i](x,n);
		}
		for (j=0; j<p; j++)
		{
			a += (1.0 / r) * g[i](x,n) * g[i](x,n);
		}
		
		return a;
	}
	
	while ( r * R(x,n) > epsilon_p ) 
	{
	
		// Qoshma qradient usulu ucun parametrler
		double epsilon	= 0.001;		//dovrun sona catma meyari
		double grad_eps	= 0.005;		//gradient
		double line_eps	= 0.1;			//parcani bolme
		double gold_eps	= 0.0001;		//qizil qayda ucun
		conjugate_gradient_method(P, x, n, line_eps, gold_eps, grad_eps, epsilon);
		puts("****************************************************************************************************************************");
		
		r = r * 0.1;
	}
}