#include "penalty.h"

void penalty_method(RnFunction f, double *x, int n, RnFunction* h, int m, RnFunction* g, int p, double r)
{
	double epsilon_p = 0.00001;
	
	double G(RnFunction g, double *x, int n)
	{
		return 1.0 / g(x,n);
	}
	
	double H(RnFunction h, double *x, int n)
	{
		return h(x,n) * h(x,n);
	}
	
	double R(double *x, int n)
	{
		int i,j;
		double a = 0.0;
		
		for (i=0; i<m; i++)
		{
			a = a + H(h[i], x, n) / r;
		}
		
		for (j=0; j<p; j++)
		{
			a = a + r * G(g[j], x, n);
		}
		
		return a;
	}
	
	double P(double *x, int n)
	{	
		return f(x,n) + R(x,n);
	}
	
	while ( r * R(x,n) > epsilon_p ) 
	{
	
		// Qoshma qradient usulu ucun parametrler
		double epsilon     = 0.001;       //dovrun sona catma meyari
		double grad_step   = 0.005;       //gradient
		double line_step   = 0.1;         //parcani bolme
		double gold_step   = 0.0001;      //qizil qayda ucun
		conjugate_gradient_method1(P, x, n, line_step, gold_step, grad_step, epsilon);
		//puts("****************************************************************************************************************************");
		
		r = r * 0.1;
	}
}
