#include "penalty.h"

void penalty_method(RnFunction f, double *x, int n, RnFunction* h, int m, RnFunction* g, int p, double r1, double r2, double epsilon)
{
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
		int i;
		double sum = 0.0;
		for (i=0; i<p; i++) sum = sum + r1 * G(g[i], x, n);
		for (i=0; i<m; i++) sum = sum + r2 * H(h[i], x, n);
		return sum;
	}
	
	double P(double *x, int n)
	{	
		return f(x,n) + R(x,n);
	}

	// Qoshma qradient usulu ucun parametrler
	double min_epsilon = 0.001;       //dovrun sona catma meyari
	double grad_step   = 0.005;       //gradient
	double line_step   = 0.1;         //parcani bolme
	double gold_step   = 0.0001;      //qizil qayda ucun
	
	double* x1 = (double*) malloc( sizeof(double) * n );
	do
	{
		memcpy( x1, x, sizeof(double) * n );
		printf("\nr1 = %.10f\nr2 = %.10f\n", r1, r2);
		printf("Minimization...\n");
		conjugate_gradient_method(P, x, n, line_step, gold_step, grad_step, min_epsilon, printer2);
		printf("Minimized...\nx1 = %.10f\nx2 = %.10f\n", x[0], x[1]);
		printf("********************************************************\n");
		r1 = r1 * 0.10;
		r2 = r2 * 10.0;
		double dist = 0.0;
		int i=0;
		for (i=0; i<n; i++) dist = dist + (x[i]-x1[i])*(x[i]-x1[i]);
		dist = sqrt(dist);
		if (dist < epsilon) break;		
	} while (1);
	free(x1);

	/*
	while ( r1 * R(x,n) > epsilon_p ) 
	{
		printf("\nr1 = %.10f\nr2 = %.10f\n", r1, r2);
		printf("Minimization...\n");
		conjugate_gradient_method(P, x, n, line_step, gold_step, grad_step, min_epsilon, printer2);
		printf("Minimized...\nx1 = %.10f\nx2 = %.10f\n", x[0], x[1]);
		printf("********************************************************\n");
		
		r1 = r1 * 0.10;
		r2 = r2 * 10.0;
	}
	*/
}
