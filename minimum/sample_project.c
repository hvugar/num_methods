#include "methods.h"

void sample_project(RnFunction f, double *x, int n, double line_eps, double gold_eps, double grad_eps, double epsilon, Printer printer)
{
	int count = 0;
	int i;
    double grad_norm = 0.0;
	double a = +2.0;
	double b = +4.0;

    double* grads = (double*) malloc(sizeof(double) * n);
    // Saves last point coordinates
    double *x1 = (double*) malloc(sizeof(double) * n);
    // Used for one dimention minimization for stopring next point coordinates
    double *x2 = (double*) malloc(sizeof(double) * n);
    do
    {   
        gradient(f, x, n, grad_eps, grads);

        double argmin(double alpha)
        {
            int i;
			for (i=0; i<n; i++)
			{
				if ( (x[i] - alpha * grads[i]) < a ) x2[i] = a; else
				if ( (x[i] - alpha * grads[i]) > b ) x2[i] = b; else
				x2[i] = x[i] - alpha * grads[i];
			}
            return f(x2, n);
        }
		
		double alpha = R1Minimize(argmin, 0.001, 0.000001);

		grad_norm = vertor_norm(grads, n);
		
        if (printer != NULL) printer(f, x, n, count, grads, grad_norm, alpha);

        memcpy(x1, x, sizeof(double) * n);
		
        for (i=0; i<n; i++)
        {
            x[i] = x[i] - alpha * grads[i];
			if (x[i] < a) { x[i] = a; }
			if (x[i] > b) { x[i] = b; }
        }
		
		count++;
    }
    while ( grad_norm > epsilon && distance(x1, x, n) > epsilon /*&& fabs(f(x1,n) - f(x,n)) > epsilon*/ );
	if (printer != NULL) printer(f, x, n, count, grads, grad_norm, 0);
    free(x1);
    free(x2);
    free(grads);
}