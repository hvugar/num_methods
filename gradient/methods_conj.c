#include "methods.h"

double minimize(RnFunction f, double *x, double *grad, int n, double alpha0, double epsilon);

void conjugate_gradient_method(RnFunction f, double *x, int n, double dx, double epsilon)
{
    int step = 0;
    int i = 0;
    double *s = (double*) malloc(sizeof(double) * n);

	printf("%4d %10.8f %10.8f %10.8f\n", step, x[0], x[1], f(x,n));
    do
    {
        // First iteration
        double* gr1 = (double*) malloc(sizeof(double) * n);
        gradient(f, x, n, dx, gr1);
        for (i=0; i<n; i++) s[i] = -gr1[i];

        double gr1_mod = 0.0;
        for (i=0; i<n; i++) gr1_mod += gr1[i]*gr1[i];

        int k = 0;
        do
        {
            double alpha0 = 0.0;
			double alpha = minimize(f, x, s, n, alpha0, epsilon);

            for (i=0; i<n; i++) x[i] = x[i] + alpha * s[i];

			double* gr2 = (double*) malloc(sizeof(double) * n);
            gradient(f, x, n, dx, gr2);

            double gr2_mod = 0.0;
            for (i=0; i<n; i++) gr2_mod += gr2[i]*gr2[i];

            double w = gr2_mod / gr1_mod;
            gr1_mod = gr2_mod;

			for (i=0; i<n; i++) s[i] = -gr2[i] + s[i] * w;

            free(gr2);
            gr2 = NULL;
			step++;
            k++;
			printf("%4d %10.8f %10.8f %10.8f\n", step, x[0], x[1], f(x,n));
        } while ( k < n );
        free(gr1);
        gr1 = NULL;
        double mod_s = s[0]*s[0] + s[1]*s[1];
        if ( mod_s < epsilon ) break;
    } while ( 1 );
    free(s);
    s = NULL;
}

double minimize(RnFunction f, double *x, double *grad, int n, double alpha0, double epsilon)
{
	double argmin(double alpha)
	{
		int j;
		for (j=0; j<n; j++) x[j] = x[j] + alpha * grad[j];
		double result = f(x, n);
		for (j=0; j<n; j++) x[j] = x[j] - alpha * grad[j];
		return result;
	}
	
	double a,b;
	straight_line_search_metod(argmin, alpha0, 0.01, &a, &b);
	double min = golden_section_search_min(argmin, a, b, epsilon);
	
	return min; 
}
