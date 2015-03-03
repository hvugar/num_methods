#include "methods.h"

double minimize(RnFunction f, double *x, double *grad, int n, double alpha0, double step, double epsilon);

/**
 * @brief Метод сопряженных градиентов Флетчера — Ривса
 * @param f Целевая функция
 * @param x Независимые переменные n - измерение
 * @param n
 * @param line_eps
 * @param gold_eps
 * @param grad_eps
 * @param epsilon
 */
void conjugate_gradient_method(RnFunction f, double *x, int n, double line_eps, double gold_eps, double grad_eps, double epsilon)
{
    int iter = 0;
    int i = 0;
    double *s = (double*) malloc(sizeof(double) * n);
    double *s1 = (double*) malloc(sizeof(double) * n);
    int count = 0;
    
    for (i=0; i<n; i++)
    {
        s[i] = 0.0;
        s1[i] = 0.0;
    }

    do
    {
        // First iteration
        double* gr1 = (double*) malloc(sizeof(double) * n);
        gradient(f, x, n, grad_eps, gr1);

        for (i=0; i<n; i++)
            s[i] = -gr1[i];

        double ss = grad_module(s, n);

        for (i=0; i<n; i++)
            s1[i] = s[i] / ss;

        double gr1_mod = 0.0;
        for (i=0; i<n; i++)
            gr1_mod += gr1[i]*gr1[i];

        int k = 0;
        do
        {
            print(iter, x, s, s1, n, f, count);

            double alpha0 = 0.0;
            double alpha = minimize(f, x, s1, n, alpha0, line_eps, gold_eps);
            line_eps /= 1.2;

            for (i=0; i<n; i++)
            {
                x[i] = x[i] + alpha * s1[i];
            }

            double* gr2 = (double*) malloc(sizeof(double) * n);
            gradient(f, x, n, grad_eps, gr2);
            count += 2*n+1;

            double gr2_mod = 0.0;
            for (i=0; i<n; i++)
                gr2_mod += gr2[i]*gr2[i];

            double w = gr2_mod / gr1_mod;
            gr1_mod = gr2_mod;

            for (i=0; i<n; i++)
                s[i] = -gr2[i] + s[i] * w;

            ss = grad_module(s, n);
            for (i=0; i<n; i++)
                s1[i] = s[i] / ss;

            free(gr2);
            gr2 = NULL;

            iter++;
            k++;
        } while ( k < n );
        
        free(gr1);
        gr1 = NULL;
        
        double mod_s = s[0]*s[0] + s[1]*s[1];
        if ( mod_s < epsilon ) break;
    } while ( 1 );
    
    free(s1);
    s1 = NULL;
    free(s);
    s = NULL;
}

void conjugate_gradient_method2(RnFunction f, double *x, int n, double line_eps, double gold_eps, double grad_eps, double epsilon)
{
    int i = 0;
	int k = 0;
    int iter = 0;
	double mod_s = 0.0;
    double *s =  (double*) malloc(sizeof(double) * n);
    double *s1 = (double*) malloc(sizeof(double) * n);
    int count = 0;
    
    for (i=0; i<n; i++)
    {
        s[i]  = 0.0;
        s1[i] = 0.0;
    }
	
	double* gr1 = (double*) malloc(sizeof(double) * n);
	double* gr2 = (double*) malloc(sizeof(double) * n);
	
	do
	{
		double ss = 0.0;
		double gr1_mod = 0.0;
		double gr2_mod = 0.0;
		// First iteration
		if (k == 0)
		{
			// First direction is gradient direction
			gradient(f, x, n, grad_eps, gr1);
			
			for (i=0; i<n; i++) s[i] = -gr1[i];
			
			// Module of direction
			ss = grad_module(s, n);
			
			// Divide direction to its module
			for (i=0; i<n; i++) s1[i] = s[i] / ss;
			
			// Module of gradient
	        for (i=0; i<n; i++) gr1_mod += gr1[i]*gr1[i];
		}
		
		print(iter, x, s, s1, n, f, count);
		iter++;
		
		// Minimization in one dimensional direction
		double alpha0 = 0.0;
		double alpha = minimize(f, x, s1, n, alpha0, line_eps, gold_eps);
		line_eps /= 1.2;
		
		// Calculating next coordinates
	    for (i=0; i<n; i++) 
		{
			x[i] = x[i] + alpha * s1[i];
		}
		
		// Calculating gradient in next coordinates
		gradient(f, x, n, grad_eps, gr2);
		
		// Module of next gradient
        for (i=0; i<n; i++) gr2_mod += gr2[i]*gr2[i];
		
        double w = gr2_mod / gr1_mod;
        gr1_mod = gr2_mod;
		
		// Calculating direction for next (k+1) iteration
        for (i=0; i<n; i++) s[i] = -gr2[i] + s[i] * w;
		
		// Module of direction
		ss = grad_module(s, n);
		
		// Divide direction to its module
		for (i=0; i<n; i++) s1[i] = s[i] / ss;
		
		if ( k == n ) k = 0;
		
		mod_s = s[0]*s[0] + s[1]*s[1];
		
	} while (mod_s > epsilon);
	
	free(gr1);
	free(gr2);
	free(s1);
	free(s);
	
	gr1 = gr2 = s1 = s = NULL;
}

double minimize(RnFunction f, double *x, double *grad, int n, double alpha0, double line_eps, double gold_eps)
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
    straight_line_search_metod(argmin, alpha0, line_eps, &a, &b);
    double min = golden_section_search_min(argmin, a, b, gold_eps);
    return min;
}

void conjugate_gradient_method1(RnFunction f, double *x, int n, double line_eps, double gold_eps, double grad_eps, double epsilon)
{
    double* s = (double*) malloc(sizeof(double) * n);
    double* g1 = (double*) malloc(sizeof(double) * n);
    double* g2 = (double*) malloc(sizeof(double) * n);
    double module_s = 0.0;
    
	printf("%8.4f %8.4f %8.4f\n", x[0], x[1], f(x,n));
	
    do
	{
        int i=0;
        int k = 0;

        gradient(f, x, n, grad_eps, g1);

        if (k == 0)
		{
            for (i=0; i<n; i++)
            {
                s[i] = -g1[i];
            }
		}

        double alpha0 = 0.0;
        double alpha = minimize(f, x, s, n, alpha0, line_eps, gold_eps);
//		double alpha = 0.1;

        for (i=0; i<n; i++)
        {
            x[i] = x[i] + alpha * s[i];
        }

        gradient(f, x, n, grad_eps, g2);

        double w = 0;
        double w1 = 0.0;
        double w2 = 0.0;
        for (i=0; i<n; i++)
        {
            w1 = w1 + g1[i]*g1[i];
            w2 = w2 + g2[i]*g2[i];
        }
        w = w2/w1;

        for (i=0; i<n; i++)
        {
            s[i] = g2[i] + s[i] * w;
        }
        k += 1;

        if ( k == n )
        {
            k = 0;
        }

		printf("%8.4f %8.4f %8.4f\n", x[0], x[1], f(x,n));
		
		module_s = 0.0;
		for (i=0; i<n; i++)
        {
            module_s = module_s + s[i] * s[i];
        }
        module_s = sqrt(module_s);
		
    } while ( module_s > epsilon);
	puts("***");

    free(g1);
    free(g2);
    free(s);
}
