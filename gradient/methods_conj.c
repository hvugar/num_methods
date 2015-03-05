#include "methods.h"
#include "print.h"

double minimize(RnFunction f, double *x, double *grad, int n, double alpha0, double line_step, double gold_step);

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
 void conjugate_gradient_method(RnFunction f, double *x, int n, double line_step, double gold_step, double grad_step, double epsilon, Printer printer)
{
    int i = 0;
    int k = 0;
    int iter = 0;
    double mod_s = 0.0;
	double dist = 0.0;
    double *s  = (double*) malloc(sizeof(double) * n);
    double *s1 = (double*) malloc(sizeof(double) * n);
    double *x1 = (double*) malloc(sizeof(double) * n);
    int count = 0;
    
    for (i=0; i<n; i++)
    {
        s[i]  = 0.0;
        s1[i] = 0.0;
    }

    double* gr1 = (double*) malloc(sizeof(double) * n);
    double* gr2 = (double*) malloc(sizeof(double) * n);

    double ss = 0.0;
    double gr1_mod = 0.0;
    double gr2_mod = 0.0;
    do
    {
        // First iteration
        if (k == 0)
        {
            // First direction is gradient direction
            gradient(f, x, n, grad_step, gr1);

            for (i=0; i<n; i++) s[i] = -gr1[i];

            // Module of direction
            ss = grad_module(s, n);
			
            // Divide direction to its module
            for (i=0; i<n; i++) s1[i] = s[i] / ss;

            // Module of gradient
            gr1_mod = 0.0;
            for (i=0; i<n; i++) gr1_mod = gr1_mod + gr1[i]*gr1[i];
        }
        else
        {
            // Calculating gradient in next coordinates
            gradient(f, x, n, grad_step, gr2);

            // Module of next gradient
            gr2_mod = 0.0;
            for (i=0; i<n; i++) gr2_mod = gr2_mod + gr2[i]*gr2[i];

            double w = gr2_mod / gr1_mod;
            gr1_mod = gr2_mod;

            // Calculating direction for next (k+1) iteration
            for (i=0; i<n; i++) s[i] = -gr2[i] + s[i] * w;

            // Module of direction
            ss = grad_module(s, n);
			
            // Divide direction to its module
            for (i=0; i<n; i++) s1[i] = s[i] / ss;
        }

		if (printer != NULL) printer(f, x, n, iter, count, s, s1, gr1, gr2);
        iter++;

        // Minimization in one dimensional direction
        double alpha0 = 0.0;
        double alpha = minimize(f, x, s1, n, alpha0, line_step, gold_step);
        line_step /= 1.2;

		memcpy( x1, x, sizeof(double) * n);
		
        // Calculating next coordinates
        for (i=0; i<n; i++)
        {
            x[i] = x[i] + alpha * s1[i];
        }

        mod_s = 0.0;
        //for (i=0; i<n; i++) mod_s = mod_s + s[i]*s[i];
		mod_s = vertor_norm(s, n);
		dist = distance(x1, x, n);
		
        if ( k == n ) { k = 0; } else { k++; }
		
    } while ( mod_s > epsilon && dist > epsilon );
	
	if (printer != NULL) printer(f, x, n, iter, count, s, s1, gr1, gr2);

    free(gr1);
    free(gr2);
    free(s1);
    free(s);
	free(x1);

    gr1 = gr2 = s1 = s = NULL;
}
 
void conjugate_gradient_method1(RnFunction f, double *x, int n, double line_step, double gold_step, double grad_step, double epsilon, Printer printer)
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
        gradient(f, x, n, grad_step, gr1);

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
			if (printer != NULL) printer(f, x, n, iter, count, s, s1);

            double alpha0 = 0.0;
            double alpha = minimize(f, x, s1, n, alpha0, line_step, gold_step);
            line_step /= 1.2;

            for (i=0; i<n; i++)
            {
                x[i] = x[i] + alpha * s1[i];
            }

            double* gr2 = (double*) malloc(sizeof(double) * n);
            gradient(f, x, n, grad_step, gr2);
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

double minimize(RnFunction f, double *x, double *grad, int n, double alpha0, double line_step, double gold_step)
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
    straight_line_search_metod(argmin, alpha0, line_step, &a, &b);
    double min = golden_section_search_min(argmin, a, b, gold_step);
    return min;
}
