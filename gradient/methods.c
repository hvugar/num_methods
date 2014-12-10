#include "methods.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "minimum.h"
#include "gradient.h"

//Метод наискорейшего спуска
void fast_proximal_gradient_method(RnFunction f, R1Function g, double *x, int n, double dx, double epsilon)
{
    int i = 0;
    double module_grad = 0;

    do
    {
        i++;

        // minimum yerleshen [a, b]
        double a,b;
        double alpha0 = 0.0;

        straight_line_search_metod(g, alpha0, 0.001, &a, &b);

        // tapilmish [a, b] parcasinda minimum alpha axtaririq
        // Funksiyanin minimumuniu tapmaq ucun qizil bolgu qaydasinda istifade edib alphani tapiriq
        double alpha = golden_section_search_min(g, a, b, epsilon);

        double* grads = (double*) malloc(sizeof(double) * n);
        gradient(f, x, n, dx, grads);

        module_grad = grad_module(grads, n);

        printf("%4d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", i, x[0], x[1], f(x, n), grads[0], grads[1], module_grad, alpha);

        int j;
        for (j=0; j<n; j++)
        {
            x[j] = x[j] - alpha * grads[j];
        }

        gradient(f, x, n, dx, grads);
        //printf("%4d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", i, x[0], x[1], f(x, n), grads[0], grads[1], grad_module(grads, n), alpha);

        free(grads);

    } while ( module_grad > epsilon );
}

void conjugate_gradient_method(RnFunction f, R1Function g, double *x, int n, double dx, double epsilon)
{
    int step = 0;
    int i = 0;
    double *s = (double*) malloc(sizeof(double) * n);

	//printf("%4d %8.6f %8.6f\n", step, x[0], x[1]);
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
            double a,b;
            double alpha0 = 0.0;
			
			double min(double alpha)
			{
				int j;
				for (j=0; j<n; j++) x[j] = x[j] - alpha * gr1[j];
				double result = f(x, n);
				for (j=0; j<n; j++) x[j] = x[j] + alpha * gr1[j];
				return result;
			}
			
            straight_line_search_metod(min, alpha0, 0.01, &a, &b);
            double alpha = golden_section_search_min(min, a, b, epsilon);

            for (i=0; i<n; i++)
                x[i] = x[i] + alpha * s[i];

			double* gr2 = (double*) malloc(sizeof(double) * n);
            gradient(f, x, n, dx, gr2);

            double gr2_mod = 0.0;
            for (i=0; i<n; i++) gr2_mod += gr2[i]*gr2[i];

            double w = gr2_mod / gr1_mod;
            gr1_mod = gr2_mod;

			for (i=0; i<n; i++)
				s[i] = -gr2[i] + s[i] * w;

            free(gr2);
            gr2 = NULL;
			step++;
            k++;
			printf("%4d %8.6f %8.6f\n", step, x[0], x[1]);
        } while ( k < n );
        free(gr1);
        gr1 = NULL;
        double mod_s = s[0]*s[0] + s[1]*s[1];
        if ( mod_s < epsilon ) break;
    } while ( 1 );
    free(s);
    s = NULL;
}