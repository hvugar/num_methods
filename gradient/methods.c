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


void print(double* x, int step, double *gr, RnFunction f)
{
    printf("%6d   x1=%8.6f   x2=%8.6f   f=%8.6f   grad1=%8.6f   grad2=%8.6f\n", step, x[0], x[1], f(x,2), gr[0], gr[1]);
}

void conjugate_gradient_method(RnFunction f, R1Function g, double *x, int n, double dx, double epsilon)
{
    int step = 0;
    int i = 0;

    double *s = (double*) malloc(sizeof(double) * n);

	printf("%4d %8.6f %8.6f\n", step, x[0], x[1]);
    do
    {
        // First iteration
        double* gr1 = (double*) malloc(sizeof(double) * n);
        gradient(f, x, n, dx, gr1);
        for (i=0; i<n; i++) s[i] = -gr1[i];

        double gr1_mod = 0.0;
        for (i=0; i<n; i++) gr1_mod += gr1[i]*gr1[i];

        //print(x, step, gr1, f);

        int k = 0;
        do
        {
            double a,b;
            double alpha0 = 0.0;
            straight_line_search_metod(g, alpha0, 0.01, &a, &b);
            double alpha = golden_section_search_min(g, a, b, epsilon);
            //double alpha = (a+b)/2;

            for (i=0; i<n; i++)
            {
                x[i] = x[i] + alpha * s[i];
            }

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
			printf("%4d %8.6f %8.6f\n", step, x[0], x[1]);

            k++;
        } while ( k < n );

        free(gr1);
        gr1 = NULL;

        double mod_s = s[0]*s[0] + s[1]*s[1];

        if ( mod_s < epsilon ) break;

    } while ( 1 );

    free(s);
    s = NULL;
}

void conjugate_gradient_method2(RnFunction f, R1Function g, double *x, int n, double dx, double epsilon)
{
    int step = 0;
    int i = 0;
    double *s = (double*) malloc(sizeof(double) * n);

    double mod_s = 0.0;
    do
    {
        int k = 0;

        do {
            k++;
            step++;

            printf("%4d %8.4f %8.4f\n", step, x[0], x[1]);

            double* gr1 = (double*) malloc(sizeof(double) * n);
            gradient(f, x, n, dx, gr1);

            for (i=0; i<n; i++) s[i] = -gr1[i];

            double gr_mod1 = 0.0;
            for (i=0; i<n; i++) gr_mod1 += gr1[i]*gr1[i];

            double a,b;
            double alpha0 = 0.0;

            straight_line_search_metod(g, alpha0, 0.001, &a, &b);

            double alpha = golden_section_search_min(g, a, b, epsilon);

            for (i=0; i<n; i++)
            {
                x[i] = x[i] + alpha * s[i];
            }

            double* gr2 = (double*) malloc(sizeof(double) * n);
            gradient(f, x, n, dx, gr2);
            double gr_mod2 = 0.0;
            for (i=0; i<n; i++) gr_mod2 += gr2[i]*gr2[i];

            double w = gr_mod2 / gr_mod1;

            for (i=0; i<n; i++)
                s[i] = -gr2[i] + w * s[i];

            free(gr1);
            gr1 = NULL;

            free(gr2);
            gr2 = NULL;

            //if (step < 100)
            //			printf("%4d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", step, x[0], x[1], f(x, n), mod_s, s[0], s[1]);
            //printf("%4d %8.4f %8.4f\n", step, x[0], x[1]);
        } while( k < n );

        mod_s = s[0]*s[0] + s[1]*s[1];
        //printf("%8.4f < %8.4f\n", mod_s, epsilon);

        if (mod_s < epsilon) break;

    } while (1);

    free(s);
}

// Метод сопряженных градиентов Флетчера — Ривса 
void conjugate_gradient_method1(RnFunction f, R1Function g, double *x, int n, double dx, double epsilon)
{
    int step = 0;
    int i=0;

    do
    {
        double *s = (double*) malloc(sizeof(double) * n);

        double* grads1 = (double*) malloc(sizeof(double) * n);
        gradient(f, x, n, dx, grads1);

        for (i=0; i<n; i++) s[i] = -grads1[i];

        print(x, step, grads1, f);
        return;

        double grad_mod1;
        double grad_mod2;

        grad_mod1 = grads1[0]*grads1[0] + grads1[1]*grads1[1];


        int k=0;
        double mod_s = 0.0;

        do
        {
            step++;

            double a,b;
            double alpha0 = 0.0;

            double func1(double alpha1)
            {
                int j;
                for (j=0; j<n; j++)
                {
                    x[j] = x[j] - alpha1 * grads1[j];
                }

                double result = f(x, n);

                for (j=0; j<n; j++)
                {
                    x[j] = x[j] + alpha1 * grads1[j];
                }

                return result;
            };

            straight_line_search_metod(func1, alpha0, 0.001, &a, &b);

            double alpha = golden_section_search_min(func1, a, b, epsilon);

            for (i=0; i<n; i++)
            {
                x[i] = x[i] + alpha * s[i];
            }

            double* grads2 = (double*) malloc(sizeof(double) * n);
            gradient(f, x, n, dx, grads2);
            grad_mod2 = grads2[0]*grads2[0] + grads2[1]*grads2[1];
            
            double w = grad_mod2 / grad_mod1;

            for (i=0; i<n; i++)
                s[i] = -grads2[i] + w * s[i];

            free(grads2);

            grad_mod1 = grad_mod2;

            k++;
        } while (k<n);

        mod_s = s[0]*s[0] + s[1]*s[1];

        //if (step%20==0 || mod_s < epsilon)
        //printf("%4d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4\n", step, x[0], x[1], f(x, n), mod_s, s[0], s[1]);

        free(s);
        s = NULL;
        free(grads1);
        grads1 = NULL;

        if (mod_s < epsilon) break;

        int c;
        //if (step % 100 == 0) scanf("%d", &c);
        //printf("%4d %10.6f %10.6f\n", step, x[0], x[1]);

    } while (1);
}
