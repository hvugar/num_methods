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

//Метод сопряженных градиентов
void conjugate_gradient_method(RnFunction f, R1Function g, double *x, int n, double dx, double epsilon)
{
    int i = 0;
	int j = 0;
    int k = 1;

    double module_grad1 = 0.0;
    double module_grad2 = 0.0;
    double module_grad  = 0.0;
	
	double *s = (double*) malloc(sizeof(double)*n);
    double w = 0.0;

    int ends = 0;

    do {
        j = 0;

        double* grads1 = (double*) malloc(sizeof(double) * n);
        gradient(f, x, n, dx, grads1);
        module_grad1 = grad_module(grads1, n);

        for (i=0; i<n; i++)
            s[i] = -grads1[i];

        do {
            double a,b;
            double alpha0 = 0.0;

            straight_line_search_metod(g, alpha0, 0.001, &a, &b);
            double alpha = golden_section_search_min(g, a, b, epsilon);

            //printf("%4d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", j, x[0], x[1], f(x, n), grads1[0], grads1[1], module_grad1, alpha);

            for (i=0; i<n; i++)
            {
                x[i] = x[i] + alpha * s[i];
            }
			
            double* grads2 = (double*) malloc(sizeof(double) * n);
            gradient(f, x, n, dx, grads2);
            module_grad2 = grad_module(grads2, n);
			
			//printf("%4d %4d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", k, j, x[0], x[1], f(x, n), grads2[0], grads2[1], module_grad2, alpha);

            w = (module_grad2*module_grad2) / (module_grad1*module_grad1);


            for (i=0; i<n; i++)
                s[i] = -grads2[i] + w * s[i];

            module_grad = grad_module(s, n);
			
			printf("%4d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", j, x[0], x[1], f(x, n), grads1[0], grads1[1], module_grad, alpha);
			
			free(grads2);
			
            if (module_grad < epsilon)
            {
                ends = 1;
                break;

            } else {
                if ( (j+1) < n )
                    j = j + 1;
                else
                    break;
            }
        } while (1);
		
		free(grads1);

        if (ends == 1)
            break;
			
		k++;

    } while (1);
	
	free(s);

}
