#include "methods.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "minimum.h"
#include "gradient.h"

//Метод наискорейшего спуска
void fast_proximal_gradient_method(RnFunction f, R1Function g, double* x, int n, double dx, double epsilon)
{
    int i = 0;
    double module_grad = 0;

    do
    {
        i++;

        // minimum yerleshen [a, b]
        double a,b;
        double alpha0 = -100.0;

		puts("1");
//        straight_line_search_metod(g, alpha0, 0.01, &a, &b);
		a = -50;
		b = +50;
		puts("2");

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
		printf("%4d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", i, x[0], x[1], f(x, n), grads[0], grads[1], grad_module(grads, n), alpha);

        free(grads);

    } while ( module_grad > epsilon );
}

//   1  -0.5000  -1.0000   7.3750  -2.2501  -7.9998   8.3102   0.2805
//   2   0.1311   1.2436  -2.2723  -2.9484   0.9744   3.1053   0.2919
//   3   0.9917   0.9591  -3.9965  -0.0491  -0.1633   0.1706   0.2405
//   4   1.0035   0.9984  -4.0000   0.0216  -0.0062   0.0225   0.1705
//   5   0.9999   0.9995  -4.0000  -0.0005  -0.0020   0.0021   0.2737
//   6   1.0000   1.0000  -4.0000   0.0003   0.0002   0.0004   0.0134

//Метод сопряженных градиентов
void conjugate_gradient_method(RnFunction f, R1Function g, double* x, int n, double dx, double epsilon)
{
    int i=0;
	int j = 0;
    int k = 1;

    double module_grad1 = 0;
    double module_grad2 = 0;
    double module_grad  = 0.0;

	
//    double **S = (double**) malloc(sizeof(double*)*n);
//    for (i=0; i<n; i++) S[i] = (double*) malloc(sizeof(double)*n);
	
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

            straight_line_search_metod(g, alpha0, 0.01, &a, &b);
            double alpha = golden_section_search_min(g, a, b, epsilon);

            //printf("%4d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", j, x[0], x[1], f(x, n), grads1[0], grads1[1], module_grad1, alpha);

            for (i=0; i<n; i++)
            {
                x[i] = x[i] + alpha * s[i];
            }
			
            double* grads2 = (double*) malloc(sizeof(double) * n);
            gradient(f, x, n, dx, grads2);
            module_grad2 = grad_module(grads2, n);
			
			printf("%4d %4d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", k, j, x[0], x[1], f(x, n), grads2[0], grads2[1], module_grad2, alpha);

            w = (module_grad2*module_grad2) / (module_grad1*module_grad1);


            for (i=0; i<n; i++)
                s[i] = -grads2[i] + w * s[i];

            module_grad = grad_module(s, n);
			
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

        if (ends == 1)
            break;
			
		k++;

    } while (1);

}
