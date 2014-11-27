#include "methods.h"
#include "minimum.h"
#include "gradient.h"

//Метод наискорейшего спуска
void fast_proximal_gradient_method(RnFunction f, R1Function g, double* x, int N)
{
    int i = 0;
    double module_grad = 0;
	
	do
    {
        i++;

		// minimum yerleshen [a, b]
        double a,b;
        double alpha0 = 0.0;
        straight_line_search_metod(g, alpha0, 0.01, &a, &b);

		// tapilmish [a, b] parcasinda minimum alpha axtaririq
		// Funksiyanin minimumuniu tapmaq ucun qizil bolgu qaydasinda istifade edib alphani tapiriq
		double alpha = golden_section_search_min(g, a, b, epsilon);

		double* grads = (double*) malloc( sizeof(double) * N );
		gradient(f, x, N, h, grads);

        module_grad = grad_module(grads, N);

        printf("%4d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", i, x[0], x[1], f(x, N), grads[0], grads[1], module_grad, alpha);

		int i;
        for (i=0; i<N; i++)
        {
            x[i] = x[i] - alpha * grads[i];
        }

        free(grads);

    } while ( module_grad > epsilon );	
}