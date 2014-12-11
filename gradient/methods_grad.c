#include "methods.h"

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