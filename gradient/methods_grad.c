#include "methods.h"

double minimize1(RnFunction f, double *x, double *grad, int n,  double alpha0, double epsilon);

//Метод наискорейшего спуска
void fast_proximal_gradient_method(RnFunction f, double *x, int n, double dx, double epsilon)
{
    int i = 0;
    double module_grad = 0;

    do
    {
        i++;

        double* grads = (double*) malloc(sizeof(double) * n);
        gradient(f, x, n, dx, grads);
        double alpha0 = 0.0;
		
        // Funksiyanin minimumuniu tapmaq ucun qizil bolgu qaydasinda istifade edib alphani tapiriq
        double alpha = minimize1(f, x, grads, n, alpha0, epsilon);

        module_grad = grad_module(grads, n);

        printf("%4d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", i, x[0], x[1], f(x, n), grads[0], grads[1], module_grad, alpha);

        int j;
        for (j=0; j<n; j++)
        {
            x[j] = x[j] - alpha * grads[j];
        }

        //printf("%4d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", i, x[0], x[1], f(x, n), grads[0], grads[1], grad_module(grads, n), alpha);

        free(grads);

    } while ( module_grad > epsilon );
}

double minimize1(RnFunction f, double *x, double *grad, int n, double alpha0, double epsilon)
{
	double argmin(double alpha)
	{
		int j;
		for (j=0; j<n; j++) x[j] = x[j] - alpha * grad[j];
		double result = f(x, n);
		for (j=0; j<n; j++) x[j] = x[j] + alpha * grad[j];
		return result;
	}
	
	double a,b;
	straight_line_search_metod(argmin, alpha0, 0.01, &a, &b);
	double min = golden_section_search_min(argmin, a, b, epsilon);
	
	return min; 
}