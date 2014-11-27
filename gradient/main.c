#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "minimum.h"
#include "gradient.h"

double f(double *x, int n);
double g(double alpha);
double g1(double alpha);
void fast_proximal_gradient_method();

double  *x;
int     N;
double  epsilon;
double  h;

double **S;
double **Xj;
int j=0;
int k=0;

int main(int argc, char** argv)
{
    epsilon = 0.001;
    h       = 0.0001;
    
	N = 2;
    x  = (double*) malloc( sizeof(double) * N );

    x[0]    = -0.5;
    x[1]    = -1.0;
	
	fast_proximal_gradient_method();

	return 0;
}

//Метод наискорейшего спуска
//Fast proximal gradient method
void fast_proximal_gradient_method()
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

double f(double *x, int n)
{
    return pow(x[0],3) + 2*pow(x[1],2) - 3*x[0] - 4*x[1];
}

double g(double alpha)
{
    double* _x = (double*) malloc( sizeof(double) * N );

    double* gr = (double*) malloc(sizeof(double) * N);
    gradient(f, x, N, h, gr);

	int i;
    for (i=0; i<N; i++)
    {
        _x[i] = x[i] - alpha * gr[i];
    }

    double result = f(_x, N);

    free(gr);
    free(_x);

    return result;
}

double g1(double alpha)
{
    double* gr = (double*) malloc( sizeof(double) * N );
    gradient(f, x, N, h, gr);

    double* _x = (double*) malloc( sizeof(double) * N );
	
	int i;
    for (i=0; i<N; i++)
    {
        _x[i] = x[i] - alpha * gr[i];
    }
    double result = f(_x, N);

    free(gr);
    free(_x);

    return result;
}