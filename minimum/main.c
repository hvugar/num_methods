#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "methods.h"
#include "print.h"
#include "optimal.h"

int N = 0;
int n = 2;
double h2 = 0.0;
double *x1 = NULL;
double *x2 = NULL;
double *t  = NULL;
double *u  = NULL;

void printX2(char *label, double *x, int n)
{
	int m = n / 10;
	int i;
	printf("double %s[] = \t{", label);
	for (i=0; i<n; i++)
	{
		if (i%m==0)
			printf("%12.8f", x[i]);
		//if (i != n-1 )
		//	printf(", ");
	}
	printf("};");
	printf("\n");
}


double F(double *x, int n)
{
    return (x[1] - 1.0) * (x[1] - 1.0);
}

double f0(double t, double *x, int n, double u)
{
    return (x[0]-t*t*t)*(x[0]-t*t*t) + x[1]*x[1] - t*t + (2*u-1.0)*(2*u-1.0);
}

double f1(double t, double *x, int n, double u)
{
    return 3.0*x[1]*x[1];
}

double f2(double t, double *x, int n, double u)
{
    return x[0] + x[1] - 2.0*u - t*t*t + 1.0;
}

double JSum(double *u, int N)
{
    //printX2("u", u, N);
    //printX2("t", t, N);
    x1[0] = 0.0;
    x2[0] = 0.0;
    {
        double *k1 = (double*) malloc( sizeof(double) * n );
        double *k2 = (double*) malloc( sizeof(double) * n );
        double *k3 = (double*) malloc( sizeof(double) * n );
        double *k4 = (double*) malloc( sizeof(double) * n );
        double *xc = (double*) malloc( sizeof(double) * n );

        double _t0 = 0.0;
        double _t1 = 1.0;
        int i=0;
        while ((_t1-_t0) >= h2)
        {
            double x[] = {x1[i], x2[i]};

            // Calculating k1 vector
            k1[0] = smp1_f1(_t0, x, n, u[i]);
            k1[1] = smp1_f2(_t0, x, n, u[i]);

            // Calculating k2 vector
            xc[0] = x[0] + (h2/2.0) * k1[0];
            xc[1] = x[1] + (h2/2.0) * k1[1];

            k2[0] = smp1_f1(_t0+h2/2.0, xc, n, u[i]);
            k2[1] = smp1_f2(_t0+h2/2.0, xc, n, u[i]);

            // Calculating k3 vector
            xc[0] = x[0] + (h2/2.0) * k2[0];
            xc[1] = x[1] + (h2/2.0) * k2[1];

            k3[0] = smp1_f1(_t0+h2/2.0, xc, n, u[i]);
            k3[1] = smp1_f2(_t0+h2/2.0, xc, n, u[i]);

            // Calculating k1 vector
            xc[0] = x[0] + h2 * k3[0];
            xc[1] = x[1] + h2 * k3[1];

            k4[0] = smp1_f1(_t0+h2, xc, n, u[i]);
            k4[1] = smp1_f2(_t0+h2, xc, n, u[i]);

            // Calculating y
            x1[i+1] = x1[i] + (h2/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
            x2[i+1] = x2[i] + (h2/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);

            _t0 = _t0 + h2;
            i++;
        }

        free(xc);
        free(k1);
        free(k2);
        free(k3);
        free(k4);
    }

    double sum = 0.0;
    int i=0;
    for (i=0; i<(N-1); i++)
    {
        int j=i+1;
        double fj = (x1[j]-t[j]*t[j]*t[j])*(x1[j]-t[j]*t[j]*t[j]) + x2[j]*x2[j] - t[j]*t[j] + (2*u[j] - 1.0)*(2*u[j] - 1.0);
        double fi = (x1[i]-t[i]*t[i]*t[i])*(x1[i]-t[i]*t[i]*t[i]) + x2[i]*x2[i] - t[i]*t[i] + (2*u[i] - 1.0)*(2*u[i] - 1.0);
        sum = sum + 0.5 * (fj+fi) * (t[j]-t[i]);
    }
    double x[] = { x1[N-1], x2[N-1] };
    sum = sum + smp1_F(x, n);
    return sum;
}



int main(int argc, char** argv)
{
/*
    N = 1001;
    h2 = 0.001;
    x1 = (double*) malloc( sizeof(double) * N );
    x2 = (double*) malloc( sizeof(double) * N );
    t  = (double*) malloc( sizeof(double) * N );
    u  = (double*) malloc( sizeof(double) * N );

    double epsilon	= 0.00001;		//dovrun sona catma meyari
    double grad_eps	= 0.0001;		//gradient
    double line_eps	= 0.01;			//parcani bolme
    double gold_eps	= 0.00001;		//qizil qayda ucun

    int i=0;
    for (i=0; i<N; i++)
    {
        t[i] = i*h2;
        u[i] = sin(t[i]);
    }
    
    conjugate_gradient_method(JSum, u, N, line_eps, gold_eps, grad_eps, epsilon, NULL, NULL);
	
	puts("----");
	for (i=0; i<N; i++)
	{
		printf("%9.5f", u[i]);
		if ((i+1)%20==0)
		puts("");
	}
	printX2("u", u, N);
	puts("-----");

    free(u);
*/
    smp1_control();
    return 0;
}
