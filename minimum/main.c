#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "methods.h"
#include "print.h"

typedef double (*RfFunction)(double t, double *x, int n, double *u, int r);
#define M_E1 2.7182818284590452353602874713527

double f0(double t, double *x, int n);
double f1(double t, double* x, int n);
double f2(double t, double* x, int n);
double pf1(double t, double *p, int n);
double pf2(double t, double *p, int n);
double HamiltonPantragen(RmFunction f0, RmFunction *f, double t, double *x, int n, double *u, int r, double *p);

void test1();
/*
double *x0 = NULL;
double *xT = NULL;
//double *x1 = NULL;
double *u  = NULL;
double *p0 = NULL;
//double *p  = NULL;

double *x1 = NULL;
double *x10 = NULL;
double *x2 = NULL;
double *x20 = NULL;
double *u0  = NULL;

int j=0;
double  *t = NULL;
double **x = NULL;
double **p = NULL;
*/
double JSum(double *u, double **x, double *t, int n);
double J(double *u, int N);

int main(int argc, char** argv)
{
    int i=0;
	
	double epsilon	= 0.01;		//dovrun sona catma meyari
	double grad_eps	= 0.005;		//gradient
	double line_eps	= 0.1;			//parcani bolme
	double gold_eps	= 0.0001;		//qizil qayda ucun
    
	int n = 11;
    double* u  = (double*) malloc( sizeof(double) * n );
	for (i=0; i<n; i++) u[i] = 0.2;
	fast_proximal_gradient_method(J, u, n, line_eps, gold_eps, grad_eps, epsilon, printer2);

	free(u);

    return 0;
    /*
    double t0 = 0.0;
    double t1 = 1.0;
    int N = 10;

    // step for time
    double h1 = (t1-t0) / N;
    // step for runga-kutta
    double h2 = 0.00001;

    int n = 2;
    //int r = 1;

    t  = (double*) malloc( sizeof(double)  * (N+1) );
    x = (double**) malloc( sizeof(double*) * (N+1) );
    p = (double**) malloc( sizeof(double*) * (N+1) );
    u  = (double*) malloc( sizeof(double)  * (N+1) );

    for (i=0; i<=N; i++)
    {
        t[i] = t0 + i*h1;
        x[i] = (double*) malloc( sizeof(double) * n );
        p[i] = (double*) malloc( sizeof(double) * n );
        u[i] = 0.9;
    }

    double *x0 = (double*) malloc( sizeof(double*) * n );
    x[0][0] = x0[0] = 1.0;
    x[0][1] = x0[1] = 2.0;
    p0 = (double*) malloc( sizeof(double) * n );

    RmFunction *f = (RmFunction*) malloc( sizeof(RmFunction*) * n );
    f[0] = f1;
    f[1] = f2;
    RmFunction *p_f = (RmFunction*) malloc( sizeof(RmFunction*) * n );
    p_f[0] = pf1;
    p_f[1] = pf2;

    //int c=0;

    double J1, J2;
    do
    {
        puts("---");
        for (i=0; i<=N; i++)
        {
            j = i;
            puts("---1");
            double t2 = t[i];
            runga_kutta2(f, x0, t0, x[i], t2, n, h2);
            puts("---2");
            printf("x1(%.2f) = %.10f x2(%.2f) = %.10f\n", t[i], x[i][0], t[i], x[i][1]);
        }

        puts("---");

        p[N][0] = p0[0] = 2*x[N][0];
        p[N][1] = p0[1] = 2*x[N][1];

        for (i=0; i<=N; i++)
        {
            j = N-i;
            runga_kutta2(p_f, p0, t1, p[j], t[j], n, h2);
            printf("p1(%.2f) = %.10f p2(%.2f) = %.10f\n", fabs(t[j]), p[j][0], fabs(t[j]), p[j][1]);
        }

        J1 = J_Sum(u, N+1);

        double* grad_J = (double*) malloc( sizeof(double) * (N+1) );

        gradient(J_Sum, u, N+1, 0.000001, grad_J);
        for (i=0; i<=N; i++)
            printf("%.6f %.6f\n", u[i], grad_J[i]);

        double *u2 = (double*) malloc(sizeof(double) * n);
        double argmin(double alpha1)
        {
            int j;
            for (j=0; j<N+1; j++) u2[j] = u[j] - alpha1 * grad_J[j];
            return J_Sum(u2, N+1);
        }

        double a,b;
        double alpha0 = 0.0;
        //double grad_eps	= 0.005;		//gradient
        double line_eps	= 0.1;			//parcani bolme
        double gold_eps	= 0.0001;		//qizil qayda ucun

        straight_line_search_metod(argmin, alpha0, line_eps, &a, &b);
        double alpha = golden_section_search_min(argmin, a, b, gold_eps);

        printf("*** %f %f\n", alpha0, alpha);

        for (i=0; i<=N; i++) u[i] = u[i] - grad_J[i]*alpha;

        free(grad_J);
        grad_J =  NULL;

        J2 = J_Sum(u, N+1);
        printf("+++ %f %f\n", J1, J2);

    } while (1);
        */

    return 0;
}

void gradient_J(double *u, double *grad_u, int n, double **p)
{
	int i=0;
	for (i=0; i<n; i++)
	{
		grad_u[i] = 2*u[i] - p[i][1];
	}
}

double J(double *u, int N)
{
    double t0 = 0.0;
    double t1 = 1.0;

    int i=0;
	int j=0;
    int n = 2;
    //int r = 1;

    // step for time
    double h1 = 0.1;
    // step for runga-kutta
    double h2 = 0.00001;

    double  *t = (double*)  malloc( sizeof(double) * N );
    double **x = (double**) malloc( sizeof(double*) * N );
    double **p = (double**) malloc( sizeof(double*) * N );

    double *x0 = (double*) malloc( sizeof(double) * n );
	double *x1 = (double*) malloc( sizeof(double) * n );
    double *p1 = (double*) malloc( sizeof(double) * n );
    double *gu = (double*) malloc( sizeof(double) * N );

    double f1(double t, double* x, int n)
    {
        return 2 * x[0];
    }

    double f2(double t, double* x, int n)
    {
        return x[0] + x[1] + u[j];
    }

    double pf1(double t, double *p, int n)
    {
        return 2*x[j][0] - p[0] - 2*p[1];
    }

    double pf2(double t, double *p, int n)
    {
        return 2*x[j][1] - 2*p[1];
    }

    RmFunction *x_f = (RmFunction*) malloc( sizeof(RmFunction*) * n );
    x_f[0] = f1;
    x_f[1] = f2;

    RmFunction *p_f = (RmFunction*) malloc( sizeof(RmFunction*) * n );
    p_f[0] = pf1;
    p_f[1] = pf2;

    x0[0] = 1.0;
    x0[1] = 2.0;
	
    for (i=0; i<N; i++)
    {
        x[i] = (double*) malloc( sizeof(double) * n );
        p[i] = (double*) malloc( sizeof(double) * n );
        t[i] = t0 + i*h1;
    }

    x[0][0] = x0[0];
    x[0][1] = x0[1];
    for (i=0; i<N; i++)
    {
        j = i;
		double _x[2];
		_x[0] = x[i][0];
		_x[1] = x[i][1];
		double _t = t[i];
        runga_kutta2(x_f, x0, t0, _x, _t, n, h2);
		x[i][0] = _x[0];
		x[i][1] = _x[1];
		runga_kutta2(x_f, x0, t0, x1  , t1, n, h2);
       // printf("t = %.2f x1[%2d] = %.10f x2[%2d] = %.10f\n", fabs(t[i]), i, x[i][0], i, x[i][1]);
       // printf("t = %.2f x1[%2d] = %.10f x2[%2d] = %.10f\n", fabs(t1), N-1, x1[0], N-1, x1[1]);
    }
	
	//puts("---");

    p[N-1][0] = p1[0] = 2*x1[0];
    p[N-1][1] = p1[1] = 2*x1[1];

	for (i=0; i<N; i++)
    {
        j = N-i-1;
        runga_kutta2(p_f, p1, t1, p[j], t[j], n, h2);
       // printf("t = %.2f p1[%2d] = %.10f p2[%2d] = %.10f\n", fabs(t[j]), j, p[j][0], j, p[j][1]);
    }
	
	double j_uk1 = JSum(u, x, t, N);
	
	for (i=0; i<n; i++)
	{
		u[i] = u[i] - 0.1*g[u];
	}
	
	double j_uk2 = JSum(u, x, t, N);
	
	//printf("%f\n", j_uk);

    for (i=0; i<N; i++)
    {
        free(x[i]);
        free(p[i]);
    }

    free(x);
    free(p);

    free(x_f);
    free(p_f);
    free(x0);
    free(p1);

    return j_uk1;
}

double JSum(double *u, double **x, double *t, int n)
{
    double sum=0.0;

    int i=0;
    for (i=0; i<n-1; i++)
    {
		int i1 = i;
		int i2 = i+1;
        double x1_1 = x[i1][0];
        double x1_2 = x[i2][0];

        double x2_1 = x[i1][1];
        double x2_2 = x[i2][1];

        double u_1 = u[i1];
        double u_2 = u[i2];

        sum += (1.0/2.0) * ((x1_1*x1_1 + x2_1*x2_1 + u_1*u_1) + (x1_2*x1_2 + x2_2*x2_2 + u_2*u_2)) * (t[i2]-t[i1]);
    }
    sum += x[n-1][0]*x[n-1][0] + x[n-1][1]*x[n-1][1];

    return sum;
}
/*
double f0(double t, double *x, int n)
{
    return x[0]*x[0] + x[1]*x[1] + u[0]*u[0];
}

double f1(double t, double* x, int n)
{
    return 2 * x[0];
}

double f2(double t, double* x, int n)
{
    return x[0] + x[1] + u[j];
}

double pf1(double t, double *p, int n)
{
    return 2*x[j][0] - p[0] - 2*p[1];
}

double pf2(double t, double *p, int n)
{
    return 2*x[j][1] - 2*p[1];
}

double HamiltonPantragen(RmFunction f0, RmFunction *f, double t, double *x, int n, double *u, int r, double *p)
{
    double a = -f0(t, x, n);
    int i=0;
    for (i=0; i<n; i++)
    {
        a = a + f[i](t, x, n) * p[i];
    }
    return a;
}
*/
void test1()
{
    /*
    int n = 2;
    int r = 1;

    x0 = (double*) malloc( sizeof(double) * n );
    xT = (double*) malloc( sizeof(double) * n );
    x1 = (double*) malloc( sizeof(double) * n );
    u  = (double*) malloc( sizeof(double) * r );
    p  = (double*) malloc( sizeof(double) * n );
    p0 = (double*) malloc( sizeof(double) * n );

    u[0] = 0.0;
    double h  = 0.00001;
    double t0 = 0.0;
    double T  = +1.0;
    double t1 = 0.1;

    x0[0] = 1.0;
    x0[1] = 2.0;

    RmFunction *f = (RmFunction*) malloc( sizeof(RmFunction*) * n );
    f[0] = f1;
    f[1] = f2;

    runga_kutta2(f, x0, t0, xT, T,  n, h);
    runga_kutta2(f, x0, t0, x1, t1, n, h);

    RmFunction *p_f = (RmFunction*) malloc( sizeof(RmFunction*) * n );
    p_f[0] = pf1;
    p_f[1] = pf2;

    p0[0]  = 2*xT[0];
    p0[1]  = 2*xT[1];
    runga_kutta2(p_f, p0, T, p, t1, n, h);

    printf("x1(1.0) = %.10f x2(1.0) = %.10f\n", xT[0], xT[1]);
    printf("x1(0.1) = %.10f x2(0.1) = %.10f\n", x1[0], x1[1]);
    printf("p1(0.1) = %.10f p2(0.1) = %.10f\n", p[0], p[1]);

    double a = HamiltonPantragen(f0, p_f, t1, x1, n, u, r, p);
    printf("%f\n", a);
*/
}
