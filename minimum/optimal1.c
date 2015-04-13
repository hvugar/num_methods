#include "optimal.h"

void sample2()
{
    double epsilon	= 0.0001;		//dovrun sona catma meyari
    double grad_eps	= 0.005;		//gradient
    double line_eps	= 0.1;			//parcani bolme
    double gold_eps	= 0.0001;		//qizil qayda ucun
    
    int n = 11;
    int i = 0;
    double* u  = (double*) malloc( sizeof(double) * n );
    for (i=0; i<n; i++) u[i] = 0.5;
    conjugate_gradient_method(Integral, u, n, line_eps, gold_eps, grad_eps, epsilon, printer4, NULL);
    free(u);
}

double Integral(double *u, int N)
{
    int i;
    int n = 2;
    double h1 = 0.1;		//step for time
    double h2 = 0.0001;	//step for runga_kutta or euler method
    //double t0 = 0.0;
    //double t1 = 1.0;
    double x10 = 1.0;
    double x20 = 2.0;
    //int M = (t1-t0)/h1;
    int M = 10;
    //int N = M + 1;
    //printf("%d %d\n", N, M);

    double *t  = (double*) malloc( sizeof(double) * N );
    //double *u  = (double*) malloc( sizeof(double) * N );
    double *x1 = (double*) malloc( sizeof(double) * N );
    double *x2 = (double*) malloc( sizeof(double) * N );
    double *p1 = (double*) malloc( sizeof(double) * N );
    double *p2 = (double*) malloc( sizeof(double) * N );
    double *gr = (double*) malloc( sizeof(double) * N );

    for (i=0; i<N; i++)
    {
        t[i] = i*h1;
    }

    //printX("u", u, N);
    x1[0] = x10;
    x2[0] = x20;
    for (i=0; i<N-1; i++)
    {
        double x0[] = { x1[i], x2[i] };
        double _x[] = { 0.0, 0.0 };
        _RungaKuttaSystem1(t[i], x0, t[i+1], _x,  n, h2, u[i]);
        x1[i+1] = _x[0];
        x2[i+1] = _x[1];
    }
    //printX("x1", x1, N);
    //printX("x2", x2, N);
	
	double j1 = JSum(t, x1, x2, n, u, N);
	
    free(gr);
    free(p2);
    free(p1);
    free(x2);
    free(x1);
//    free(u);
	
	return j1;
}

double _p1(double t, double *x, int n, double *p, double u)
{
    return 2*x[0] - 2*p[0] - p[1];
}

double _p2(double t, double *x, int n, double *p, double u)
{
    return 2*x[1] - p[1];
}
