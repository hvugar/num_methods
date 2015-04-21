#include "optimal.h"

double __fx0(double t, double* x, int n, double u)
{
    double x1 = x[0];
    double x2 = x[1];
    double a1 = x1-t*t*t;
    double a2 = x2-t;
    double a3 = 2*u-t;
    return  a1*a1 + a2*a2 + a3*a3;
}

double __fx1(double t, double *x, int n, double u)
{
    double x2 = x[1];
    return 3.0*x2*x2;
}

double __fx2(double t, double *x, int n, double u)
{
    double x1 = x[0];
    double x2 = x[1];
    return x1 + x2 - 2.0*u - t*t*t + 1.0;
}

double H(double t, double *x, int n, double u, int r, double *psi)
{
    return -__fx0(t, x, n, u) + psi[0] * __fx1(t, x, n, u) + psi[1] * __fx2(t, x, n, u);
}

double __fp1(double t, double *x, int n, double *psi, double u)
{
    double h =  0.000001;
    double x1[] = { x[0] + h, x[1] };
    double x2[] = { x[0] - h, x[1] };
    return (H(t, x1, n, u, 1, psi) - H(t, x2, n, u, 1, psi)) / (2 * h);
}

double __fp2(double t, double *x, int n, double *psi, double u)
{
    double h =  0.000001;
    double x1[] = { x[0], x[1] + h };
    double x2[] = { x[0], x[1] - h };
    return (H(t, x1, n, u, 1, psi) - H(t, x2, n, u, 1, psi)) / (2 * h);
}

double __du(double t, double *x, int n, double *psi, double u)
{
    double h =  0.000001;
    double u1 = u + h;
    double u2 = u - h;
    return (H(t, x, n, u1, 1, psi) - H(t, x, n, u2, 1, psi)) / (2 * h);
}

double __JSum(double *t, double *x1, double *x2, int n, double *u, int N);

void __calculate()
{
    double t0 = 0.0;
    double t1 = 1.0;
    double h = 0.000001;
    int N = (int)ceil((t1-t0)/h) + 1;
    int n = 2;
    //int r = 1;
    int i;

    double x0[] = { 0.0, 0.0 };

    double **x = (double**) malloc ( sizeof(double*) * n );
    x[0] = (double*) malloc( sizeof(double) * N );
    x[1] = (double*) malloc( sizeof(double) * N );

    double **p = (double**) malloc ( sizeof(double*) * n );
    p[0] = (double*) malloc( sizeof(double) * N );
    p[1] = (double*) malloc( sizeof(double) * N );

    double *u = (double*) malloc( sizeof(double) * N );
    double *t = (double*) malloc( sizeof(double) * N );

    double *gr = (double*) malloc( sizeof(double) * N );

    for (i=0; i<N; i++)
    {
        t[i] = i*h;
        u[i] = t[i]/2.0;
    }

    double k1[] = {0.0, 0.0};
    double k2[] = {0.0, 0.0};
    double k3[] = {0.0, 0.0};
    double k4[] = {0.0, 0.0};
	
    x[0][0] = x0[0];
    x[1][0] = x0[1];
	double _x[] = { x[0][i], x[1][i] };
    for (i=0; i<N-1; i++)
    {
		_x[0] = x[0][i];
		_x[1] = x[1][i];
        k1[0] = __fx1(t[i], _x, n, u[i]);
        k1[1] = __fx2(t[i], _x, n, u[i]);
		
		_x[0] = x[0][i] + (h/2.0) * k1[0];
		_x[1] = x[1][i] + (h/2.0) * k1[1];
        k2[0] = __fx1(t[i]+h/2.0, _x, n, u[i]);
        k2[1] = __fx2(t[i]+h/2.0, _x, n, u[i]);
        
		_x[0] = x[0][i] + (h/2.0) * k2[0];
		_x[1] = x[1][i] + (h/2.0) * k2[1];		
		k3[0] = __fx1(t[i]+h/2.0, _x, n, u[i]);
        k3[1] = __fx2(t[i]+h/2.0, _x, n, u[i]);
		
		_x[0] = x[0][i] + h * k3[0];
		_x[1] = x[1][i] + h * k3[1];		
        k4[0] = __fx1(t[i]+h, _x, n, u[i]);
        k4[1] = __fx2(t[i]+h, _x, n, u[i]);

        x[0][i+1] = x[0][i] + (h/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        x[1][i+1] = x[1][i] + (h/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
    }
	
	_print1("x1", x[0], N);
	_print1("x2", x[1], N);
	
	free(gr);
	free(t);
	free(u);
	free(p[0]);
	free(p[1]);
	free(p);
	free(x[0]);
	free(x[1]);
	free(x);	
}


