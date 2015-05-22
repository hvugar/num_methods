#include "optimal.h"

#define _USE_CONJUGATE_GRADIENT_

double fx0(double t, double* x, int n, double u)
{
    double x1 = x[0];
    double x2 = x[1];
    double a1 = x1-t*t*t;
    double a2 = x2-t;
    double a3 = 2*u-t;
    return  a1*a1 + a2*a2 + a3*a3;
}

double fx1(double t, double *x, int n, double u)
{
    double x2 = x[1];
    return 3.0*x2*x2;
}

double fx2(double t, double *x, int n, double u)
{
    double x1 = x[0];
    double x2 = x[1];
    return x1 + x2 - 2.0*u - t*t*t + 1.0;
}

double H(double t, double *x, int n, double u, int r, double *psi)
{
    return -1.0 * fx0(t, x, n, u) + psi[0] * fx1(t, x, n, u) + psi[1] * fx2(t, x, n, u);
}

double fp1(double t, double *x, int n, double *psi, double u)
{
    double h =  0.000001;
    double x1[] = { x[0] + h, x[1] };
    double x2[] = { x[0] - h, x[1] };
    return -1.0 * (H(t, x1, n, u, 1, psi) - H(t, x2, n, u, 1, psi)) / (2 * h);
    /*
    double x1 = x[0];
    double p2 = psi[1];
    return 2.0 * (x1 - t*t*t) - p2;
*/
}

double fp2(double t, double *x, int n, double *psi, double u)
{
    double h =  0.000001;
    double x1[] = { x[0], x[1] + h };
    double x2[] = { x[0], x[1] - h };
    return -1.0 * (H(t, x1, n, u, 1, psi) - H(t, x2, n, u, 1, psi)) / (2 * h);
    /*
    double x2 = x[1];
    double p1 = psi[0];
    double p2 = psi[1];
    return 2.0 * (x2 - t) - 6.0 * x2 * p1 - p2;
*/
}

double __du(double t, double *x, int n, double *psi, double u)
{
    double h =  0.000001;
    double u1 = u + h;
    double u2 = u - h;
    return (H(t, x, n, u1, 1, psi) - H(t, x, n, u2, 1, psi)) / (2.0 * h);
}

double __JSum(double *t, double **x, int n, double *u, int N)
{
    double sum = 0.0;
    int i=0;
    for (i=0; i<(N-1); i++)
    {
        int j=i+1;
        double fj = (x[0][j]-t[j]*t[j]*t[j])*(x[0][j]-t[j]*t[j]*t[j]) + (x[1][j] - t[j])*(x[1][j] - t[j]) + (2.0*u[j] - t[j])*(2.0*u[j] - t[j]);
        double fi = (x[0][i]-t[i]*t[i]*t[i])*(x[0][i]-t[i]*t[i]*t[i]) + (x[1][i] - t[i])*(x[1][i] - t[i]) + (2.0*u[i] - t[i])*(2.0*u[i] - t[i]);
        sum = sum + 0.5 * (fj+fi) * (t[j]-t[i]);
    }
    sum = sum + (x[1][N-1] - 1.0) * (x[1][N-1] - 1.0);
    return sum;
}

void __calculate()
{
    double t0 = 0.0;
    double t1 = 1.0;
    double h = 0.001;
    int N = (int)ceil((t1-t0)/h) + 1;
    int n = 2;
    //int r = 1;
    int i,j;
    double j1;
	printf("%d\n", N);

    double **x = (double**) malloc ( sizeof(double*) * n );
    x[0] = (double*) malloc( sizeof(double) * N );
    x[1] = (double*) malloc( sizeof(double) * N );

    double **p = (double**) malloc ( sizeof(double*) * n );
    p[0] = (double*) malloc( sizeof(double) * N );
    p[1] = (double*) malloc( sizeof(double) * N );

    double *u = (double*) malloc( sizeof(double) * N );
    double *u1 = (double*) malloc( sizeof(double) * N );
    double *t = (double*) malloc( sizeof(double) * N );

    double *gr = (double*) malloc( sizeof(double) * N );
    double *s  = (double*) malloc( sizeof(double) * N );

    for (i=0; i<N; i++)
    {
        t[i] = i*h;
        u[i] = sin(t[i]);
        x[0][i] = x[1][i] = p[0][i] = p[1][i] = 0.0;
    }

    j1 = 0.0;
    do
    {
        double k1[] = {0.0, 0.0};
        double k2[] = {0.0, 0.0};
        double k3[] = {0.0, 0.0};
        double k4[] = {0.0, 0.0};

        typedef double (*RR1Function)(double t, double *x, int n, double u);
        RR1Function fx[] = { fx1, fx2 };

        double x0[] = { 0.0, 0.0 };

        i = 0;
        x[0][i] = x0[0];
        x[1][i] = x0[1];
        double _x[] = { x[0][0], x[1][0] };
        h = +fabs(h);
        for (i=0; i<N-1; i++)
        {
            double k1[] = {0.0, 0.0};
            double k2[] = {0.0, 0.0};
            double k3[] = {0.0, 0.0};
            double k4[] = {0.0, 0.0};

            for (j=0; j<n; j++) _x[j] = x[j][i];
            for (j=0; j<n; j++) k1[j] = fx[j](t[i], _x, n, u[i]);
            for (j=0; j<n; j++) _x[j] = x[j][i] + (h/2.0) * k1[j];
            for (j=0; j<n; j++) k2[j] = fx[j](t[i]+h/2.0, _x, n, u[i]);
            for (j=0; j<n; j++) _x[j] = x[j][i] + (h/2.0) * k2[j];
            for (j=0; j<n; j++) k3[j] = fx[j](t[i]+h/2.0, _x, n, u[i]);
            for (j=0; j<n; j++) _x[j] = x[j][i] + h * k3[j];
            for (j=0; j<n; j++) k4[j] = fx[j](t[i]+h, _x, n, u[i]);
            for (j=0; j<n; j++) x[j][i+1] = x[j][i] + (h/6.0) * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
        }

        typedef double (*RR2Function)(double t, double *x, int n, double *psi, double uu);
        RR2Function fp[] = { fp1, fp2 };

        i=N-1;
        p[0][i] = 0.0;
        p[1][i] = -2.0 * (x[1][i] - 1.0);
        double _p[] = { p[0][i], p[1][i] };
        h = -fabs(h);
        for (i=N-1; i>0; i--)
        {
            for (j=0; j<n; j++) _x[j] = x[j][i];
            for (j=0; j<n; j++) _p[j] = p[j][i];
            for (j=0; j<n; j++) k1[j] = fp[j](t[i], _x, n, _p, u[i]);
            for (j=0; j<n; j++) _p[j] = p[j][i] + (h/2.0) * k1[j];
            for (j=0; j<n; j++) k2[j] = fp[j](t[i]+h/2.0, _x, n, _p, u[i]);
            for (j=0; j<n; j++) _p[j] = p[j][i] + (h/2.0) * k2[j];
            for (j=0; j<n; j++) k3[j] = fp[j](t[i]+h/2.0, _x, n, _p, u[i]);
            for (j=0; j<n; j++) _p[j] = p[j][i] + h * k3[j];
            for (j=0; j<n; j++) k4[j] = fp[j](t[i]+h, _x, n, _p, u[i]);
            for (j=0; j<n; j++) p[j][i-1] = p[j][i] + (h/6.0) * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
        }
        j1 = __JSum(t, x, n, u, N);
		
        for (i=0; i<N; i++)
        {
            _x[0] = x[0][i];
            _x[1] = x[1][i];
            _p[0] = p[0][i];
            _p[1] = p[1][i];
            gr[i] = __du(t[i], _x, n, _p, u[i]);
        }
        //_print1("u", u, N);
        //_print1("x1", x[0], N);
        //_print1("x2", x[1], N);
        //_print1("p1", p[0], N);
        //_print1("p2", p[1], N);
        //_print1("gr", gr, N);
        printf("J(u[k])    = %.10f\n",j1);
		//_seperator();
		
		double argmin(double alpha)
        {
            int i;
            double *u2  = (double*) malloc( sizeof(double) * N );
			for (i=0; i<N; i++) { u2[i] = u[i] - alpha * gr[i]; }
            double J = __JSum(t, x, n, u2, N);
            free(u2);
            return J;
        }
		
		double alpha = R1Minimize(argmin, 0.01, 0.000001);
		
		memcpy(u1, u, sizeof(double) * N);
		for (i=0; i<N; i++) { u[i] = u[i] - alpha * gr[i]; }
    } while ( distance(u1, u, N) > 0.0000001 );

    free(gr);
    free(s);
    free(t);
    free(u);
    free(u1);
    free(p[0]);
    free(p[1]);
    free(p);
    free(x[0]);
    free(x[1]);
    free(x);
}

void RungaKuttaSystem1(double t0, const double *x0, double t1, double *x1, const int n, double h, double u)
{
    double ____fx1(double _t, double *_x, int _n) { return fx1(_t, _x, _n, u); }
    double ____fx2(double _t, double *_x, int _n) { return fx2(_t, _x, _n, u); }
    RmFunction fx[] = { ____fx1, ____fx2 };

    RungaKuttaSystem(fx, t0, x0, t1, x1, n, h);
}
