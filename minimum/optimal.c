#include "optimal.h"

void _RungaKuttaSystem1(double x0, double *y0, double x, double *y, int n, double h, double u);
void _RungaKuttaSystem2(double x0, double *y0, double x1, double *y1, int n, double h, double *x);

double _f0(double t, double *x, int n, double u)
{
    return x[0]*x[0] + x[1]*x[1] + u*u;
}

double _f1(double t, double *x, int n, double u)
{
    return 2 * x[0];
}

double _f2(double t, double *x, int n, double u)
{
    return x[0] + x[1] + u;
}

double _p1(double t, double *p, int n, double *x)
{
    return 2*x[0] - 2*p[0] - p[1];
}

double _p2(double t, double *p, int n, double *x)
{
    return 2*x[1] - p[1];
}

double _H1(double t, double *x, int n, double u, double *p)
{
    return -_f0(t, x, n, u) + p[0] * _f1(t, x, n, u) + p[1] * _f2(t, x, n, u);
}

double _p11(double t, double *x, int n, double u, double *p)
{
    double x1[n];
    double x2[n];
    memcpy(x1, x, sizeof(double)*n);
    memcpy(x2, x, sizeof(double)*n);

    double h = 0.00001;
    x1[0] = x1[0] + h;
    x2[0] = x2[0] - h;

    //double d = (H(t, x1, n, u, p) - H(t, x2, n, u, p)) / (2 * h);
    return 0.0;
}

double J(double *t, double *x1, double *x2, double n, double *u, int N)
{
    double sum = 0.0;
    int i=0;
    for (i=0; i<N-1; i++)
    {
        int j=i+1;
        double f1 = x1[i]*x1[i] + x2[i]*x2[i] + u[i]*u[i];
        double f2 = x1[j]*x1[j] + x2[j]*x2[j] + u[j]*u[j];
        sum = sum + (0.5 * (f2+f1)*(t[j]-t[i]));
    }
    sum = sum + (x1[N-1]*x1[N-1] + x2[N-1]*x2[N-1]);
    return sum;
}

void sample2()
{
    int i;
    int n = 2;
    double h1 = 0.1;		//step for time
    double h2 = 0.000001;	//step for runga_kutta or euler method
    //double t0 = 0.0;
    //double t1 = 1.0;
    double x10 = 1.0;
    double x20 = 2.0;
    //int M = (t1-t0)/h1;
    int M = 10;
    int N = M + 1;
    printf("%d %d\n", N, M);

    double *t  = (double*) malloc( sizeof(double) * N );
    double *u  = (double*) malloc( sizeof(double) * N );
    double *x1 = (double*) malloc( sizeof(double) * N );
    double *x2 = (double*) malloc( sizeof(double) * N );
    double *p1 = (double*) malloc( sizeof(double) * N );
    double *p2 = (double*) malloc( sizeof(double) * N );
    double *gr = (double*) malloc( sizeof(double) * N );


    for (i=0; i<N; i++)
    {
        t[i] = i*h1;
        u[i] = 0.5;
    }

    double J1(double *u, int N)
    {
        return J(t, x1, x2, n, u, N);
    }

    double j1,j2;
    while (1)
    {
        printX("u", u, N);
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
        printX("x1", x1, N);
        printX("x2", x2, N);

        p1[N-1] = 2 * x1[N-1];
        p2[N-1] = 2 * x2[N-1];
        for (i=(N-1); i>0; i--)
        {
            double p0[] = { p1[i], p2[i] };
            double _x[] = { x1[i], x2[i] };
            double _p[] = { 0.0, 0.0 };
            _RungaKuttaSystem2(t[i], p0, t[i-1], _p,  n, h2, _x);
            p1[i-1] = _p[0];
            p2[i-1] = _p[1];
        }
        printX("p1", p1, N);
        printX("p2", p2, N);

        for (i=0; i<N; i++)
        {
            gr[i] = -2*u[i]+p2[i];
        }
        printX("gr", gr, N);

        double u2[N];
        double argmin(double alpha)
        {
            int j;
            for (j=0; j<N; j++) u2[j] = u[j] - alpha * gr[j];
            return J1(u2, N);
        }

        double a,b;
        double alpha0 = 0.0;
        straight_line_search_metod(argmin, alpha0, 0.1, &a, &b);
        double alpha = golden_section_search_min(argmin, a, b, 0.0001);

        if ( argmin(alpha) > argmin(alpha0) ) alpha = alpha0;

        printf("alpha\t%14.10f\n", alpha);

        //printf("J1\t%14.10f %14.10f\n", J1, JSum(u, N));
        j1 = J(t, x1, x2, n, u, N);
        for (i=0; i<N; i++)
        {
            u[i] = u[i] - alpha*gr[i];
        }
        j2 = J(t, x1, x2, n, u, N);
        printf("J1=%.10f J2=%.10f\n", j1, j2);
        puts("--------------------------------");
        //break;
    }

    free(gr);
    free(p2);
    free(p1);
    free(x2);
    free(x1);
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
	
	double j1 = J(t, x1, x2, n, u, N);
	
/*
    p1[N-1] = 2 * x1[N-1];
    p2[N-1] = 2 * x2[N-1];
    for (i=(N-1); i>0; i--)
    {
        double p0[] = { p1[i], p2[i] };
        double _x[] = { x1[i], x2[i] };
        double _p[] = { 0.0, 0.0 };
        _RungaKuttaSystem2(t[i], p0, t[i-1], _p,  n, h2, _x);
        p1[i-1] = _p[0];
        p2[i-1] = _p[1];
    }
    printX("p1", p1, N);
    printX("p2", p2, N);

    for (i=0; i<N; i++)
    {
        gr[i] = -2*u[i]+p2[i];
    }
    printX("gr", gr, N);

    double u2[N];
    double argmin(double alpha)
    {
        int j;
        for (j=0; j<N; j++) u2[j] = u[j] - alpha * gr[j];
        return J1(u2, N);
    }

    double a,b;
    double alpha0 = 0.0;
    straight_line_search_metod(argmin, alpha0, 0.1, &a, &b);
    double alpha = golden_section_search_min(argmin, a, b, 0.0001);

    if ( argmin(alpha) > argmin(alpha0) ) alpha = alpha0;

    printf("alpha\t%14.10f\n", alpha);

    //printf("J1\t%14.10f %14.10f\n", J1, JSum(u, N));
    j1 = J(t, x1, x2, n, u, N);
    for (i=0; i<N; i++)
    {
        u[i] = u[i] - alpha*gr[i];
    }
    j2 = J(t, x1, x2, n, u, N);
    printf("J1=%.10f J2=%.10f\n", j1, j2);
    puts("--------------------------------");
*/
    free(gr);
    free(p2);
    free(p1);
    free(x2);
    free(x1);
//    free(u);
	
	return j1;
}

void sample1()
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

void _RungaKuttaSystem1(double x0, double *y0, double x, double *y, int n, double h, double u)
{
    double __f1(double _t, double *_x, int _n) { return _f1(_t, _x, _n, u); }
    double __f2(double _t, double *_x, int _n) { return _f2(_t, _x, _n, u); }
    RmFunction f[n];
    f[0] = __f1;
    f[1] = __f2;

    return RungaKuttaSystem(f, y0, x0, y, x, n, h);
}

void _RungaKuttaSystem2(double x0, double *y0, double x1, double *y1, int n, double h, double *x)
{
    double __p1(double _t, double *_p, int _n) { return _p1(_t, _p, _n, x); }
    double __p2(double _t, double *_p, int _n) { return _p2(_t, _p, _n, x); }
    RmFunction p[n];
    p[0] = __p1;
    p[1] = __p2;

    return RungaKuttaSystem(p, y0, x0, y1, x1, n, h);
}
