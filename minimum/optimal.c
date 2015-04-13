#include "optimal.h"

double f0(double t, double *x, int n, double u)
{
    return x[0]*x[0] + x[1]*x[1] + u*u;
}

double f1(double t, double *x, int n, double u)
{
    return 2 * x[0];
}

double f2(double t, double *x, int n, double u)
{
    return x[0] + x[1] + u;
}

double Hamilton(double t, double *x, int n, double u, double *p)
{
    return -f0(t, x, n, u) + p[0] * f1(t, x, n, u) + p[1] * f2(t, x, n, u);
}

double dPhi1(double t, double *x, int n, double *p, double u)
{
    double _x1[n];
    double _x2[n];

    double h = 0.000001;
    _x1[0] = x[0] + h;
    _x1[1] = x[1];
    _x2[0] = x[0] - h;
    _x2[1] = x[1];

    double dH = -(Hamilton(t, _x1, n, u, p) - Hamilton(t, _x2, n, u, p)) / (2.0 * h);
    return dH;
}

double dPhi2(double t, double *x, int n, double *p, double u)
{
    double _x1[n];
    double _x2[n];

    double h = 0.000001;
    _x1[0] = x[0];
    _x1[1] = x[1] + h;
    _x2[0] = x[0];
    _x2[1] = x[1] - h;

    double dH = -(Hamilton(t, _x1, n, u, p) - Hamilton(t, _x2, n, u, p)) / (2.0 * h);
    return dH;
}

double dU(double t, double *x, int n, double u, double *p)
{
    double h = 0.000001;
    double u1 = u + h;
    double u2 = u - h;
    double dH = (Hamilton(t, x, n, u1, p) - Hamilton(t, x, n, u2, p)) / (2.0 * h);
    return dH;
}

double JSum(double *t, double *x1, double *x2, double n, double *u, int N)
{
    double sum = 0.0;
    int i=0;
    for (i=0; i<(N-1); i++)
    {
        int j=i+1;
        double f1 = x1[i]*x1[i] + x2[i]*x2[i] + u[i]*u[i];
        double f2 = x1[j]*x1[j] + x2[j]*x2[j] + u[j]*u[j];
        sum = sum + (0.5 * (f2+f1)*(t[j]-t[i]));
    }
    sum = sum + (x1[N-1]*x1[N-1] + x2[N-1]*x2[N-1]);
    return sum;
}

void control1()
{
    int i;
    int n = 2;
    double h1 = 0.1;		//step for time
    double h2 = 0.000001;	//step for runga_kutta or euler method
    double t0 = 0.0;
    double t1 = 1.0;
    double x10 = 1.0;
    double x20 = 2.0;
    int M = 10;
    int N = M + 1;

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
        u[i] = 1.0;
    }

    double j1,j2;
    do
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
            _RungaKuttaSystem2(t[i], p0, t[i-1], _p, n, h2, _x, u[i]);
            p1[i-1] = _p[0];
            p2[i-1] = _p[1];
        }
        printX("p1", p1, N);
        printX("p2", p2, N);
		
        for (i=0; i<N; i++)
        {
            double _x[] = { x1[i], x2[i] };
			double _p[] = { p1[i], p2[i] };
			gr[i] = dU(t[i], _x, n, u[i], _p);
         }
        printX("gr", gr, N);

        double argmin(double alpha)
        {
            int j;
			double u2[N];
            for (j=0; j<N; j++) u2[j] = u[j] - alpha * gr[j];
            return JSum(t, x1, x2, n, u2, N);
        }

        double a,b;
        double alpha0 = 0.0;
        straight_line_search_metod(argmin, alpha0, 0.1, &a, &b);
        double alpha = golden_section_search_min(argmin, a, b, 0.000001);
        if ( argmin(alpha) > argmin(alpha0) ) alpha = alpha0;

        printf("alpha\t%14.10f\n", alpha);

        //printf("J1\t%14.10f %14.10f\n", J1, JSum(u, N));
        j1 = JSum(t, x1, x2, n, u, N);
        for (i=0; i<N; i++)
        {
            u[i] = u[i] - alpha*gr[i];
        }
        j2 = JSum(t, x1, x2, n, u, N);
        printf("J1=%.10f J2=%.10f\n", j1, j2);
        puts("--------------------------------");
		if (fabs( j1 - j2 ) < 0.00001) break;
    } while (j1 > j2);

    free(gr);
    free(p2);
    free(p1);
    free(x2);
    free(x1);
    free(u);
}

void _RungaKuttaSystem1(double x0, double *y0, double x, double *y, int n, double h, double u)
{
    double xf1(double _t, double *_x, int _n) { return f1(_t, _x, _n, u); }
    double xf2(double _t, double *_x, int _n) { return f2(_t, _x, _n, u); }
    RmFunction xf[n];
    xf[0] = xf1;
    xf[1] = xf2;

    return RungaKuttaSystem(xf, x0, y0, x, y, n, h);
}

void _RungaKuttaSystem2(double t0, double *p0, double t, double *p, int n, double h, double *x, double u)
{
//    double __p1(double _t, double *_p, int _n) { return _p1(_t, x, _n, _p, u); }
//    double __p2(double _t, double *_p, int _n) { return _p2(_t, x, _n, _p, u); }
    double __phi1(double _t, double *_p, int _n) { return dPhi1(_t, x, _n, _p, u); }
    double __phi2(double _t, double *_p, int _n) { return dPhi2(_t, x, _n, _p, u); }
    RmFunction fp[n];
    fp[0] = __phi1;
    fp[1] = __phi2;

    return RungaKuttaSystem(fp, t0, p0, t, p, n, h);
}
