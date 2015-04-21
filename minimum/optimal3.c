#include "optimal.h"

double _fx0(double t, double x1, double x2, double u)
{
    double a1 = x1-t*t*t;
    double a2 = x2-t;
    double a3 = 2*u-t;
    return  a1*a1 + a2*a2 + a3*a3;
}

double _fx1(double t, double x1, double x2, double u)
{
    return 3.0*x2*x2;
}

double _fx2(double t, double x1, double x2, double u)
{
    return x1 + x2 - 2.0*u - t*t*t + 1.0;
}

double _fp1(double t, double x1, double x2, double p1, double p2, double u)
{
    return 2.0 * (x1 - t*t*t) - p2;
}

double _fp2(double t, double x1, double x2, double p1, double p2, double u)
{
    return 2.0 * (x2 - t) - 6.0 * x2 * p1 - p2;
}

double du(double t, double x1, double x2, double p1, double p2, double u)
{
    return -4.0 * ( 2.0 * u - t ) - 2.0 * p2;
}

double _JSum(double *t, double *x1, double *x2, int n, double *u, int N)
{
    double sum = 0.0;
    int i=0;
    for (i=0; i<(N-1); i++)
    {
        int j=i+1;
        double fj = (x1[j]-t[j]*t[j]*t[j])*(x1[j]-t[j]*t[j]*t[j]) + (x2[j] - t[j])*(x2[j] - t[j]) + (2*u[j] - t[j])*(2*u[j] - t[j]);
        double fi = (x1[i]-t[i]*t[i]*t[i])*(x1[i]-t[i]*t[i]*t[i]) + (x2[i] - t[i])*(x2[i] - t[i]) + (2*u[i] - t[i])*(2*u[i] - t[i]);
        sum = sum + 0.5 * (fj+fi) * (t[j]-t[i]);
    }
    //double x[] = { x1[N-1], x2[N-1] };
    return sum;
}

void _print1(char *s, double *a, int n)
{
    int i;
    printf("double %s[] =\t{", s);
    for (i=0; i<n; i++)
    {
        if ( i%((n-1)/10) == 0 )
            printf("%12.8f", a[i]);
        if ( i%((n-1)/10) == 0 && i != n-1 )
            printf(", ");
    }
    printf("};");
    printf("\n");
}

void _seperator()
{
    puts("---------------------------------------------------------------------------------------------------------------------------------------------");
}

void calculate()
{
    double h = 0.000001;
    double t0 = 0.0;
    double t1 = 1.0;
    int N = (int)ceil((t1-t0)/h) + 1;
    int i;

    double *t  = (double*) malloc( sizeof(double) * N );
    double *u  = (double*) malloc( sizeof(double) * N );
    double *x1 = (double*) malloc( sizeof(double) * N );
    double *x2 = (double*) malloc( sizeof(double) * N );
    double *p1 = (double*) malloc( sizeof(double) * N );
    double *p2 = (double*) malloc( sizeof(double) * N );
    double *gr = (double*) malloc( sizeof(double) * N );

    for (i=0; i<N; i++)
    {
        t[i] = i*h;
        u[i] = sin(t[i]);
        x1[i] = x2[i] = p1[i] = p2[i] = 0.0;
    }
    _print1("t", t, N);
    _print1("u", u, N);
    _print1("x1", x1, N);
    _print1("x2", x2, N);
    _print1("p1", p1, N);
    _print1("p2", p2, N);

    double k11;
    double k12;
    double k21;
    double k22;
    double k31;
    double k32;
    double k41;
    double k42;

    double j1, j2;
    do
    {
        //print1("t", t, N);
        _print1("u", u, N);
        //_seperator();
        i=0;
        x1[i] = 0.0;
        x2[i] = 0.0;
        h = +0.000001;
        for (i=0; i<N-1; i++)
            //        while (i<N-1)
        {
            k11 = _fx1(t[i],       x1[i],                 x2[i],                 u[i]);
            k12 = _fx2(t[i],       x1[i],                 x2[i],                 u[i]);
            k21 = _fx1(t[i]+h/2.0, x1[i] + (h/2.0) * k11, x2[i] + (h/2.0) * k12, u[i]);
            k22 = _fx2(t[i]+h/2.0, x1[i] + (h/2.0) * k11, x2[i] + (h/2.0) * k12, u[i]);
            k31 = _fx1(t[i]+h/2.0, x1[i] + (h/2.0) * k21, x2[i] + (h/2.0) * k22, u[i]);
            k32 = _fx2(t[i]+h/2.0, x1[i] + (h/2.0) * k21, x2[i] + (h/2.0) * k22, u[i]);
            k41 = _fx1(t[i]+h,     x1[i] + (h) * k31,     x2[i] + (h) * k32,     u[i]);
            k42 = _fx2(t[i]+h,     x1[i] + (h) * k31,     x2[i] + (h) * k32,     u[i]);

            x1[i+1] = x1[i] + (h/6.0) * (k11 + 2*k21 + 2*k31 + k41);
            x2[i+1] = x2[i] + (h/6.0) * (k12 + 2*k22 + 2*k32 + k42);

            //t0 = t0 + h;
            //i++;
            //		x1[i] = x10;
            //		x2[i] = x20;
        }

        //_seperator();
        _print1("x1", x1, N);
        _print1("x2", x2, N);
        //_seperator();

        //if (t1<t0) h = -fabs(h);

        //double p10 = 0.0;
        //double p20 = -2.0 * (x2[i] - 1.0);
        //t0 = 1.0;
        //t1 = 0.0;

        i=N-1;
        p1[i] = 0.0;
        p2[i] = -2.0 * (x2[i] - 1.0);
        h = -0.000001;
        //while (fabs(t1-t0) >= fabs(h/2))
        while (i>0)
        {
            k11 = _fp1(t[i],       x1[i], x2[i], p1[i],                 p2[i],                 u[i]);
            k12 = _fp2(t[i],       x1[i], x2[i], p1[i],                 p2[i],                 u[i]);

            k21 = _fp1(t[i]+h/2.0, x1[i], x2[i], p1[i] + (h/2.0) * k11, p2[i] + (h/2.0) * k12, u[i]);
            k22 = _fp2(t[i]+h/2.0, x1[i], x2[i], p1[i] + (h/2.0) * k11, p2[i] + (h/2.0) * k12, u[i]);
            
            k31 = _fp1(t[i]+h/2.0, x1[i], x2[i], p1[i] + (h/2.0) * k21, p2[i] + (h/2.0) * k22, u[i]);
            k32 = _fp2(t[i]+h/2.0, x1[i], x2[i], p1[i] + (h/2.0) * k21, p2[i] + (h/2.0) * k22, u[i]);
            
            k41 = _fp1(t[i]+h,     x1[i], x2[i], p1[i] + (h) * k31,     p2[i] + (h) * k32,     u[i]);
            k42 = _fp2(t[i]+h,     x1[i], x2[i], p1[i] + (h) * k31,     p2[i] + (h) * k32,     u[i]);

            p1[i-1] = p1[i] + (h/6.0) * (k11 + 2*k21 + 2*k31 + k41);
            p2[i-1] = p2[i] + (h/6.0) * (k12 + 2*k22 + 2*k32 + k42);

            //t0 = t0 + h;
            i--;
            //p1[i] = p10;
            //p2[i] = p20;
        }
        _print1("p1", p1, N);
        _print1("p2", p2, N);
        //_seperator();
        //printf("%.10f\n", JSum(t, x1, x2, 2, u, N));

        for (i=0; i<N; i++)
        {
            gr[i] = du(t[i], x1[i], x2[i], p1[i], p2[i], u[i]);
        }
        _print1("gr", gr, N);

        double argmin1(double alpha)
        {
            int i;
            //double u1[N];
            double *u1  = (double*) malloc( sizeof(double) * N );
            for (i=0; i<N; i++) u1[i] = u[i] - alpha * gr[i];
            double J = _JSum(t, x1, x2, 2, u1, N);
            free(u1);
            return J;
        }

        double a,b;
        double alpha0 = 0.0;
        straight_line_search_metod(argmin1, alpha0, 0.01, &a, &b);
        double alpha = golden_section_search_min(argmin1, a, b, 0.0001);
        if ( argmin1(alpha) > argmin1(alpha0) ) alpha = alpha0;
        printf("alpha\t%14.10f\n", alpha);

        j1 = _JSum(t, x1, x2, 2, u, N) + (x2[N-1] - 1.0)*(x2[N-1] - 1.0);
        for (i=0; i<N; i++)
        {
            u[i] = u[i] - alpha*gr[i];
        }
        j2 = _JSum(t, x1, x2, 2, u, N) + (x2[N-1] - 1.0)*(x2[N-1] - 1.0);

        //double s=0.0;
        //for (i=0; i<N; i++)
        //    s += (t[i]-1.0)*(t[i]-1.0)*(t[i]-1.0)/3.0;
        printf("%.10f %.10f\n", j1, j2);
        _seperator();
    } while ((j1-j2) > 0.00000001);

    free(t);
    free(u);
    free(x1);
    free(x2);
    free(p1);
    free(p2);
    free(gr);
}
