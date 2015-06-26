#include "optimal1.h"

double _y(double x) { return x*x*x*x; }
double _u(double x, double t) { return x*x*x*x + t*t*t - 1.0; }

double _m11(double t) { return t*t*t - 1.0; }
double _m21(double t) { return t*t*t; }
double _fi1(double x) { return x*x*x*x - 1.0 ; }
double _f(double x, double t) { return 3.0*t*t -12.0*x*x; }

double _m12(double t) { return 0.0; }
double _m22(double t) { return 0.0; }
double _fi2(double x) { return -2.0*(_u(x,1.0) - _y(x)); }
double _f2(double x, double t) { return 0.0; }

Process1 p;

double fxt1(double x, double t)
{
    int j = (int)(ceil(t/p.dt));
    int i = (int)(ceil(x/p.dx));
    return p.f[j][i];
}

double fxt2(double x, double t)
{
    return 0.0;
}

double __fi(double x)
{
    int i = (int)(ceil(x/p.dx));
    int j = (int)(ceil(p.t1/p.dt));
    return -2.0 * (p.u[j][i] - _y(x));
}

void calculate_u(Process1 *p)
{
    int j;
    Grid g;
    implicit_difference_scheme1(fxt1, _fi1, _m11, _m21, p->alpha, p->dx, p->dt, p->x0, p->x1, p->t0, p->t1, &g);
    for (j=0; j<g.m; j++)
    {
        memcpy(p->u[j], g.u[j], sizeof(double) * g.n);
        free(g.u[j]);
    }
    free(g.u);
}

void calculate_p(Process1 *p)
{
    int j;
    Grid g;
    implicit_difference_scheme1(fxt2, __fi, _m12, _m22, p->alpha, p->dx, p->dt, p->x0, p->x1, p->t0, p->t1, &g);
    for (j=0; j<g.m; j++)
    {
        memcpy(p->p[g.m-1-j], g.u[j], sizeof(double) * g.n);
        free(g.u[j]);
    }
    free(g.u);
}

double norm_f(Process1 *p)
{
    int j,i;
    double norm = 0.0;
    for (j=0; j<p->m; j++)
    {
        for (i=0; i<p->n; i++)
        {
            norm = norm + p->f[j][i]*p->f[j][i];
        }
    }
    return sqrt(norm);
}

double _JSum(Process1 *p) 
{
    calculate_u(p);

    double sum = 0.0;
    int n = p->n;
    int m = p->m;
    int i, j;

    for (i=0; i<(n-1); i++)
    {
        int i1=i+1;
        double fj = p->u[m-1][i1] - _y(p->dx*i1);
        double fi = p->u[m-1][i]; - _y(p->dx*i);
        double f0 = (fj*fj+fi*fi);
        sum = sum + 0.5 * f0 * (p->dx);
    }

    double f = 0.0;
    for (j=0; j<m; j++)
    {
        for (i=0; i<n; i++)
        {
            f += p->f[j][i]*p->f[j][i];
        }
    }
    sum = sum + f;
    return sum;
}

void init_process(Process1 *p)
{
    p->t0 = 0.0;
    p->t1 = 1.0;
    p->x0 = 0.0;
    p->x1 = 1.0;
    p->dx = 0.01;
    p->dt = 0.01;

    p->n = (int)(ceil((p->x1-p->x0)/p->dx)) + 1;
    p->m = (int)(ceil((p->t1-p->x0)/p->dt)) + 1;
    p->alpha = 1.0;

    int i,j;
    p->u = (double **) malloc( sizeof(double*) * p->m );
    p->f = (double **) malloc( sizeof(double*) * p->m );
    p->p = (double **) malloc( sizeof(double*) * p->m );
    for (j=0; j<p->m; j++)
    {
        p->u[j] = (double*) malloc( sizeof(double) * p->n );
        p->f[j] = (double*) malloc( sizeof(double) * p->n );
        p->p[j] = (double*) malloc( sizeof(double) * p->n );
        for (i=0; i<p->n; i++)
        {
            p->u[j][i] = 0.0;
            p->p[j][i] = 0.0;
            p->f[j][i] = sin(j*p->dt)*cos(i*p->dx);
        }
    }
}

void _calculate()
{
    int i,j;

    init_process(&p);

    double dstnc = 0.0;
    double step = 0.01;
    double gold_epsilon = 0.0000001;
    double dist_epsilon = 0.0000001;
    double norm_epsilon = 0.0000001;

    double J1, J2;

    do
    {
	_printM(p.f, p.m, p.n);
	_seperator();
        calculate_u(&p);
        calculate_p(&p);

        J1 = _JSum(&p);
        printf("J1 = %.18f\n", J1);


        double grad_norm = norm_f(&p);
        //for (j=0; j<p.m; j++) { for (i=0; i<p.n; i++) { p.p[j][i] = p.p[j][i] / grad_norm; } }

        if (grad_norm < norm_epsilon)
        {
            puts("Iteration breaks. Because norm of gradient is less than given epsilon...");
            break;
        }

        double argmin(double alpha)
        {
            int i,j;
            double **f  = (double**) malloc( sizeof(double*) * p.m );
            for (j=0; j<p.m; j++)
            {
                f[j] = (double*) malloc( sizeof(double) * p.n );
                memcpy(f[j], p.f[j], sizeof(double) * p.n);
                for (i=0; i<p.n; i++)
                {
                    p.f[j][i] = p.f[j][i] - alpha * p.p[j][i];
                }
            }
            double sum = _JSum(&p);
            for (j=0; j<p.m; j++)
            {
                memcpy(p.f[j], f[j], sizeof(double) * p.n);
                free(f[j]);
            }
            free(f);

            return sum;
        }

        double alpha = R1Minimize(argmin, step, gold_epsilon);

        dstnc = 0.0;
        for (j=0; j<p.m; j++)
        {
            for (i=0; i<p.n; i++)
            {
                double f = p.f[j][i] - alpha * p.p[j][i];
                double df = f - p.f[j][i];
                p.f[j][i] = f;
                dstnc = dstnc + df*df;
            }
        }
        dstnc = sqrt(dstnc);

        J2 = _JSum(&p);
        if (J2 > J1)
        {
            puts("Error accure. Then next value of functional is greater than last value...");
            exit(-1);
        }

        if (dstnc < dist_epsilon)
        {
            puts("Iteration breaks. Because distance between last and current point is less than given epsilon...");
            break;
        }

    } while (1);

    _printM(p.f, p.m, p.n);
}
