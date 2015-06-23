#include "optimal1.h"

double _y(double x) { return x*x; }
double _u(double x, double t) { return x*x*x*x + t*t*t - 1.0; }

double _m1(double t) { return t*t*t - 1.0; }
double _m2(double t) { return t*t*t; }
double _fi(double x) { return x*x*x*x - 1.0 ; }
double _f(double x, double t) { return 3.0*t*t -12.0*x*x; }

double _JSum() 
{	
}

void _printM(double **x, int m, int n)
{
    int i,j;
    for (j=0; j<m; j++)
    {
        if (j%(m/10)==0)
        {
            for (i=0; i<n; i++)
            {
                if (i%(n/10)==0)
                    printf("%12.6f", x[j][i]);
            }
            puts("");
        }
    }
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
    p->f = (double **) malloc( sizeof(double*) * p->m );
    p->g = (double **) malloc( sizeof(double*) * p->m );
    for (j=0; j<p->m; j++)
    {
        p->f[j] = (double*) malloc( sizeof(double) * p->n );
        p->g[j] = (double*) malloc( sizeof(double) * p->n );
        for (i=0; i<p->n; i++) p->f[j][i] = p->g[j][i] = 0.0;
    }
}

void _calculate()
{
    int i,j;

    Process1 p;
    init_process(&p);

    Grid g1;
    implicit_difference_scheme(_f, _fi, _m1, _m2, p.alpha, p.dx, p.dt, p.x0, p.x1, p.t0, p.t1, &g1);

    _printM(g1.u, g1.m, g1.n);
    puts("-----------------------------");
    _printM(p.f, p.m, p.n);
    puts("-----------------------------");
    _printM(p.g, p.m, p.n);
    puts("-----------------------------");

    double **psi = (double**) malloc( sizeof(double*) * g1.m );
    for (j=0; j<g1.m; j++)
    {
        psi[j] = (double*) malloc( sizeof(double) * g1.n );
        for (i=0; i<g1.n; i++)
        {
            psi[j][i] = 0.0;
        }
    }

    for (i=0; i<g1.n; i++)
    {
        psi[g1.m-1][i] = -2.0*(g1.u[g1.m-1][i] - _y(p.dx*i));
    }

    _printM(psi, g1.m, g1.n);
}
