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

void _calculate()
{
    int i,j;

    Process1 p;

    p.t0 = 0.0;
    p.t1 = 1.0;

    p.x0 = 0.0;
    p.x1 = 1.0;

    p.dx = 0.01;
    p.dt = 0.01;

    p.alpha = 1.0;

    Grid g1;
    implicit_difference_scheme(_f, _fi, _m1, _m2, p.alpha, p.dx, p.dt, p.x0, p.x1, p.t0, p.t1, &g1);
	
	p.f = (double **) malloc( sizeof(double*) * g1.m );
	for (j=0; j<g1.m; j++)
	{
		p.f[j] = (double*) malloc( sizeof(double) * g1.n );
		for (i=0; i<g1.n; i++)
		{
			p.f[j][i] = 0.0;
		}
	}

    printf("%d %d\n", g1.n, g1.m);

    for (j=0; j<g1.m; j++)
    {
        if (j%(g1.m/10)==0)
        {
            for (i=0; i<g1.n; i++)
            {
                if (i%(g1.n/10)==0)
                    printf("%12.6f", g1.u[j][i]);
            }
            puts("");
        }
    }
}
