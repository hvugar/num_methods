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

double _JSum(Process1 *p) 
{	
}

void _printM(double **x, int m, int n)
{
    int i,j;
    for (j=0; j<m; j++)
    {
        if (j%(m/10)==0)
        {
			printf("%6d|", j);
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
        for (i=0; i<p->n; i++) p->f[j][i] = sin(j*p->dt)*cos(i*p->dx);
		for (i=0; i<p->n; i++) p->g[j][i] = 0.0;
    }
}

void _calculate_gradient(Process1 *p)
{
	
}

void _calculate()
{
    int i,j;

    Process1 p;
    init_process(&p);

	Grid g1;
	Grid g2;

	
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
	
	do
    {
		implicit_difference_scheme(fxt1, _fi1, _m11, _m21, p.alpha, p.dx, p.dt, p.x0, p.x1, p.t0, p.t1, &g1);
		p.u = (double **) malloc ( sizeof(double*) * g1.m);
		for (j=0; j<g1.m; j++)
		{
			p.u[j] = (double *) malloc( sizeof(double) * g1.n );
			memcpy(p.u[j], g1.u[j], sizeof(double) * g1.n);
			free(g1.u[j]);
		}
		free(g1.u);
		
		implicit_difference_scheme(fxt2, __fi, _m12, _m22, p.alpha, p.dx, p.dt, p.x0, p.x1, p.t0, p.t1, &g2);
		p.p = (double **) malloc ( sizeof(double*) * g2.m);
		for (j=0; j<g2.m; j++)
		{
			p.p[g2.m-1-j] = (double *) malloc( sizeof(double) * g2.n );
			memcpy(p.p[g2.m-1-j], g2.u[j], sizeof(double) * g2.n);
			free(g2.u[j]);
		}
		free(g2.u);
		
		puts("---");
		_printM(p.u, g1.m, g1.n);
		puts("---");
		_printM(p.p, g1.m, g1.n);
/*
		double grad_norm = vertor_norm(p.gradJ, p.n);
        for (i=0; i<p.n; i++) p.gradJ[i] = p.gradJ[i] / grad_norm;

        if (grad_norm < norm_epsilon)
        {
            puts("Iteration breaks. Because norm of gradient is less than given epsilon...");
            break;
        }
*/
        double argmin1(double alpha)
        {
            int i;
            double *u  = (double*) malloc( sizeof(double) * p.n );
            memcpy(u, p.u, sizeof(double)*p.n);
			for (j=0; j<p.m; j++)
			{
				
			
            for (i=0; i<p.n; i++)
            {
                p.f[j][i] = p.f[j][i] - alpha * p.p[j][i];
            }
			}
            double j = JSum(&p);
            memcpy(p.u, u, sizeof(double)*p.n);
            free(u);
            return j;
        }

        double alpha = R1Minimize(argmin1, step, gold_epsilon);

        dstnc = 0.0;
        for (i=0; i<p.n; i++)
        {
            double u = p.u[i] - alpha * p.gradJ[i];
            //			if (u < a) u = a;
            //			if (u > b) u = b;
            double du = u - p.u[i];
            p.u[i] = u;
            dstnc = dstnc + du*du;
        }
        dstnc = sqrt(dstnc);

        J2 = JSum(&p);
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
*/
	break;
    } while (1);
	
/*

    _printM(g1.u, g1.m, g1.n);
    puts("-----------------------------");
    _printM(g2.u, g2.m, g2.n);
    puts("-----------------------------");
	
	printf("%f\n", __f(0.1, 0.3));

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

    //_printM(psi, g1.m, g1.n);
*/
}
