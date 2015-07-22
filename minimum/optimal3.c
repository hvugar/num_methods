#include "optimal2.h"
#include "method_grid.h"
#include "print.h"
#include <math.h>

typedef struct {
	double t0;
	double t1;
	double x0;
	double x1;
	double dx;
	double dt;
	unsigned int n;
	unsigned int m;
	
	double *x;
	double *t;
	
	double epsilon;
	
	double a;
	
	double *u;
	double *U;
	double *f;
	double *p;
	
	double alpha;
} Process3;

Process3 p;

double delta(double x)
{
	double s = 0.0;
	if ( fabs(x - p.a) <= (p.epsilon / 2.0 + 0.000001))
		s = 1.0 / p.epsilon;
	return s;
}

double f1(double t) { return t*t; }
double u(double x, double t) { return x*x + 2.0*x + f1(t)*delta(x); }
double U(x) { return x*x + 2.0*x + delta(x); }
double fi(double x) { return x*x + 2.0*x; }
double m1(double t) { return f1(t)*delta(p.x0); }
double m2(double t) { return f1(t)*delta(p.x1) + 3.0; }
//double F(double x, double t) { return 2.0*t*delta(x) - 2.0; }
double F(double x, double t) { return (2.0*t - 2.0*(t*t)/((p.epsilon)*(p.epsilon))) * delta(x) - 2.0; }

double fxt(double x, double t) 
{
    int j = (int)(ceil(t/p.dt));
    int i = (int)(ceil(x/p.dx));
    return p.f[j*p.m+i];
}

void process_init(Process3 *p)
{
	p->x0 = 0.0;  p->x1 = 1.0;
	p->t0 = 0.0;  p->t1 = 1.0;
	p->dx = 0.001; p->dt = 0.001;
	p->n = 1001;   p->m  = 1001;
	p->epsilon = 0.2;
	p->alpha = 1.0;
	p->a = 0.3;
	p->x = (double*) malloc(sizeof(double) * p->n);
	p->t = (double*) malloc(sizeof(double) * p->m);
	
	int i;
	for (i=0; i<p->m; i++) p->t[i] = p->dt * i;
	for (i=0; i<p->n; i++) p->x[i] = p->dx * i;
	
	unsigned int size = p->m * p->n;
	p->u = (double*) malloc(sizeof(double) * size);
	p->f = (double*) malloc(sizeof(double) * size);	
	p->p = (double*) malloc(sizeof(double) * size);
	p->U = (double*) malloc(sizeof(double) * size);
	
	for (i=0; i<size; i++)
    {
		double x = p->dx*(i%p->n);
		double t = p->dt*(i/p->m);
        p->f[i] = F(x, t);

        p->u[i] = 0.0;
        p->p[i] = 0.0;
		
		p->U[i] = u(x, t);
    }
	_printV(p->f, p->m, p->n);
}

void calculate_u(Process3 *p)
{
    int i;
    Grid g;
    implicit_difference_scheme1(F, fi, m1, m2, p->alpha, p->dx, p->dt, p->x0, p->x1, p->t0, p->t1, &g);
    for (i=0; i<g.m; i++)
    {
        memcpy(p->u+(i*g.n), g.u[i], sizeof(double) * g.n);
        free(g.u[i]);
    }
    free(g.u);
}

void process_destroy(Process3 *p)
{
	free(p->x); free(p->t); free(p->u); free(p->f); free(p->p);
}

void _calculate()
{
	process_init(&p);
	
	calculate_u(&p);
	_seperator();
	_printV(p.u, p.m, p.n);
	_seperator();
	_printV(p.U, p.m, p.n);
	_seperator();
	
	int i;
	for (i=0; i<p.m; i++)
	{
		double t = p.dt * i;
		if (i%100==0)
		printf("%f %f %f %f %f %f\n", t, u(0.2, t), u(0.3, t), u(0.4, t), delta(0.2), delta(0.4));
	}
	
	process_destroy(&p);
}