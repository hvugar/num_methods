#include "method_grid.h"
#include "print.h"

typedef struct {
    double *x;
    double *t;

    double x0;
    double x1;

    double t0;
    double t1;

    double dx;
    double dt;

    double alpha;
    int n;
    int m;

    double *f;
    double *u;
    double *p;
    double *g;
    double *u1;

} Process3;

double u(double x, double t) { return x*x + t*t + 2.0*x; }

double U(double x) { return x*x + 2.0*x + 1.0; }
double f(double x, double t) { return 2.0*t - 2.0; }

double fi(double x) { return x*x + 2.0*x ; }
double m1(double t) { return t*t; }
double m2(double t) { return t*t + 3.0; }

//double p_fi(double x) { return -2.0*(u(x,1.0) - y(x)); }
double p_m1(double t) { return 0.0; }
double p_m2(double t) { return 0.0; }
double p_f(double x, double t) { return 0.0; }

Process3 p;

double fxt1(double x, double t)
{
    int j = (int)(ceil(t/p.dt));
    int i = (int)(ceil(x/p.dx));
    return p.f[j*p.m+i];
}

double fxt2(double x, double t)
{
    return 0.0;
}

double p_fi(double x)
{
    int i = (int)(ceil(x/p.dx));
    int j = (int)(ceil(p.t1/p.dt));
    return -2.0 * (p.u[j*p.m+i] - U(x));
}

void calculate_u(Process3 *p)
{
    int j;
    Grid g;
    implicit_difference_scheme1(fxt1, fi, m1, m2, p->alpha, p->dx, p->dt, p->x0, p->x1, p->t0, p->t1, &g);
    for (j=0; j<g.m; j++)
    {
        memcpy(p->u+(j*g.n), g.u[j], sizeof(double) * g.n);
        free(g.u[j]);
    }
    free(g.u);
}

void calculate_p(Process3 *p)
{
    int j;
    Grid g;
    implicit_difference_scheme1(fxt2, p_fi, p_m1, p_m2, p->alpha, p->dx, p->dt, p->x0, p->x1, p->t0, p->t1, &g);
    for (j=0; j<g.m; j++)
    {
        memcpy(p->p+(g.n*(g.m-1-j)), g.u[j], sizeof(double) * g.n);
        free(g.u[j]);
    }
    free(g.u);
	
	int size = p->m * p->n;
	for (j=0; j<size; j++)
    {
		double x = p->dx*(j%p->n);
		double t = p->dt*(j/p->m);
        p->p[j] = -p->p[j] + 2.0 * (p->f[j] - f(x, t));
    }
}

double norm_f(Process3 *p)
{
    int j;
    double norm = 0.0;
    for (j=0; j<p->m*p->n; j++)
    {
        norm = norm + p->f[j]*p->f[j];
    }
    return sqrt(norm);
}

double _JSum(Process3 *p) 
{
    calculate_u(p);

    double sum = 0.0;
    int n = p->n;
    int m = p->m;
    int i, j;

    for (i=0; i<(n-1); i++)
    {
		int i0=i;
        int i1=i+1;
        double fj = p->u[(m-1)*n+i1] - U(p->dx*i1);
        double fi = p->u[(m-1)*n+i0] - U(p->dx*i0);
        double f0 = fj*fj + fi*fi;
        sum = sum + 0.5 * f0 * (p->dx);
    }

    double f_sum = 0.0;
    for (j=0; j<m*n; j++)
    {
		double x = p->dx*(j%p->n);
		double t = p->dt*(j/p->m) ;
		
        f_sum += (p->f[j]-f(x,t))*(p->f[j]-f(x,t));
    }
    sum = sum + f_sum;
    return sum;
}

void init_process(Process3 *p)
{
    p->t0 = 0.0;
    p->t1 = 1.0;
    p->x0 = 0.0;
    p->x1 = 1.0;
    p->dx = 0.01;
    p->dt = 0.01;

    p->n = (int)(ceil((p->x1-p->x0)/p->dx)) + 1;
    p->m = (int)(ceil((p->t1-p->t0)/p->dt)) + 1;
    p->alpha = 1.0;

    int j;
	int size = p->m*p->n;
    p->u  = (double*) malloc(sizeof(double) * size);
    p->f  = (double*) malloc(sizeof(double) * size);
    p->p  = (double*) malloc(sizeof(double) * size);
	p->u1 = (double*) malloc(sizeof(double) * size);
	
	p->t = (double*) malloc(sizeof(double) * p->m);
	for (j=0; j<p->m; j++) p->t[j] = p->dt*j;
	
    for (j=0; j<size; j++)
    {
        p->u[j] = 0.0;
        p->p[j] = 0.0;
        p->f[j] = 2.0;//2.0*p->t[j] - 2.0;
    }
}

void calculate()
{
    int j;

    init_process(&p);

    double dstnc = 0.0;
    double step = 0.01;
    double gold_epsilon = 0.0000001;
    double dist_epsilon = 0.0000001;
    double norm_epsilon = 0.0000001;

    double J1, J2;

    do
    {
		_seperator();
		//printX("t", p.t, p.m);
        calculate_u(&p);
        calculate_p(&p);

		//_printM(p.f, p.m, p.n);
		puts("---");
		//_printM(p.u, p.m, p.n);
		
		//for (j=0; j<p.m; j++)
		//{
		//	for (i=0; i<p.n; i++)
		//	{
		//		p.u1[j][i] = u(p.dx*i, p.dt*j);				
		//	}
		//}		
		//puts("---");
		//_printM(p.u1, p.m, p.n);

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
            int j;
            double *_f  = (double*) malloc(sizeof(double) * p.m*p.n);
            memcpy(_f, p.f, sizeof(double) * p.m*p.n);
            for (j=0; j<p.m*p.n; j++)
            {
                p.f[j] = p.f[j] - alpha * p.p[j];
            }
            double sum = _JSum(&p);
            memcpy(p.f, _f, sizeof(double) * p.m*p.n);
            free(_f);
            return sum;
        }

        double alpha = R1Minimize(argmin, step, gold_epsilon);

        dstnc = 0.0;
        for (j=0; j<p.m*p.n; j++)
        {
            double f = p.f[j]- alpha * p.p[j];
            double df = f - p.f[j];
            p.f[j] = f;
            dstnc = dstnc + df*df;
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

	_printV(p.u, p.m, p.n);
	puts("");
	_printV(p.p, p.m, p.n);
    puts("");
	_printV(p.f, p.m, p.n);
}
