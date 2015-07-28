#include "optimal4.h"
#include "methods.h"

Process4 p;

double _RungaKutta(R2Function f, double y0, double x0, double x, double h)
{
    while (x0 <= x)
    {
        double k1 = f(x0, y0);
        double k2 = f(x0+h/2.0, y0+(h/2.0)*k1);
        double k3 = f(x0+h/2.0, y0+(h/2.0)*k2);
        double k4 = f(x0+h, y0+h*k3);

        y0 = y0 + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
        x0 = x0 + h;
		
		puts("---");
    }
    return y0;
}

void init_process(Process4 *p)
{
	p->dx = 0.01;
	p->dt = 0.01;

	p->x0 = 0.0;
	p->x1 = 1.2;
	
	p->n = (unsigned int)(ceil((p->x1-p->x0)/p->dx))+1;	
	
	p->t0 = 0.0;
	p->t1 = 0.5;
	
	p->T1 = 0.3;
	p->epsilon = 0.002;
}

double delta(double t)
{
	double d = 0.0;
	if ( fabs(t - 0.3) <= (p.epsilon / 2.0 + 0.000001))
		d = 1.0 / p.epsilon;
//	return d;
	
	if (fabs(t-0.3) <= 0.000001) return 1;
	else return 0.0;
}

double fx(double t, double x)
{
	double d = delta(t);
	printf("time: %8.6f x: %8.6f delta: %f\n", t, x, d);
	return 3*t*t+ t*t*t - x + 0.2 * d;
}

void calculate()
{
	init_process(&p);
	double x1 = _RungaKutta(fx, p.x0, p.t0, p.t1, p.dt);
	printf("%f\n", x1);
}