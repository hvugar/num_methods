#include "optimal4.h"
#include "methods.h"

Process4 p;

double _RungaKutta1(R2Function f, double y0, double x0, double x, double h)
{
	int i=0;
	p.x[i] = y0;
	p.t[i] = x0;

    while (x0 <= x)
    {
        double k1 = f(x0, y0);
        double k2 = f(x0+h/2.0, y0+(h/2.0)*k1);
        double k3 = f(x0+h/2.0, y0+(h/2.0)*k2);
        double k4 = f(x0+h, y0+h*k3);

        y0 = y0 + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
        x0 = x0 + h;

		i++;
		p.x[i] = y0;
		p.t[i] = x0;
    }
	//printf("t: %f %f\n", x0, x);
    return y0;
}

double _RungaKutta2(R2Function f, double psi1, double t0, double t1, double h)
{
	int i=p.N-1;
	p.psi[i] = psi1;
	p.t[i] = t1;

    while (t0 <= t1)
    {
        double k1 = f(t1, psi1);
        double k2 = f(t1+h/2.0, psi1+(h/2.0)*k1);
        double k3 = f(t1+h/2.0, psi1+(h/2.0)*k2);
        double k4 = f(t1+h, psi1+h*k3);

        psi1 = psi1 + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
        t1 = t1 + h;

		i--;
		p.psi[i] = psi1;
    }
    return psi1;
}

void init_process(Process4 *p)
{
	p->t0 = 0.0;
	p->t1 = 1.0;

	p->x0 = 0.0;
	p->x1 = 1.7391017563;
	
	p->psi0 = 0.0;
	p->psi1 = -1.0;

	p->dx = 0.00001;
	p->dt = 0.00001;
	
	p->N = (unsigned int)(ceil((p->t1 - p->t0)/p->dt))+1;

	p->epsilon = 0.02;

	p->x   = (double*) malloc( sizeof(double)* p->N);
	p->t   = (double*) malloc( sizeof(double)* p->N);	
	p->psi = (double*) malloc( sizeof(double)* p->N);	
	
	p->T[0] = 0.2;
	p->T[1] = 0.5;

	p->p[0] = 0.3;
	p->p[1] = 0.4;
}

double delta(double t)
{
	double dlt = 0.0;
	if (fabs(t - p.T[0]) <= ((p.epsilon / 2.0) + 0.000001)) dlt = 1.0 / p.epsilon;
	if (fabs(t - p.T[1]) <= ((p.epsilon / 2.0) + 0.000001)) dlt = 1.0 / p.epsilon;
	return dlt;
}

double f(double t, double x)
{
	double d = delta(t);
	return 3*t*t+ t*t*t - x;
}

double px(double t, double psi)
{
	return psi;
}

double dx(double t, double x)
{
	return f(t, x) + p.p[0]*delta(t) + p.p[1]*delta(t);
}

void calculate_x(Process4 *p)
{
	double x1 = _RungaKutta1(dx, p->x0, p->t0, p->t1, p->dt);
}

void calculate_psi(Process4 *p)
{
	double psi0 = _RungaKutta2(px, p->psi1, p->t0, p->t1, -p->dt);
}

double J(Process4 *p)
{
	calculate_x(p);
	
	int n = p->N-1;
	return (p->x[n] - p->x1)*(p->x[n] - p->x1);
}

void calculate()
{
	init_process(&p);
		
	p.p[0] = 10;
	p.p[1] = +0.7925423850;
	
	int i = 0;
	int k = 0;
    double grad_norm = 0.0;
	double alpha = 0.0;
	
    double dstnc = 0.0;
    double step = 0.01;
    double gold_epsilon = 0.0000001;
    double dist_epsilon = 0.0000001;
    double norm_epsilon = 0.0000001;

    double J1, J2;

	do
    {
		calculate_x(&p);
		calculate_psi(&p);
		
		_print2("x:\t", p.x, p.N);
		_print2("psi:\t", p.psi, p.N);
		printf("p[0]: %.10f p[1]: %.10f\n", p.p[0], p.p[1]);

        J1 = J(&p);
        //printf("J1 = %.18f\n", J1);
		
		//printX("x", p.x, p.N);
		
		p.grad[0] = p.psi[30000];
		p.grad[1] = p.psi[40000];
	
		/* if gradinet norm at current point is less than epsilon then break. no minimize */
		grad_norm = vertor_norm(p.grad, 2);
		if (grad_norm < norm_epsilon)
        {
            puts("Iteration breaks. Because norm of gradient is less than given epsilon...");
            break;
        }
		
		k++;
		
		/* calculating unit vectors */
        for (i=0; i<2; i++) p.grad[i] = p.grad[i] / grad_norm;
		
        double argmin(double alpha)
        {
            int i;
            double *pc  = (double*) malloc(sizeof(double) * 2);
            memcpy(pc, p.p, sizeof(double) * 2);
            for (i=0; i<2; i++)
            {
                p.p[i] = p.p[i] - alpha * p.grad[i];
            }
            double sum = J(&p);
            memcpy(p.p, pc, sizeof(double) * 2);
            free(pc);
            return sum;
        }

        double alpha = R1Minimize(argmin, step, gold_epsilon);
		
        dstnc = 0.0;
        for (i=0; i<2; i++)
        {
            p.p[i] = p.p[i]- alpha * p.grad[i];
            dstnc = dstnc + (alpha * p.grad[i])*(alpha * p.grad[i]);
        }
        dstnc = sqrt(dstnc);

        J2 = J(&p);
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
    }
    while ( 1 );
	
	_print2("x:\t", p.x, p.N);
	_print2("psi:\t", p.psi, p.N);
	printf("p[0]: %.10f p[1]: %.10f\n", p.p[0], p.p[1]);
	
/*	
	for (i=0; i<p.N; i++)
	{
		if ((i%10000)==0)
		printf("%8d \tt: %0.10f \tx: %.10f \tp: %.10f\n", i, p.t[i], p.x[i], p.psi[i]);		
	}
*/
}