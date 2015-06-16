#include "optimal.h"
#include "print.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct 
{
    double *t;
    double *u;
    double *x1;
    double *x2;
    double *psi1;
    double *psi2;
    double *gradJ;

    double t1;
    double t2;
    double h;
    double n;

    double x01;
    double x02;

    double *s;
} Process;

double __fx0(double t, double x1, double x2, double u)
{
    return (x1-t*t*t)*(x1-t*t*t) + (x2-t)*(x2-t) + (2*u-t)*(2*u-t);
}

double __fx1(double t, double x1, double x2, double u)
{
    return 3.0*x2*x2;
}

double __fx2(double t, double x1, double x2, double u)
{
    return x1 + x2 - 2.0*u - t*t*t + 1.0;
}

double _H(double t, double x1, double x2, double u, double psi1, double psi2)
{
    return -1.0 * __fx0(t, x1, x2, u) + psi1 * __fx1(t, x1, x2, u) + psi2 * __fx2(t, x1, x2, u);
}

double __fp1(double t, double x1, double x2, double psi1, double psi2, double u)
{
    return 2.0 * (x1 - t*t*t) - psi2;
}

double __fp2(double t, double x1, double x2, double psi1, double psi2, double u)
{
    return 2.0 * (x2 - t) - 6.0 * x2 * psi1 - psi2;
}

double gradJ(double t, double x1, double x2, double psi1, double psi2, double u)
{
    double h =  0.000001;
    return (_H(t, x1, x2, u + h, psi1, psi2) - _H(t, x1, x2, u - h, psi1, psi2)) / (2 * h);
}

double ___JSum(double *t, double *x1, double *x2, double *u, int N)
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
    sum = sum + (x2[N-1] - 1.0) * (x2[N-1] - 1.0);
    return sum;
}

void init_process(Process *p)
{
    p->t1 = 0.0;
    p->t2 = 1.0;
    p->h  = 0.001;
    p->n  = (int)(ceil((p->t2-p->t1)/p->h)) + 1;

    p->t = (double*) malloc( sizeof(double) * p->n );
    p->u = (double*) malloc( sizeof(double) * p->n );

    p->x1 = (double*) malloc( sizeof(double) * p->n );
    p->x2 = (double*) malloc( sizeof(double) * p->n );

    p->psi1 = (double*) malloc( sizeof(double) * p->n );
    p->psi2 = (double*) malloc( sizeof(double) * p->n );

    p->gradJ = (double*) malloc( sizeof(double) * p->n );
    p->s     = (double*) malloc( sizeof(double) * p->n );

    int i=0;
    for (i=0; i<p->n; i++)
    {
        p->t[i]    = i*p->h;
        p->u[i]    = 0.01;
        //        p->u[i]    = p->t[i]/2.0;
        p->x1[i]   = 0.0;
        p->x2[i]   = 0.0;
        p->psi1[i] = 0.0;
        p->psi2[i] = 0.0;
    }

    p->x01 = 0.0;
    p->x02 = 0.0;
}

void free_process(Process *p)
{
    free(p->t);
    free(p->u);
    free(p->x1);
    free(p->x2);
    free(p->psi1);
    free(p->psi2);
    free(p->gradJ);
}

void calculate_params(Process *p)
{

    double k1[] = {0.0, 0.0};
    double k2[] = {0.0, 0.0};
    double k3[] = {0.0, 0.0};
    double k4[] = {0.0, 0.0};

    double h = +fabs(p->h);
    int i = 0;
    p->x1[i] = p->x01;
    p->x2[i] = p->x02;
    for (i=0; i<p->n-1; i++)
    {
        double t  = p->t[i];
        double u  = p->u[i];
        double x1 = p->x1[i];
        double x2 = p->x2[i];

        k1[0] = __fx1(t, x1, x2, p->u[i]);
        k1[1] = __fx2(t, x1, x2, p->u[i]);
        
        x1 = p->x1[i] + (h/2.0) * k1[0];
        x2 = p->x2[i] + (h/2.0) * k1[1];
        k2[0] = __fx1(t+h/2.0, x1, x2, u);
        k2[1] = __fx2(t+h/2.0, x1, x2, u);

        x1 = p->x1[i] + (h/2.0) * k2[0];
        x2 = p->x2[i] + (h/2.0) * k2[1];
        k3[0] = __fx1(t+h/2.0, x1, x2, u);
        k3[1] = __fx2(t+h/2.0, x1, x2, u);
        
        x1 = p->x1[i] + h * k3[0];
        x2 = p->x2[i] + h * k3[1];
        k4[0] = __fx1(t+h, x1, x2, u);
        k4[1] = __fx2(t+h, x1, x2, u);
        
        p->x1[i+1] = p->x1[i] + (h/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        p->x2[i+1] = p->x2[i] + (h/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
    }


    h = -fabs(p->h);
    i=p->n-1;
    p->psi1[i] = 0.0;
    p->psi2[i] = -2.0 * (p->x2[i] - 1.0);
    for (i=p->n-1; i>0; i--)
    {
        double t = p->t[i];
        double u = p->u[i];
        double x1 = p->x1[i];
        double x2 = p->x2[i];
        double psi1 = p->psi1[i];
        double psi2 = p->psi2[i];

        k1[0] = __fp1(t, x1, x2, psi1, psi2, u);
        k1[1] = __fp2(t, x1, x2, psi1, psi2, u);

        psi1 = p->psi1[i] + (h/2.0) * k1[0];
        psi2 = p->psi2[i] + (h/2.0) * k1[1];
        k2[0] = __fp1(t+h/2.0, x1, x2, psi1, psi2, u);
        k2[1] = __fp2(t+h/2.0, x1, x2, psi1, psi2, u);

        psi1 = p->psi1[i] + (h/2.0) * k2[0];
        psi2 = p->psi2[i] + (h/2.0) * k2[1];
        k3[0] = __fp1(t+h/2.0, x1, x2, psi1, psi2, u);
        k3[1] = __fp2(t+h/2.0, x1, x2, psi1, psi2, u);

        psi1 = p->psi1[i] + h * k3[0];
        psi2 = p->psi2[i] + h * k3[1];
        k4[0] = __fp1(t+h, x1, x2, psi1, psi2, u);
        k4[1] = __fp2(t+h, x1, x2, psi1, psi2, u);

        p->psi1[i-1] = p->psi1[i] + (h/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        p->psi2[i-1] = p->psi2[i] + (h/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);

    }
}

void calculate_gradient(Process *p)
{
    int i;
    for (i=0; i<p->n; i++)
    {
        p->gradJ[i] = gradJ(p->t[i], p->x1[i], p->x2[i], p->psi1[i], p->psi2[i], p->u[i]);
    }
}

void __calculate()
{
    Process p;
    init_process(&p);

    int k = 0;
    double gr0_mod = 0.0;
    double gr1_mod = 0.0;
    double gr2_mod = 0.0;
    double dstnc = 0.0;
    int i=0;

    do
    {
        calculate_params(&p);
        calculate_gradient(&p);

        _print1("t", p.t, p.n);
        _print1("u", p.u, p.n);
        _print1("x1", p.x1, p.n);
        _print1("x2", p.x2, p.n);
        _print1("p1", p.psi1, p.n);
        _print1("p2", p.psi2, p.n);
        _print1("gr", p.gradJ, p.n);

        double J = ___JSum(p.t, p.x1, p.x2, p.u, p.n);
        printf("J = %.10f\n", J);
        //_seperator();

        // Module of gradient
        gr0_mod = 0.0;
        for (i=0; i<p.n; i++) gr0_mod = gr0_mod + p.gradJ[i]*p.gradJ[i];

        if (k == 0)
        {
            gr1_mod = gr0_mod;
            // First direction is antigradient
            for (i=0; i<p.n; i++) p.s[i] = -p.gradJ[i];
        }
        else
        {
            gr2_mod = gr0_mod;
            double w = gr2_mod / gr1_mod;
            gr1_mod = gr2_mod;
            // Direction in next (k+1) iteration
            for (i=0; i<p.n; i++) p.s[i] = -p.gradJ[i] + p.s[i] * w;
        }
		
		double grad_norm = vertor_norm(p.s, p.n);
        for (i=0; i<p.n; i++) p.s[i] = p.s[i] / grad_norm;

        double argmin1(double alpha)
        {
            int i;
            double *u1  = (double*) malloc( sizeof(double) * p.n );
            for (i=0; i<p.n; i++)
            {
                u1[i] = p.u[i] + alpha * p.s[i];
            }
            double J = ___JSum(p.t, p.x1, p.x2, u1, p.n);
            free(u1);
            return J;
        }

        double alpha = R1Minimize(argmin1, 0.01, 0.000001);

        dstnc = 0.0;
        for (i=0; i<p.n; i++)
        {
            p.u[i] = p.u[i] + alpha * p.s[i];
            dstnc = dstnc + (alpha * p.s[i]) * (alpha * p.s[i]);
        }
        dstnc = sqrt(dstnc);

        if ( k == p.n ) { k = 0; } else { k++; }

    } while (dstnc > 0.0000001);


    free_process(&p);
}
