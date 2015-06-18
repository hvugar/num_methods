#include "optimal.h"
#include "print.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double fx0(double t, double x1, double x2, double u)
{
    return (x1-t*t*t)*(x1-t*t*t) + (x2-t)*(x2-t) + (2*u-t)*(2*u-t);
    //return (x1-(t*t)/2.0)*(x1-(t*t)/2.0) + (x2-t)*(x2-t);
    //return (x1-cos(t))*(x1-cos(t)) + (x2-sin(t))*(x2-sin(t)) + (2*u-t)*(2*u-t);
    //return (x1-exp(t))*(x1-exp(t)) + (x2-exp(3*t))*(x2-exp(3*t)) + (u-t)*(u-t);
}

double T(double t, double x1, double x2, double u)
{
    return (x2 - 1.0) * (x2 - 1.0);
    //return 0.0;
    //return (x2-sin(1.0)) * (x2-sin(1.0));
    //return (x2-exp(3.0)) * (x2-exp(3.0));
}

double fx1(double t, double x1, double x2, double u)
{
    return 3.0*x2*x2;
    //return x2;
    //return -x2;
    //return x1+x2-exp(3.0*t)+u-t;
}

double fx2(double t, double x1, double x2, double u)
{
    return x1 + x2 - 2.0*u - t*t*t + 1.0;
    //return -6.0*x1 + x2 + u - t + 1.0;
    //return x1 + x2 - sin(t)- 2.0*u + t;
    //return 3.0*x1*x1*x1;
}

//-------------------------------------------------------------------------------

double H(double t, double x1, double x2, double u, double psi1, double psi2)
{
    return -1.0 * fx0(t, x1, x2, u) + psi1 * fx1(t, x1, x2, u) + psi2 * fx2(t, x1, x2, u);
}

double fp1(double t, double x1, double x2, double psi1, double psi2, double u)
{
    return -1.0 * ( H(t, x1 + dx, x2, u, psi1, psi2) - H(t, x1 - dx, x2, u, psi1, psi2) ) / (2 * dx);
    //return 2.0 * (x1 - t*t*t) - psi2;
    //return 2.0 * x1 - t * t + 6 * psi2;
    //return 2.0 * (x1-cos(t)) - psi2;
    //return 2.0 * (x1-exp(t)) - psi1 - 9.0*(x1*x1)*psi2;
}

double fp2(double t, double x1, double x2, double psi1, double psi2, double u)
{
    return -1.0 * ( H(t, x1, x2 + dx, u, psi1, psi2) - H(t, x1, x2 - dx, u, psi1, psi2) ) / (2 * dx);
    //return 2.0 * (x2 - t) - 6.0 * x2 * psi1 - psi2;
    //return 2.0 * (x2 - t) - psi1 - psi2;
    //return 2.0 * (x2-sin(t)) + psi1 - psi2;
    //return 2.0 * (x2-exp(3.0*t)) - psi1;
}

double gradJ(double t, double x1, double x2, double psi1, double psi2, double u)
{
    return ( H(t, x1, x2, u + dx, psi1, psi2) - H(t, x1, x2, u - dx, psi1, psi2) ) / (2 * dx);
    //return -4.0*(2.0*u - t) - 2.0*psi2;
    //return psi2;
    //return -4.0*(2.0*u - t) - psi2;
    //return -2.0*(u-t) + psi1;
}

double JSum(Process *p)
{
    calculate_x(p);

    int n = p->n;
    double sum = 0.0;
    int i=0;
    for (i=0; i<(n-1); i++)
    {
        int j=i+1;
        double fj = fx0(p->t[j], p->x1[j], p->x2[j], p->u[j]);
        double fi = fx0(p->t[i], p->x1[i], p->x2[i], p->u[i]);
        sum = sum + 0.5 * (fj+fi) * (p->h);
    }
    sum = sum + T(0.0, p->x1[n-1], p->x2[n-1], 0.0);
    //sum = sum + (p->x2[n-1] - 1.0) * (p->x2[n-1] - 1.0);
    //sum = sum + (x2[n-1] - sin(1.0)) * (x2[n-1] - sin(1.0));
    //sum = sum + (x2[n-1]-exp(3.0)) * (x2[n-1]-exp(3.0));
    return sum;
}

void init_process(Process *p)
{
    p->t1    = 0.0;
    p->t2    = 1.0;
    p->h     = 0.001;
    p->n     = (int)(ceil((p->t2-p->t1)/p->h)) + 1;
    
    p->x01   = 0.0;
    p->x02   = 0.0;

    p->t     = (double*) malloc( sizeof(double) * p->n );
    p->u     = (double*) malloc( sizeof(double) * p->n );

    p->x1    = (double*) malloc( sizeof(double) * p->n );
    p->x2    = (double*) malloc( sizeof(double) * p->n );

    p->psi1  = (double*) malloc( sizeof(double) * p->n );
    p->psi2  = (double*) malloc( sizeof(double) * p->n );

    p->gradJ = (double*) malloc( sizeof(double) * p->n );
    p->s     = (double*) malloc( sizeof(double) * p->n );

    int i=0;
    for (i=0; i<p->n; i++)
    {
        p->t[i]     = i*p->h;
        p->u[i]     = sin(p->t[i]);//0.02;
        p->x1[i]    = 0.0;
        p->x2[i]    = 0.0;
        p->psi1[i]  = 0.0;
        p->psi2[i]  = 0.0;
        p->gradJ[i] = 0.0;
        p->s[i]     = 0.0;
    }
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

void calculate_x(Process *p)
{
    double k1[] = {0.0, 0.0};
    double k2[] = {0.0, 0.0};
    double k3[] = {0.0, 0.0};
    double k4[] = {0.0, 0.0};

    double h = +fabs(p->h);
    int n = p->n;
    int i = 0;

    p->x1[i] = p->x01;
    p->x2[i] = p->x02;
    for (i=0; i<n-1; i++)
    {
        double t  = p->t[i];
        double u  = p->u[i];
        double x1 = p->x1[i];
        double x2 = p->x2[i];

        k1[0] = fx1(t, x1, x2, u);
        k1[1] = fx2(t, x1, x2, u);
        
        x1 = p->x1[i] + (h/2.0) * k1[0];
        x2 = p->x2[i] + (h/2.0) * k1[1];
        k2[0] = fx1(t+h/2.0, x1, x2, u);
        k2[1] = fx2(t+h/2.0, x1, x2, u);

        x1 = p->x1[i] + (h/2.0) * k2[0];
        x2 = p->x2[i] + (h/2.0) * k2[1];
        k3[0] = fx1(t+h/2.0, x1, x2, u);
        k3[1] = fx2(t+h/2.0, x1, x2, u);
        
        x1 = p->x1[i] + h * k3[0];
        x2 = p->x2[i] + h * k3[1];
        k4[0] = fx1(t+h, x1, x2, u);
        k4[1] = fx2(t+h, x1, x2, u);
        
        p->x1[i+1] = p->x1[i] + (h/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        p->x2[i+1] = p->x2[i] + (h/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
    }
}

void calculate_psi(Process *p)
{
    double k1[] = {0.0, 0.0};
    double k2[] = {0.0, 0.0};
    double k3[] = {0.0, 0.0};
    double k4[] = {0.0, 0.0};

    double h = -fabs(p->h);
    int n = p->n;
    int i=n-1;

    p->psi1[i] = -1.0*(T(0.0, p->x1[i] + dx, p->x2[i], 0.0)-T(0.0, p->x1[i] - dx, p->x2[i], 0.0))/(2.0*dx);
    p->psi2[i] = -1.0*(T(0.0, p->x1[i], p->x2[i] + dx, 0.0)-T(0.0, p->x1[i], p->x2[i] - dx, 0.0))/(2.0*dx);
    for (i=n-1; i>0; i--)
    {
        double t = p->t[i];
        double u = p->u[i];
        double x1 = p->x1[i];
        double x2 = p->x2[i];
        double psi1 = p->psi1[i];
        double psi2 = p->psi2[i];

        k1[0] = fp1(t, x1, x2, psi1, psi2, u);
        k1[1] = fp2(t, x1, x2, psi1, psi2, u);

        psi1 = p->psi1[i] + (h/2.0) * k1[0];
        psi2 = p->psi2[i] + (h/2.0) * k1[1];
        k2[0] = fp1(t+h/2.0, x1, x2, psi1, psi2, u);
        k2[1] = fp2(t+h/2.0, x1, x2, psi1, psi2, u);

        psi1 = p->psi1[i] + (h/2.0) * k2[0];
        psi2 = p->psi2[i] + (h/2.0) * k2[1];
        k3[0] = fp1(t+h/2.0, x1, x2, psi1, psi2, u);
        k3[1] = fp2(t+h/2.0, x1, x2, psi1, psi2, u);

        psi1 = p->psi1[i] + h * k3[0];
        psi2 = p->psi2[i] + h * k3[1];
        k4[0] = fp1(t+h, x1, x2, psi1, psi2, u);
        k4[1] = fp2(t+h, x1, x2, psi1, psi2, u);

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

void calculate()
{
    Process p;
    init_process(&p);

    double dstnc = 0.0;
    int i=0;

    double step = 0.01;
    double gold_epsilon = 0.0000001;
    double dist_epsilon = 0.0000001;
    double norm_epsilon = 0.0000001;

    double J1, J2;

    double a = -0.3;
    double b = +0.3;

    do
    {
        calculate_x(&p);
        calculate_psi(&p);
        calculate_gradient(&p);

        _print1("t", p.t, p.n);
        _print1("u", p.u, p.n);
        _print1("x1", p.x1, p.n);
        _print1("x2", p.x2, p.n);
        _print1("p1", p.psi1, p.n);
        _print1("p2", p.psi2, p.n);
        _print1("gr", p.gradJ, p.n);

        J1 = JSum(&p);
        printf("J1 = %.18f\n", J1);
        //_seperator();

        double grad_norm = vertor_norm(p.gradJ, p.n);
        for (i=0; i<p.n; i++) p.gradJ[i] = p.gradJ[i] / grad_norm;

        if (grad_norm < norm_epsilon)
        {
            puts("Iteration breaks. Because norm of gradient is less than given epsilon...");
            break;
        }

        double argmin1(double alpha)
        {
            int i;
            double *u  = (double*) malloc( sizeof(double) * p.n );
            memcpy(u, p.u, sizeof(double)*p.n);
            for (i=0; i<p.n; i++)
            {
                p.u[i] = p.u[i] - alpha * p.gradJ[i];

                //				if (p.u[i] < a) p.u[i] = a;
                //				if (p.u[i] > b) p.u[i] = b;
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

    } while (1);


    free_process(&p);
}
