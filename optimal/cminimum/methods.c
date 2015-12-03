#include "methods.h"

double derivative1(R1Function fx, double x, double h)
{
    return (fx(x+h) - fx(x)) / h;
}

double derivative2(R1Function fx, double x, double h)
{
    return (fx(x) - fx(x-h)) / h;
}

double derivative3(R1Function fx, double x, double h)
{
    return (fx(x+h) - fx(x-h)) / (2.0*h);
}

void gradient1(RnFunction fx, double *x, double *g, unsigned int n, double h)
{
    unsigned int i;
    double tx;
    double f0 = fx(x, n);
    for (i=0; i<n; i++)
    {
        tx = x[i];
        x[i] += h;
        double f = fx(x, n);
        g[i] = (f-f0)/h;
        x[i] = tx;
    }
}

void gradient2(RnFunction fx, double *x, double *g, unsigned int n, double h)
{
    unsigned int i;
    double tx;
    double f0 = fx(x, n);
    for (i=0; i<n; i++)
    {
        tx = x[i];
        x[i] -= h;
        double f = fx(x, n);
        g[i] = (f-f0)/h;
        x[i] = tx;
    }
}

void gradient3(RnFunction fx, double *x, double *g, unsigned int n, double h)
{
    unsigned int i;
    double tx;
    for (i=0; i<n; i++)
    {
        tx = x[i];
        x[i] = tx + h;
        double f1 = fx(x, n);
        x[i] = tx - h;
        double f2 = fx(x, n);
        g[i] = (f2-f1)/(2.0*h);
        x[i] = tx;
    }
}

double trapesium1(R1Function fx, unsigned int n, double a, double b)
{
    double sum = 0.0;
    double h = (a-b)/n;
    unsigned int i;
    for (i=0; i<=n-1; i++)
    {
        double x = a + i*h;
        double f1 = fx(x);
        double f2 = fx(x+h);
        sum = sum + (f1+f2);
    }
    sum = (h/2.0)*sum;
    return sum;
}

double trapesium2(R1Function fx, double h, double a, double b)
{
    unsigned int n = (unsigned int)round((b - a)/h);
    double sum = 0.0;
    unsigned int i;
    for (i=0; i<=(n-1); i++)
    {
        double x = a + i*h;
        double f1 = fx(x);
        double f2 = fx(x+h);
        sum = sum + (f1+f2);
    }
    sum = (h/2.0)*sum;
    return sum;
}

void tomasAlgorithm(const double *a, const double *b, const double *c, const double *d, double *x, unsigned int n)
{
    double *p = (double*)malloc(sizeof(double)*n);
    double *q = (double*)malloc(sizeof(double)*n);

    const unsigned int j = (unsigned)0-1;
    unsigned int i;
    for (i=0; i<n; i++)
    {
        if (i==0)
        {
            p[0] = +d[0]/b[0];
            q[0] = -c[0]/b[0];
        }
        else
        {
            if(i==(n-1))
            {
                p[n-1] = (d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
                q[n-1] = 0.0;
            }
            else
            {
                p[i] = +(d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
                q[i] = -c[i]/(b[i]+a[i]*q[i-1]);
            }
        }
    }

    for (i=n-1; i != j; i--)
    {
        if (i==(n-1))
        {
            x[i] = p[i];
        }
        else
        {
            x[i] = p[i] + q[i]*x[i+1];
        }
    }

    free(p);
    free(q);
}
