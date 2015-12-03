#include "methods.h"

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
