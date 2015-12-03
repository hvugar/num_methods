#include "tomasmethod.h"

void TomasAlgorithm(const DoubleVector &a, const DoubleVector &b, const DoubleVector &c, const DoubleVector &d, DoubleVector &x)
{
    if (x.size() != a.size() || x.size() != b.size() || x.size() != c.size() || x.size() != d.size()) return;

    int n = x.size();

    DoubleVector p(n);
    DoubleVector q(n);

    for (int i=0; i<n; i++)
    {
        if (i==0)
        {
            p[0] = d[0]/b[0];
            q[0] = -c[0]/b[0];
        } else
        {
            if(i==n-1)
            {
                p[n-1] = (d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
                q[n-1] = 0.0;
            }
            else
            {
                p[i] = (d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
                q[i] = -c[i]/(b[i]+a[i]*q[i-1]);
            }
        }
    }

    for (int i=n-1; i>=0; i--)
    {
        if (i==n-1)
        {
            x[i] = p[i];
        }
        else
        {
            x[i] = p[i] + q[i]*x[i+1];
        }
    }

    p.clear();
    q.clear();
}


void tomasAlgorithm(const double *a, const double *b, const double *c, const double *d, double *x, unsigned int n)
{
    double *p = (double*)malloc(sizeof(double)*n);
    double *q = (double*)malloc(sizeof(double)*n);

    for (unsigned int i=0; i<n; i++)
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

    for (unsigned int i=n-1, j=((unsigned)0-1); i != j; i--)
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
