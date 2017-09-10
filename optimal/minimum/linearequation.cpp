#include "linearequation.h"
#include <math.h>
#include <float.h>

void GaussianElimination1(const DoubleMatrix& m, const DoubleVector& b, DoubleVector &x)
{
    DoubleMatrix M = m;
    DoubleVector B = b;

    const unsigned int ui = (unsigned)0-1;

    unsigned int n = x.length();

    for (unsigned k=0; k<n; k++)
    {
        if (fabs(M.at(k,k)) <= DBL_EPSILON)
        {
            for (unsigned int p=k+1; p<n; p++)
            {
                if (fabs(M[p][k]) > DBL_EPSILON)
                {
                    M.switchRows(k, p);
                    double bk = B[k];
                    B[k] = B[p];
                    B[p] = bk;
                    break;
                }
            }
        }

        for (unsigned int j=(k+1); j<n; j++)
        {
            double c = M.at(j,k)/M.at(k,k);
            for (unsigned int i=k; i<n; i++)
            {
                M.at(j,i) = M.at(j,i) - M.at(k,i) * c;
            }
            B[j] = B[j] - B[k] *c;
        }
    }

    for (unsigned int i=(n-1); i!=ui; i--)
    {
        for (unsigned int j=(n-1); j>i; j--) B[i] -= (M.at(i,j) * x[j]);
        x[i] = B[i] / M.at(i,i);
    }
}

void LinearEquation::GaussianElimination(const DoubleMatrix& m, const DoubleVector& b, DoubleVector& x)
{
    GaussianElimination1(m,b,x);
}

void LinearEquation::FirstRowLoaded(const double *e, double f, const double *a, const double *b, const double *c, const double *d, unsigned int N)
{
    double *p = (double*) malloc(sizeof(double)*N);
    double *q = (double*) malloc(sizeof(double)*N);
    double *r = (double*) malloc(sizeof(double)*N);
    double *x = (double*) malloc(sizeof(double)*N);

    p[0] = e[0];
    q[0] = e[1];
    p[0] = f;

    for (unsigned int n=1; n<N-1; n++)
    {
        p[n] = -p[n-1]*(b[n]/a[n]) + q[n-1];
        q[n] = -p[n-1]*(c[n]/a[n]) + e[n+1];
        r[n] = -p[n-1]*(d[n]/a[n]) + r[n-1];
    }

    p[N-1] = -p[N-2]*(b[N-1]/a[N-1]) + q[N-2];
    q[N-1] = 0.0;
    r[N-1] = -p[N-2]*(d[N-1]/a[N-1]) + r[N-2];

    x[N-1] = r[N-1]/p[N-1];

    printf("%f\n", x[N-1]);

    free(r);
    free(q);
    free(p);
}
