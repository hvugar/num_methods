#include "linearequation.h"
#include <math.h>
#include <float.h>

void GaussianElimination1(const DoubleMatrix& m, const DoubleVector& b, DoubleVector &x)
{
    DoubleMatrix M = m;
    DoubleVector B = b;

    const unsigned int ui = (unsigned)0-1;

    unsigned int n = x.size();

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
