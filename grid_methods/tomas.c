#include <math.h>
#include <stdlib.h>

void tomas_algorithm(double *a[], double* b, double* x, int size)
{
    double* p = (double*) malloc(sizeof(double)*size);
    double* q = (double*) malloc(sizeof(double)*size);

    int i=0;
    for (i=0; i<size; i++)
    {
        if (i==0)
        {
            p[i] = -a[i][i+1] / a[i][i];
            q[i] = b[i] / a[i][i];
        }
        else if (i==size-1)
        {
            p[i] = 0.0;
            q[i] = (b[i] - a[i][i-1]*q[i-1])/(a[i][i]+a[i][i-1]*p[i-1]);
        }
        else
        {
            p[i] = -(a[i][i+1])/(a[i][i]+a[i][i-1]*p[i-1]);
            q[i] = (b[i] - a[i][i-1]*q[i-1])/(a[i][i]+a[i][i-1]*p[i-1]);
        }
    }

    for (i=size-1; i>=0; i--)
    {
        if (i==size-1)
            x[i] = q[i];
        else
            x[i] = p[i]*x[i+1] + q[i];
    }

    free(p);
	p = NULL;
    free(q);
	q = NULL;
}