#include "method_grid.h"

void tomas_algorithm(double **a, double* b, double* x, int size)
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

void implicit_difference_scheme(R2Function f, R1Function fi, R1Function m1, R1Function m2, double alpha, double dx, double dt, double x0, double x1, double t0, double t1, grid *g)
{
    double x = x1 - x0;
    double t = t1 - t0;

    int n    = (int)(ceil(x/dx)) + 1;
    int m    = (int)(ceil(t/dt)) + 1;
    int k    = n - 2;

    g->u = (double**) malloc( sizeof(double*) * m );
    g->n = n;
    g->m = m;

    double **a = (double**) malloc( sizeof(double*) * k );
    double  *b = (double*)  malloc( sizeof(double)  * k );

    int i,j;
    for (j=0; j<m; j++)
    {
        g->u[j] = (double*) malloc( sizeof(double) * n );

        for (i=0; i<n; i++)
        {
            g->u[j][i] = 0.0;
            if (j==0)   g->u[j][i] = fi( dx * i );
            if (i==0)   g->u[j][i] = m1( dt * j );
            if (i==n-1) g->u[j][i] = m2( dt * j );
        }
    }

    for (i=0; i<k; i++)
    {
        a[i] = (double*)  malloc( sizeof(double)  * k );
        for (j=0; j<k; j++) a[i][j] = 0.0;
    }

    double c1 = -alpha * (dt / (dx*dx));
    double c2 = 1.0 + (2.0*alpha) * (dt / (dx*dx));

    for (j=1; j<m; j++)
    {
        for (i=1; i<=k; i++)
        {
            b[i-1] = 0.0;
            b[i-1] = g->u[j-1][i] + dt * f(i*dx, j*dt);
            if ( i==1 ) b[i-1] = -c1*g->u[j][0] + g->u[j-1][1] + dt * f(i*dx, j*dt);
            if ( i==k ) b[i-1] = -c1*g->u[j][k+1] + g->u[j-1][k] + dt * f(i*dx, j*dt);
        }

        for (i=0; i<k; i++)
        {
            if ( i == 0 )
            {
                a[i][0] = c2;
                a[i][1] = c1;
            }
            else if ( i == (k-1) )
            {
                a[i][i-1] = c1;
                a[i][i+0] = c2;
            }
            else
            {
                a[i][i+1] = -dt/(dx*dx);
                a[i][i+0] = 1 + (2*dt)/(dx*dx);
                a[i][i-1] = -dt/(dx*dx);
            }
        }

        tomas_algorithm(a, b, (g->u[j]+1), k);
    }

    for (i=0; i<k; i++)
    {
        free(a[i]);
    }
    free(a);
    free(b);
}
