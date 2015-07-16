#include "gridmethod.h"
#include <math.h>

GridMethod::GridMethod()
{
}

void GridMethod::implicitDifferenceScheme()
{
    double x = x1 - x0;
    double t = t1 - t0;

    unsigned int n = (unsigned int)(ceil(x/dx)) + 1;
    unsigned int m = (unsigned int)(ceil(t/dt)) + 1;
    unsigned int k = n - 2;

    // initilizing

    for (unsigned int i=0; i<m; i++) u.push_back(DoubleVector(m));

    for (unsigned int j=0; j<m; j++)
    {
        u[j][0]   = m1->fx(dt * j );
        u[j][n-1] = m2->fx( dt * j );
    }

    for (unsigned int i=0; i<n; i++)
    {
        u[0][i] = fi->fx( dx * i );
    }

    /////////////////////////////////////////////////

    DoubleVector b(k);
    std::vector<DoubleVector> a;
    for (unsigned int i=0; i<k; i++)
        a.push_back(DoubleVector(k));

    ////////////////////////////////////////////////
    double c1 = -alpha * (dt / (dx*dx));
    double c2 = 1.0 + (2.0*alpha) * (dt / (dx*dx));

    for (unsigned int j=1; j<m; j++)
    {
        for (unsigned int i=1; i<=k; i++)
        {
            if ( i == 1 )
            {
                a[i-1][0] = c2;
                a[i-1][1] = c1;
                b[i-1] = -c1*u[j][0] + u[j-1][1] + dt * f->fx(i*dx, j*dt);
            }
            else if ( i == k )
            {
                a[i-1][i-2] = c1;
                a[i-1][i-1] = c2;
                b[i-1] = -c1*u[j][k+1] + u[j-1][k] + dt * f->fx(i*dx, j*dt);
            }
            else
            {
                a[i-1][i+0] = c1;
                a[i-1][i-1] = c2;
                a[i-1][i-2] = c1;
                b[i-1] = u[j-1][i] + dt * f->fx(i*dx, j*dt);
            }
        }

        DoubleVector x(k);
        tomas_algorithm(a, b, x);
        for (unsigned int i=0; i<k; i++) u[j][i+1] = x[i];
    }

    for (unsigned int i=0; i<k; i++)
    {
        a[i].clear();
    }
    a.clear();
    b.clear();
}

void GridMethod::tomas_algorithm(std::vector<DoubleVector> &a, const DoubleVector& b, DoubleVector& x)
{
    unsigned int size = x.size();

    DoubleVector p(size);
    DoubleVector q(size);

    unsigned int i=0;
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


    p.clear();
    q.clear();
}

void GridMethod::setF(R2Function *f)
{
    this->f = f;
}
