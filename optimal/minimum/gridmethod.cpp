#include "gridmethod.h"
#include <math.h>

GridMethod::GridMethod()
{
    this->alpha = 1.0;
}

void GridMethod::setLengthInterval(double x0, double x1)
{
    this->x0 = x0;
    this->x1 = x1;
}

void GridMethod::setTimeInterval(double t0, double t1)
{
    this->t0 = t0;
    this->t1 = t1;
}

void GridMethod::setLengthTimeStep(double dx, double dt)
{
    this->dx = dx;
    this->dt = dt;
    this->n = (unsigned int)(ceil(x1-x0)/dx) + 1;
    this->m = (unsigned int)(ceil(t1-t0)/dt) + 1;
}

void GridMethod::setLengthTimeStepCount(unsigned int n, unsigned int m)
{
    this->n = n;
    this->m = m;
    this->dx = (x1 - x0) / (n - 1);
    this->dt = (t1 - t0) / (m - 1);
}

void GridMethod::setF(R2Function* f)
{
    this->f = f;
}

void GridMethod::setM1(R1Function* m1)
{
    this->m1 = m1;
}

void GridMethod::setM2(R1Function* m2)
{
    this->m2 = m2;
}

void GridMethod::setFi(R1Function* fi)
{
    this->fi = fi;
}

void GridMethod::implicitDifferenceScheme()
{
    unsigned int k = n - 2;

    // initilizing

    for (unsigned int i=0; i<m; i++) u.push_back(DoubleVector(m));


    for (unsigned int j=0; j<m; j++)
    {
        u[j][0]   = m1->fx(dt * j);
        u[j][n-1] = m2->fx(dt * j);
    }

    for (unsigned int i=0; i<n; i++)
    {
        u[0][i] = fi->fx( dx * i );
    }

    /////////////////////////////////////////////////

    DoubleVector b(k);
    std::vector<DoubleVector> a;
    for (unsigned int i=0; i<k; i++)
    {
        a.push_back(DoubleVector(k));
    }

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

    for (int i=0; i<m; i++)
    {
        if (i%100==0 || i==999) {
        for (int j=0; j<n; j++)
        {
            if (j%100==0 || j==999) printf("%f ", u[i][j]);
        }
        puts("");
        }
    }
    printf("%d %d\n", m, n);
}

void GridMethod::tomas_algorithm(std::vector<DoubleVector> &a, const DoubleVector& b, DoubleVector& x)
{
    unsigned int size = x.size();

    DoubleVector p(size);
    DoubleVector q(size);

    for (int i=0; i<size; i++)
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

    for (int i=size-1; i>=0; i--)
    {
        if (i==size-1)
            x[i] = q[i];
        else
            x[i] = p[i]*x[i+1] + q[i];
    }

    p.clear();
    q.clear();
}
