#include "discreteparabolic.h"
#include <tomasmethod.h>

void DiscreteParabolic::main()
{
    DiscreteParabolic dp;

    DoubleVector f((dp.N+1)*(dp.M+1));
    for (unsigned int i=0; i<=f.size(); i++) f[i] = dp.f()
}

DiscreteParabolic::DiscreteParabolic()
{
    x0 = 0.0;
    x1 = 1.0;
    t0 = 0.0;
    t1 = 1.0;
    N = 1000;
    M = 1000;
    hx = 0.001;
    ht = 0.001;
    a = 1.0;
    lamda = 0.0;
}

DiscreteParabolic::~DiscreteParabolic()
{

}

double DiscreteParabolic::fx(const DoubleVector& f)
{
    return 0.0;
}

void DiscreteParabolic::gradient(const DoubleVector& x, DoubleVector& g, double gradient_step)
{
}

double DiscreteParabolic::fx(double x)
{
    return 0.0;
}

double DiscreteParabolic::fi(unsigned int i)
{
    return (i*hx)*(i*hx);
}

double DiscreteParabolic::m1(unsigned int j)
{
    return (j*ht)*(j*ht);
}

double DiscreteParabolic::m2(unsigned int j)
{
    return (j*ht)*(j*ht)+1.0;
}

double DiscreteParabolic::f(unsigned int i, unsigned int j)
{
    return 2.0*(j*ht) - 2.0*a;
}

void DiscreteParabolic::calculateU(const DoubleVector &f, DoubleVector &u)
{
    u.clear();
    u.resize(N+1);

    DoubleVector da;
    DoubleVector db;
    DoubleVector dc;
    DoubleVector dd;
    DoubleVector rx;

    da.resize(N-1);
    db.resize(N-1);
    dc.resize(N-1);
    dd.resize(N-1);
    rx.resize(N-1);

    double alpha = -(a*ht)/(hx*hx);
    double beta  = 1.0 + (2.0*a*ht)/(hx*hx);

    for (unsigned int j=0; j<=M; j++)
    {
        if (j == 0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u[i] = fi(i*hx);
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = alpha;
                db[i-1] = beta;
                dc[i-1] = alpha;
                //d[i-1] = u1[i] + ht * fxt(i*h, j*dt);
                dd[i-1] = u[i] + ht * f[j*(N+1)+i];
            }

            da[0]   = 0.0;
            dc[N-2] = 0.0;
            dd[0]   -= alpha * m1(j*ht);
            dd[N-2] -= alpha * m2(j*ht);

            TomasAlgorithm(da, db, dc, dd, rx);

            u[0] = m1(j*ht);
            for (unsigned int i=1; i<=N-1; i++)
            {
                u[i] = rx[i-1];
            }
            u[N] = m2(j*ht);
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();
}

