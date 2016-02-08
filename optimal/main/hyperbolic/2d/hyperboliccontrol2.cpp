#include "hyperboliccontrol2.h"

HyperbolicControl2::HyperbolicControl2()
{
    x0 = 0.0;
    x1 = 1.0;
    t0 = 0.0;
    t1 = 1.0;
    a = 1.0;
    N = 100;
    M = 100;
    hx = (x1 - x0) / N;
    ht = (t1 - t0) / M;
    alpha0 = 1.0;
    alpha1 = 1.0;

    U0.resize(N+1);
    U1.resize(N+1);
    for (unsigned int i=0; i<=N; i++)
    {
        U0[i] = 0.0;
        U1[i] = 0.0;
    }
}

HyperbolicControl2::~HyperbolicControl2()
{
}

double HyperbolicControl2::fx(const DoubleVector &x)
{
    DoubleVector u1;
    IHyperbolicEquation::calculateU(u1, hx, ht, M, N, a);
    DoubleVector u2;
    IHyperbolicEquation::calculateU(u2, hx, ht, M-1, N, a);

    double sum = 0.0;

    double sum1 = 0.0;
    for (unsigned int i=0; i<=N; i++)
    {
        double k = 1.0;
        if (k==0 || k==N) k = 0.5;
        sum1 = sum1 + (u1[i]-U0[i])*(u1[i]-U0[i]);
    }
    sum1 = hx*sum1;

    double sum2 = 0.0;
    for (unsigned int i=0; i<=N; i++)
    {
        double k = 1.0;
        if (k==0 || k==N) k = 0.5;
        sum2 = sum2 + ((u1[i]-u2[i])/hx-U1[i])*((u1[i]-u2[i])/hx-U1[i]);
    }
    sum2 = hx*sum2;

    sum = sum1 + sum2;

    return sum;
}

void HyperbolicControl2::gradient(const DoubleVector &x, DoubleVector &g)
{
}

double HyperbolicControl2::fi1(unsigned int i) const
{
    return 2.0;
}

double HyperbolicControl2::fi2(unsigned int i) const
{
    return 0.0;
}

double HyperbolicControl2::m1(unsigned int j) const
{
    double t = j*ht;
    return t*t;
}

double HyperbolicControl2::m2(unsigned int j) const
{
    double t = j*ht;
    return t*t+1.0;
}

double HyperbolicControl2::f(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

double HyperbolicControl2::bfi1(unsigned int i) const
{
    C_UNUSED(i);
    return 0.0;
}

double HyperbolicControl2::bfi2(unsigned int i) const
{
    C_UNUSED(i);
    return 0.0;
}

double HyperbolicControl2::bm1(unsigned int j) const
{
    C_UNUSED(j);
    return 0.0;
}

double HyperbolicControl2::bm2(unsigned int j) const
{
    C_UNUSED(j);
    return 0.0;
}

double HyperbolicControl2::bf(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

double HyperbolicControl2::v1(double t) const
{
    return t;
}

double HyperbolicControl2::v2(double t) const
{
    return t;
}

double HyperbolicControl2::fxt(double x, double t) const
{
    return x+t;
}

