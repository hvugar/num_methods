#include "borderhyperbolic2d1.h"

void BorderHyperbolic21::main()
{
    DoubleMatrix m;
    BorderHyperbolic21 bh;
    bh.calculateU1(m, bh.h1, bh.h2, bh.ht, bh.N1, bh.N2, bh.M, bh.a1, bh.a2, bh.qamma);
}

BorderHyperbolic21::BorderHyperbolic21()
{
    x10 = x20 = 0.0;
    x11 = x21 = 1.0;
    t0 = 0.0;
    t1 = 10.0;

    h1 = 0.010;
    h2 = 0.010;
    ht = 0.005;
    N1 = (unsigned int)(ceil(x11-x10)/h1);
    N2 = (unsigned int)(ceil(x21-x20)/h2);
    M  = (unsigned int)(ceil(t1-t0)/ht);

    a1 = a2 = 1.0;
    qamma = 0.2;
    U = 0.0;

    alpha0 = 4.0;
    alpha1 = 5.0;
    alpha2 = 6.0;
    alpha3 = 4.0;

    e.resize(2);
    e[0] = 0.2;
    e[1] = 0.2;
}

BorderHyperbolic21::~BorderHyperbolic21()
{
}

double BorderHyperbolic21::initial1(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

double BorderHyperbolic21::initial2(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

double BorderHyperbolic21::boundary(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double BorderHyperbolic21::f(unsigned int i, unsigned int j, unsigned int k) const
{
    return fxt(i, j, k);
}

double BorderHyperbolic21::fxt(unsigned int i, unsigned int j, unsigned int k) const
{
    double sum = 0.0;
    double x1 = i*h1;
    double x2 = j*h2;
    double t = k*ht;
    sum += alpha1*exp(-alpha2*((x1-e[0])*(x1-e[0])+(x2-e[1])*(x2-e[1]))-alpha3*t);
    return sum;
}
