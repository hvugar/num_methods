#include "borderhyperbolic2d.h"

BorderHyperbolic2D::BorderHyperbolic2D()
{
    x10 = x20 = t0 = 0.0;
    x11 = x21 = t1 = 1.0;
    N1 = 1000;
    N2 = 1000;
    M  = 2000;
    h1 = x11 / N1;
    h2 = x21 / N2;
    ht = t1  / M;
    a1 = a2 = 1.0;
}

BorderHyperbolic2D::~BorderHyperbolic2D()
{
}

double BorderHyperbolic2D::u(unsigned int i, unsigned int j, unsigned int k) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    double t  = k*ht;
    return x1*x1 + x2*x2 + t*t;
}

double BorderHyperbolic2D::fi1(unsigned int i, unsigned int j) const
{
    return u(i, j, 0);
}

double BorderHyperbolic2D::fi2(unsigned int i, unsigned int j) const
{
    return 0.0;
}

double BorderHyperbolic2D::m1(unsigned int j, unsigned int k) const
{
    return u(0, j, k);
}

double BorderHyperbolic2D::m2(unsigned int j, unsigned int k) const
{
    return u(N1, j, k);
}

double BorderHyperbolic2D::m3(unsigned int i, unsigned int k) const
{
    return u(i, 0, k);
}

double BorderHyperbolic2D::m4(unsigned int i, unsigned int k) const
{
    return u(i, N2, k);
}

double BorderHyperbolic2D::f(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
    return 2.0 - 2.0*(a1*a1) - 2.0*(a2*a2);
}


