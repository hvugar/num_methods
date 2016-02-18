#include "borderparabolic2d.h"

BorderParabolic2D::BorderParabolic2D()
{
    a1 = a2 = 1.0;
    x10 = x20 = t0 = 0.0;
    x11 = x21 = t1 = 1.0;
    h1 = 0.01;
    h2 = 0.01;
    ht = 0.01;
    N1 = (unsigned int)(ceil(x11-x10)/h1);
    N2 = (unsigned int)(ceil(x21-x20)/h2);
    M  = (unsigned int)(ceil(t1-t0)/ht);
}

BorderParabolic2D::~BorderParabolic2D()
{}

double BorderParabolic2D::u(unsigned int i, unsigned int j, unsigned int k) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    double t = 0.5*k*ht;
    return x1*x1 + x2*x2 + t*t;
}

double BorderParabolic2D::fi(unsigned int i, unsigned int j) const
{
    return u(i, j, 0);
}

double BorderParabolic2D::m1(unsigned int j, unsigned int k) const
{
    return u(0, j, k);
}

double BorderParabolic2D::m2(unsigned int j, unsigned int k) const
{
    return u(N1, j, k);
}

double BorderParabolic2D::m3(unsigned int i, unsigned int k) const
{
    return u(i, 0, k);
}

double BorderParabolic2D::m4(unsigned int i, unsigned int k) const
{
    return u(i, N2, k);
}

double BorderParabolic2D::f(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    double t = 0.5*k*ht;
    return 2.0*t - 2.0*a1 - 2.0*a2;
}
