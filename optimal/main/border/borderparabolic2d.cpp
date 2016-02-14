#include "borderparabolic2d.h"

S1BorderParabolic2D::S1BorderParabolic2D()
{
    x10 = x20 = t0 = 0.0;
    x11 = x21 = t1 = 1.0;
    N1 = 100;
    N2 = 100;
    M  = 100;
    h1 = x11 / N1;
    h2 = x21 / N2;
    ht = t1  / M;
    a1 = a2 = 1.0;
    count1 = 0;
}

S1BorderParabolic2D::~S1BorderParabolic2D()
{}

double S1BorderParabolic2D::u(unsigned int i, unsigned int j, unsigned int k) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    double t = 0.5*k*ht;
    return x1*x1 + x2*x2 + t*t;
}

double S1BorderParabolic2D::fi(unsigned int i, unsigned int j) const
{
    return u(i, j, 0);
}

double S1BorderParabolic2D::m1(unsigned int j, unsigned int k) const
{
    return u(0, j, k);
}

double S1BorderParabolic2D::m2(unsigned int j, unsigned int k) const
{
    return u(N1, j, k);
}

double S1BorderParabolic2D::m3(unsigned int i, unsigned int k) const
{
    return u(i, 0, k);
}

double S1BorderParabolic2D::m4(unsigned int i, unsigned int k) const
{
    return u(i, N2, k);
}

double S1BorderParabolic2D::f(unsigned int i, unsigned int j, unsigned int k) const
{
    const_cast<S1BorderParabolic2D*>(this)->count1++;
    C_UNUSED(i);
    C_UNUSED(j);
    double t = 0.5*k*ht;
    return 2.0*t - 2.0*a1 - 2.0*a2;
}

S2BorderParabolic2D::S2BorderParabolic2D()
{
    x10 = x20 = t0 = 0.0;
    x11 = x21 = t1 = 1.0;
    N1 = 100;
    N2 = 100;
    M  = 200;
    h1 = x11 / N1;
    h2 = x21 / N2;
    ht = t1  / M;
    a1 = a2 = 1.0;
    count2 = 0;
}

S2BorderParabolic2D::~S2BorderParabolic2D()
{}

double S2BorderParabolic2D::u(unsigned int i, unsigned int j, unsigned int k) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    double t  = k*ht;
    return x1*x1 + x2*x2 + t*t;
}

double S2BorderParabolic2D::fi(unsigned int i, unsigned int j) const
{
    return u(i, j, 0);
}

double S2BorderParabolic2D::m1(unsigned int j, unsigned int k) const
{
    return u(0, j, k);
}

double S2BorderParabolic2D::m2(unsigned int j, unsigned int k) const
{
    return u(N1, j, k);
}

double S2BorderParabolic2D::m3(unsigned int i, unsigned int k) const
{
    return u(i, 0, k);
}

double S2BorderParabolic2D::m4(unsigned int i, unsigned int k) const
{
    return u(i, N2, k);
}

double S2BorderParabolic2D::f(unsigned int i, unsigned int j, unsigned int k) const
{
    const_cast<S2BorderParabolic2D*>(this)->count2++;
    C_UNUSED(i);
    C_UNUSED(j);
    double t = k*ht;
    return 2.0*t - 2.0*a1 - 2.0*a2;
}
