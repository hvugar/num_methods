#include "borderhyperbolic2d.h"

#define SAMPLE5

BorderHyperbolic2D::BorderHyperbolic2D()
{
    a1 = a2 = 1.0;
    x10 = x20 = t0 = 0.0;
    x11 = x21 = t1 = 1.0;
    h1 = 0.01;
    h2 = 0.01;
    ht = 0.0005;
    N1 = (unsigned int)(ceil(x11-x10)/h1);
    N2 = (unsigned int)(ceil(x21-x20)/h2);
    M  = (unsigned int)(ceil(t1-t0)/ht);

    qamma = 0.2;
}

BorderHyperbolic2D::~BorderHyperbolic2D()
{
}

double BorderHyperbolic2D::u(unsigned int i, unsigned int j, unsigned int k) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    double t  = k*ht;
#ifdef SAMPLE1
    return x1*x1 + x2*x2 + t*t;
#endif
#ifdef SAMPLE2
    return x1*x1 + x1*x2 + sin(x1) + cos(x2) + exp(t);
#endif
#ifdef SAMPLE3
    return sin(x1)*sin(x1) + cos(x2) + x1*x2 + t*t*t;
#endif
#ifdef SAMPLE4
    return sin(t) + exp(x1) + x2*x2;
#endif
#ifdef SAMPLE5
    return x1*x1 + x1*x2 + 5.0*sin(10.0*x1) + cos(5.0*x2) + 4.0*cos(5.0*x2*t);
#endif
}

double BorderHyperbolic2D::fi1(unsigned int i, unsigned int j) const
{
    return u(i, j, 0);
}

double BorderHyperbolic2D::fi2(unsigned int i, unsigned int j) const
{
#ifdef SAMPLE1
    return 0.0;
#endif
#ifdef SAMPLE2
    return 1.0;
#endif
#ifdef SAMPLE3
    return 0.0;
#endif
#ifdef SAMPLE4
    return 1.0;
#endif
#ifdef SAMPLE5
    return 0.0;
#endif
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
    double x1 = i*h1;
    double x2 = j*h2;
    double t  = k*ht;
    C_UNUSED(x1);
    C_UNUSED(x2);
    C_UNUSED(t);
#ifdef SAMPLE1
    return 2.0 - 2.0*(a1*a1) - 2.0*(a2*a2);
#endif
#ifdef SAMPLE2
    return exp(t) + (a1*a1)*(sin(x1) - 2.0) + (a2*a2)*cos(x2);// + qamma*exp(t);
#endif
#ifdef SAMPLE3
    return 6*t - (a1*a1)*(2.0*cos(2.0*x1)) + (a2*a2)*cos(x2) ;
#endif
#ifdef SAMPLE4
    return -sin(t) - (a1*a1)*exp(x1) - (a2*a2)*2.0;
#endif
#ifdef SAMPLE5
    return -100.0*x1*x1*cos(5.0*x1*t) - a1*a1*(2.0 - 500.0*sin(10.0*x1) - 100.0*t*t*cos(5.0*x1*t)) - a2*a2*(-75.0*cos(5.0*x2));
#endif
}


