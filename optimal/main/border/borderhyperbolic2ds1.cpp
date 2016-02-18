#include "borderhyperbolic2ds1.h"

BorderHyperbolic2DS1::BorderHyperbolic2DS1()
{
}

double BorderHyperbolic2DS1::fi1(unsigned int i, unsigned int j) const
{}

double BorderHyperbolic2DS1::fi2(unsigned int i, unsigned int j) const
{}

double BorderHyperbolic2DS1::m1(unsigned int j, unsigned int k) const
{}

double BorderHyperbolic2DS1::m2(unsigned int j, unsigned int k) const
{}

double BorderHyperbolic2DS1::m3(unsigned int i, unsigned int k) const
{}

double BorderHyperbolic2DS1::m4(unsigned int i, unsigned int k) const
{}

double BorderHyperbolic2DS1::f(unsigned int i, unsigned int j, unsigned int k) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    double t  = k*ht;
    return exp(t) + sin(x1) + cos(x2) - 2.0;
}
