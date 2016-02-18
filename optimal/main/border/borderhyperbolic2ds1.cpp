#include "borderhyperbolic2ds1.h"

void BorderHyperbolic2DS1::main()
{

}

BorderHyperbolic2DS1::BorderHyperbolic2DS1()
{
}

double BorderHyperbolic2DS1::fi1(unsigned int i, unsigned int j) const
{
    return u(i,j, 0);
}

double BorderHyperbolic2DS1::fi2(unsigned int i, unsigned int j) const
{
    return 1.0;
}

double BorderHyperbolic2DS1::m1(unsigned int j, unsigned int k) const
{
    return u(0, j, k);
}

double BorderHyperbolic2DS1::m2(unsigned int j, unsigned int k) const
{
    return u(N1, j, k);
}

double BorderHyperbolic2DS1::m3(unsigned int i, unsigned int k) const
{
    return u(i, 0, k);
}

double BorderHyperbolic2DS1::m4(unsigned int i, unsigned int k) const
{
    return u(i, N2, k);
}

double BorderHyperbolic2DS1::f(unsigned int i, unsigned int j, unsigned int k) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    double t  = k*ht;
    return exp(t) + sin(x1) + cos(x2) - 2.0;
}

double BorderHyperbolic2DS1::u(unsigned int i, unsigned int j, unsigned int k) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    double t  = k*ht;
    return x1*x2 + x1*x2 + sin(x1) + cos(x2) + exp(t);
}
