#include "heatcontrol1.h"

HeatControl1::HeatControl1() : RnFunction(), ParabolicEquation(0.0, 1.0, 0.0, 1.0, 1.0, 1000, 1000)
{
}

double HeatControl1::fx(const DoubleVector &x)
{
    DoubleVector u;
    ParabolicEquation::calculateU(u);
    return 0.0;
}

void HeatControl1::gradient(const DoubleVector &x, DoubleVector &g, double gradient_step)
{
    DoubleVector u;
    ParabolicEquation::calculateU(u);
}

double HeatControl1::fi(unsigned int i, double x) const
{
    return u(i*hx, t0);
}

double HeatControl1::m1(unsigned int j, double t) const
{
    return u(x0, j*ht);
}

double HeatControl1::m2(unsigned int j, double t) const
{
    return u(x1, j*ht);
}

double HeatControl1::f(unsigned int i, double x, unsigned int j, double t) const
{
    return 2.0*j*ht - 2.0*a;
}

double HeatControl1::u(double x, double t) const
{
    return x*x+t*t;
}

void HeatControl1::calculate(const DoubleVector &f)
{
}

