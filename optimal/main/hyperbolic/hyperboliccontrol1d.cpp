#include "hyperboliccontrol1d.h"

void HyperbolicControl1D::main()
{
    DoubleVector u;
    HyperbolicControl1D hc;
    hc.calculate(u);
    Printer::printVector(u);
}

HyperbolicControl1D::HyperbolicControl1D()
    : RnFunction(), HyperbolicEquation(1000, 1000), Printer()
{
    dt = 0.05;
}

HyperbolicControl1D::~HyperbolicControl1D()
{

}

double HyperbolicControl1D::fx(const DoubleVector& x)
{
    return 0.0;
}

void HyperbolicControl1D::gradient(const DoubleVector& v, DoubleVector& g, double gradient_step)
{

}

double HyperbolicControl1D::fi1(unsigned int i) const
{
    double x = i*hx;
    return x*x;
}

double HyperbolicControl1D::fi2(unsigned int i) const
{
    return 0.0;
}

double HyperbolicControl1D::m1(unsigned int j) const
{
    double t = j*ht;
    return t*t;
}

double HyperbolicControl1D::m2(unsigned int j) const
{
    double t = j*ht;
    return t*t + 1.0;
}

double HyperbolicControl1D::f(unsigned int i, unsigned int j) const
{
    return 0.0;
}

void HyperbolicControl1D::print(unsigned int i, const DoubleVector &u, const DoubleVector &g, double a, RnFunction *fn) const
{
}
