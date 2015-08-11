#include "heatcontrol.h"

HeatControl::HeatControl(double t0, double t1, double x0, double x1, double dt, double dx)
{
    this->t0 = t0;
    this->t1 = t1;
    this->x0 = x0;
    this->x1 = x1;
    this->dt = dt;
    this->dx = dx;

    n = (unsigned int)(ceil((x1-x0)/dx)) + 1;
    m = (unsigned int)(ceil((t1-t0)/dt)) + 1;
}

double HeatControl::fx(const DoubleVector &u)
{
    return 0.0;
}

void HeatControl::gradient(double gradient_step, const DoubleVector &u, DoubleVector &g)
{
}

double HeatControl::U(double x) const
{
    return x*x + 2.0*x + 1.0;
}

double HeatControl::F(double x, double t)
{
    return 2.0*t - 2.0;
}

double HeatControl::m1(double t)
{
    return t*t;
}

double HeatControl::m2(double t)
{
    return t*t + 3.0;
}

double HeatControl::fi(double x)
{
    return x*x + 2.0*x;
}

double HeatControl::fxt1(double x, double t)
{
    unsigned int j = (unsigned int)(ceil(t/dt));
    unsigned int i = (unsigned int)(ceil(x/dx));
    return mf[j*m+i];
}

void HeatControl::main()
{

}
